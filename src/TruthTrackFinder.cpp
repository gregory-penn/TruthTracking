#include "TruthTrackFinder.h"

#include <k4ActsTracking/KFRunner.hxx>
#include <k4ActsTracking/SourceLink.hxx>

#include <k4FWCore/GaudiChecks.h>

#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <edm4hep/MCParticle.h>
#include <edm4hep/SimTrackerHit.h>
#include <edm4hep/TrackState.h>
#include <edm4hep/TrackerHit.h>

#include <podio/ObjectID.h>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {
  /// The hits of one MCParticle: the particle as found in the input MCParticle collection, and one SeedHit (global
  /// position, transverse radius, ACTS source link) per hit.
  struct HitGroup {
    edm4hep::MCParticle                particle;
    std::vector<ACTSTracking::SeedHit> hits;
  };

  /// Transverse momentum [GeV] from the omega of an edm4hep track state and the field Bz [T] at the origin.
  /// omega = 1/R in 1/mm and pT = 0.299792458e-3 * Bz * R.
  double ptFromOmega(float omega, double bzTesla) {
    return omega == 0.f ? 0. : 0.299792458e-3 * std::abs(bzTesla) / std::abs(omega);
  }
}  // namespace

DECLARE_COMPONENT(TruthTrackFinder)

TruthTrackFinder::TruthTrackFinder(const std::string& name, ISvcLocator* svcLoc)
    : MultiTransformer(name, svcLoc,
                       {KeyValues("InputTrackerHitCollection", {"MergedTrackerHits"}),
                        KeyValues("InputTrackerHitRelationCollection", {"MergedTrackerHitsRelations"})},
                       {KeyValues("OutputTrackCollection", {"TruthTracks"}),
                        KeyValues("OutputTrackRelationCollection", {"TruthTrackRelations"})}) {}

StatusCode TruthTrackFinder::initialize() {
  m_actsGeoSvc = svcLoc()->service<IActsGeoSvc>("ActsGeoSvc");
  K4_GAUDI_CHECK(m_actsGeoSvc);

  if (m_extrapolateToCalo && m_actsGeoSvc->caloSurfaceGeoIds().empty()) {
    warning() << "ExtrapolateToCalo requested but ActsGeoSvc provides no calorimeter-face surfaces; "
                 "no AtCalorimeter track states will be produced."
              << endmsg;
  }
  if (m_addEndcapCaloState && !m_extrapolateToCalo) {
    error() << "AddEndcapCaloState requested but ExtrapolateToCalo is off" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_minHitsPerTrk < 3) {
    // The start parameters are estimated from three hits (innermost, middle, outermost)
    error() << "MinHitsPerTrack must be at least 3, got " << m_minHitsPerTrk.value() << endmsg;
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection> TruthTrackFinder::operator()(
    const edm4hep::TrackerHitPlaneCollection& trackerHits,
    const edm4hep::TrackerHitSimTrackerHitLinkCollection& hitLinks) const {
  ++m_nEvents;

  // Default-construct ACTS contexts, as the CKF (CKFTrackingAlg) does. An empty geometry context means the nominal
  // sensor positions, with no alignment corrections.
  const Acts::GeometryContext      geoCtx = Acts::GeometryContext::dangerouslyDefaultConstruct();
  const Acts::MagneticFieldContext magCtx{};
  const Acts::CalibrationContext   calCtx{};

  // Step 1: convert every tracker hit into an ACTS measurement and source link (the same function the CKF uses). For each converted hit, keep its SeedHit (position, radius, source link), keyed by the hit itself.
  // The key is the hit handle, not ACTSTracking::trackerHitKey (built from the ObjectID): collections made during the
  // job may all carry the same collection ID until they are written, so ObjectID-based keys can collide across the
  // VXD/IT/OT collections merged into MergedTrackerHits (Claude is concerned about this since we saw it in calorimeter)
  //  A hit that has no ACTS surface gets no entry.
  ACTSTracking::SourceLinkContainer                             sourceLinks;
  ACTSTracking::MeasurementContainer                            measurements;
  ACTSTracking::HitContainer                                    hitContainer;
  std::unordered_map<edm4hep::TrackerHit, ACTSTracking::SeedHit> seedHitOf;
  seedHitOf.reserve(trackerHits.size());

  ACTSTracking::prepareTrackerHits(
      *this, *m_actsGeoSvc, geoCtx, trackerHits, measurements, sourceLinks, hitContainer, /*numThreads=*/1,
      [&](const edm4hep::TrackerHitPlane& hit, const ACTSTracking::SourceLink& sl, const Acts::Vector3& globalPos,
          const Acts::Surface& /*surface*/, const Acts::SquareMatrix2& /*localCov*/) {
        seedHitOf.emplace(edm4hep::TrackerHit(hit),
                          ACTSTracking::SeedHit{globalPos, std::hypot(globalPos.x(), globalPos.y()), sl});
      });

  // Step 2: group the hits by MCParticle. Each link connects one reconstructed hit ("from") to the SimTrackerHit it
  // was made from ("to"); a SimTrackerHit records exactly one MCParticle, which is used as it is, so the output links
  // point into whichever collection holds it: MCParticles without overlay, the overlay's copy of the signal particles
  // (NewMCParticles) with overlay. With overlay, the BIB SimTrackerHits have no particle object behind their particle
  // relation (isAvailable() is false), so BIB hits are skipped here.
  // The particle's ObjectID is not used: with overlay, the copied signal particles carry collection ID ffffffff
  // ("untracked") during the job and get their real ID only when written. The groups are therefore keyed by the
  // particle handle (which compares the underlying object), and kept in the order their particle is first seen, so
  // the output order does not depend on hashing.
  std::vector<HitGroup>                                groups;
  std::unordered_map<edm4hep::MCParticle, std::size_t> groupIndexOf;
  std::size_t nHitsNotInFit(0), nSecondaryHits(0), nParticlesNotFound(0);

  for (const auto& link : hitLinks) {
    const auto seedHitIter = seedHitOf.find(link.getFrom());
    if (seedHitIter == seedHitOf.end()) {
      // Not in the input hit collection, or no ACTS surface for its cellID
      ++nHitsNotInFit;
      continue;
    }
    const edm4hep::SimTrackerHit simHit = link.getTo();
    if (m_skipSecondaryHits && simHit.isProducedBySecondary()) {
      ++nSecondaryHits;
      continue;
    }
    const edm4hep::MCParticle particle = simHit.getParticle();
    if (!particle.isAvailable()) {
      // No particle object behind the relation: the BIB particles with overlay
      ++nParticlesNotFound;
      continue;
    }

    const auto [groupIter, isNew] = groupIndexOf.try_emplace(particle, groups.size());
    if (isNew) {
      groups.push_back(HitGroup{particle, {}});
    }
    groups[groupIter->second].hits.push_back(seedHitIter->second);
  }

  debug() << fmt::format("{} tracker hits, {} with a source link; {} hit links: {} to a hit without a source link, {} "
                         "from secondaries skipped, {} with no particle object (BIB); {} MCParticles with hits",
                         trackerHits.size(), seedHitOf.size(), hitLinks.size(), nHitsNotInFit, nSecondaryHits,
                         nParticlesNotFound, groups.size())
          << endmsg;

  // Step 3: fit each group with the Kalman fitter. Unlike the CKF, no hit is rejected by a chi2 cut: every hit of the
  // group is given to the fit.
  Acts::MagneticFieldProvider::Cache magCache = m_actsGeoSvc->magneticField()->makeCache(magCtx);

  // Bz at the origin, only for the fitted-pT debug print
  double bzTesla(0.);
  if (const auto field = m_actsGeoSvc->magneticField()->getField(Acts::Vector3::Zero(), magCache); field.ok()) {
    bzTesla = field.value().z() / Acts::UnitConstants::T;
  }

  const ACTSTracking::KFRunner kfRunner(*m_actsGeoSvc, geoCtx, magCtx, calCtx, measurements, hitContainer,
                                        {.propagateBackward  = m_propagateBackward,
                                         .extrapolateToCalo  = m_extrapolateToCalo,
                                         .addEndcapCaloState = m_addEndcapCaloState});

  edm4hep::TrackCollection               tracks;
  edm4hep::TrackMCParticleLinkCollection trackLinks;
  std::size_t nNeutral(0), nLowPt(0), nTooFewHits(0), nSeedFailed(0), nFitFailed(0);

  for (auto& group : groups) {
    const edm4hep::MCParticle& particle = group.particle;
    const auto&                mom      = particle.getMomentum();
    const double               truePt   = std::hypot(mom.x, mom.y);

    if (particle.getCharge() == 0.f) {
      ++nNeutral;
      continue;
    }
    if (truePt < m_minTruePt) {
      ++nLowPt;
      continue;
    }
    if (group.hits.size() < m_minHitsPerTrk) {
      ++nTooFewHits;
      continue;
    }

    // Order by transverse radius, as k4ActsTracking's collectSeedHits does. The start parameters come from the
    // innermost, middle and outermost hit. If two hits of the particle are on the same sensor surface, the Kalman
    // fitter uses only the first one it is given, i.e. the one at smaller radius.
    std::sort(group.hits.begin(), group.hits.end(),
              [](const ACTSTracking::SeedHit& a, const ACTSTracking::SeedHit& b) { return a.r < b.r; });

    std::optional<Acts::BoundTrackParameters> startParams = ACTSTracking::estimateSeedParameters(
        *this, *m_actsGeoSvc, geoCtx, group.hits, hitContainer, magCache, m_initialTrackError_pos,
        m_initialTrackError_phi, m_initialTrackError_lambda, m_initialTrackError_relP, m_initialTrackError_time);
    if (!startParams) {
      ++nSeedFailed;
      debug() << fmt::format("MCParticle {} (PDG {}, pT {:.3f} GeV, {} hits): start-parameter estimation failed",
                             particle.getObjectID().index, particle.getPDG(), truePt, group.hits.size())
              << endmsg;
      continue;
    }

    std::vector<Acts::SourceLink> groupSourceLinks;
    groupSourceLinks.reserve(group.hits.size());
    for (const ACTSTracking::SeedHit& hit : group.hits) {
      groupSourceLinks.emplace_back(hit.sl);
    }

    std::optional<edm4hep::MutableTrack> track =
        kfRunner.fit(*this, groupSourceLinks, *startParams, magCache, &m_caloMonitor);
    if (!track) {
      ++nFitFailed;
      debug() << fmt::format("MCParticle {} (PDG {}, pT {:.3f} GeV, {} hits): Kalman fit failed",
                             particle.getObjectID().index, particle.getPDG(), truePt, group.hits.size())
              << endmsg;
      continue;
    }

    // The track is added to the collection first, so the link points to the track stored in the output
    tracks.push_back(*track);
    auto trackLink = trackLinks.create();
    trackLink.setFrom(*track);
    trackLink.setTo(particle);
    trackLink.setWeight(1.f);  // every hit of the track is from this particle, by construction

    if (msgLevel(MSG::DEBUG)) {
      double fittedPt(-1.);
      for (const auto& state : track->getTrackStates()) {
        if (state.location == edm4hep::TrackState::AtIP) {
          fittedPt = ptFromOmega(state.omega, bzTesla);
        }
      }
      debug() << fmt::format("MCParticle {} (PDG {}): true pT {:.3f} GeV, fitted pT {:.3f} GeV; {} hits given, {} on "
                             "the track; chi2/ndf {:.1f}/{}",
                             particle.getObjectID().index, particle.getPDG(), truePt, fittedPt, group.hits.size(),
                             track->getTrackerHits().size(), track->getChi2(), track->getNdf())
              << endmsg;
    }
  }

  debug() << fmt::format("{} tracks from {} MCParticles with hits (rejected: {} neutral, {} below MinTruePt, {} below "
                         "MinHitsPerTrack, {} start-parameter failures, {} fit failures)",
                         tracks.size(), groups.size(), nNeutral, nLowPt, nTooFewHits, nSeedFailed, nFitFailed)
          << endmsg;

  m_nLinks += hitLinks.size();
  m_nHitsNotInFit += nHitsNotInFit;
  m_nSecondaryHits += nSecondaryHits;
  m_nParticlesNotFound += nParticlesNotFound;
  m_nCandidates += groups.size();
  m_nNeutral += nNeutral;
  m_nLowPt += nLowPt;
  m_nTooFewHits += nTooFewHits;
  m_nSeedFailed += nSeedFailed;
  m_nFitFailed += nFitFailed;
  m_nTracks += tracks.size();

  return std::make_tuple(std::move(tracks), std::move(trackLinks));
}

StatusCode TruthTrackFinder::finalize() {
  info() << fmt::format("Summary over {} events: {} hit links ({} to a hit without a source link, {} from secondaries "
                        "skipped, {} with no particle object (BIB)); {} MCParticles with hits -> {} "
                        "tracks (rejected: {} neutral, {} below MinTruePt, {} below MinHitsPerTrack, {} start-parameter "
                        "failures, {} fit failures)",
                        m_nEvents.load(), m_nLinks.load(), m_nHitsNotInFit.load(), m_nSecondaryHits.load(),
                        m_nParticlesNotFound.load(), m_nCandidates.load(), m_nTracks.load(), m_nNeutral.load(),
                        m_nLowPt.load(), m_nTooFewHits.load(), m_nSeedFailed.load(), m_nFitFailed.load())
         << endmsg;
  if (m_extrapolateToCalo) {
    info() << m_caloMonitor.summary() << endmsg;
  }
  return StatusCode::SUCCESS;
}

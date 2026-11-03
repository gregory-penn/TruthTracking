#include "TruthTrackFinder.h"

#include <k4ActsTracking/KFRunner.hxx>
#include <k4ActsTracking/SourceLink.hxx>

#include <k4FWCore/GaudiChecks.h>

#include <GaudiKernel/GaudiException.h>

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
#include <mutex>
#include <optional>
#include <unordered_map>
#include <unordered_set>
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

  if (m_fitMode.value() == "CKF") {
    m_useCKF = true;
    // Built once with the same configuration fields as CKFTrackingAlg; the reference surface for the AtIP state is
    // the default perigee at the origin
    m_ckfRunner.emplace(*m_actsGeoSvc,
                        ACTSTracking::CKFRunner::Config{.chi2CutOff            = m_CKF_chi2CutOff,
                                                        .numMeasurementsCutOff = m_CKF_numMeasurementsCutOff,
                                                        .chi2CutOffOutlier     = m_CKF_chi2CutOffOutlier,
                                                        .propagateBackward     = m_propagateBackward,
                                                        .extrapolateToCalo     = m_extrapolateToCalo,
                                                        .addEndcapCaloState    = m_addEndcapCaloState,
                                                        .useBranchStopper      = m_useBranchStopper,
                                                        .bsMaxHoles            = m_bsMaxHoles,
                                                        .bsMaxOutliers         = m_bsMaxOutliers,
                                                        .bsMinMeasurements     = m_bsMinMeasurements,
                                                        .bsPtMin               = m_bsPtMin,
                                                        .bsPtMinMeasurements   = m_bsPtMinMeasurements});
  } else if (m_fitMode.value() == "KF") {
    m_useCKF = false;
  } else {
    error() << "FitMode must be \"CKF\" or \"KF\", got \"" << m_fitMode.value() << "\"" << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "FitMode " << m_fitMode.value() << endmsg;
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
  std::size_t nHitsNotInFit(0), nHitsNotInInput(0), nSecondaryHits(0), nParticlesNotFound(0);

  // Every input hit, to tell the two reasons a link's hit can lack a source link apart (see below)
  std::unordered_set<edm4hep::TrackerHit> inputHits;
  inputHits.reserve(trackerHits.size());
  for (const auto& hit : trackerHits) {
    inputHits.insert(edm4hep::TrackerHit(hit));
  }

  for (const auto& link : hitLinks) {
    const auto seedHitIter = seedHitOf.find(link.getFrom());
    if (seedHitIter == seedHitOf.end()) {
      if (inputHits.contains(link.getFrom())) {
        // In the input hit collection, but prepareTrackerHits refused it (no ACTS surface for its cellID, or more
        // than 0.5 um off its surface). prepareTrackerHits already printed a warning, and the CKF loses it too.
        ++nHitsNotInFit;
      } else {
        // Not in the input hit collection at all: the hit and relation collections do not belong together
        ++nHitsNotInInput;
      }
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

  // A link to a hit that is not in the input hit collection can only come from a configuration error (e.g. the
  // relations of the coned hits given together with the un-coned hits), so the job is stopped.
  if (nHitsNotInInput > 0) {
    throw GaudiException(fmt::format("{} of {} hit links point to a hit that is not in the input tracker-hit collection "
                                     "'{}': the hit collection and the relation collection '{}' do not belong together",
                                     nHitsNotInInput, hitLinks.size(), inputLocations(0).front(),
                                     inputLocations(1).front()),
                         name(), StatusCode::FAILURE);
  }

  // Step 3: turn each group into a track, with the CKF (FitMode "CKF") or the Kalman fitter (FitMode "KF")
  Acts::MagneticFieldProvider::Cache magCache = m_actsGeoSvc->magneticField()->makeCache(magCtx);

  // Bz at the origin, only for the fitted-pT debug print
  double bzTesla(0.);
  if (const auto field = m_actsGeoSvc->magneticField()->getField(Acts::Vector3::Zero(), magCache); field.ok()) {
    bzTesla = field.value().z() / Acts::UnitConstants::T;
  }

  // The Kalman fitter binds this event's measurements, so it is built per event (only for FitMode "KF")
  std::optional<ACTSTracking::KFRunner> kfRunner;
  if (!m_useCKF) {
    kfRunner.emplace(*m_actsGeoSvc, geoCtx, magCtx, calCtx, measurements, hitContainer,
                     ACTSTracking::KFRunner::Config{.propagateBackward  = m_propagateBackward,
                                                   .extrapolateToCalo  = m_extrapolateToCalo,
                                                   .addEndcapCaloState = m_addEndcapCaloState});
  }

  edm4hep::TrackCollection               tracks;
  edm4hep::TrackMCParticleLinkCollection trackLinks;
  std::size_t nNeutral(0), nLowPt(0), nTooFewHits(0), nSeedFailed(0), nFitFailed(0), nCkfCandidates(0);

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
    // innermost, middle and outermost hit; both fits start on the innermost hit's surface. FitMode "KF": if two hits
    // of the particle are on the same sensor surface, the Kalman fitter uses only the first one it is given, i.e. the
    // one at smaller radius. FitMode "CKF": the order does not matter, the CKF looks the hits up by surface.
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

    std::optional<edm4hep::MutableTrack> track;
    std::size_t                          nCandidatesThisParticle(0);
    if (m_useCKF) {
      // The CKF looks up the hits on each surface it crosses in a source-link container (sorted by surface). It is
      // given a container with only this particle's hits, so it can only choose among them, and with the one seed
      // from the start parameters. The measurements and hitContainer are the whole event's; the source links index
      // into them.
      ACTSTracking::SourceLinkContainer particleSourceLinks;
      particleSourceLinks.reserve(group.hits.size());
      for (const ACTSTracking::SeedHit& hit : group.hits) {
        particleSourceLinks.insert(hit.sl);
      }

      // The CKF can return several candidates for one seed (it branches when more than one of the particle's hits on
      // a surface passes the chi2 cut). They are collected in a scratch collection; the one with the most hits, then
      // the lowest chi2, is kept (the ordering ACTSDuplicateRemoval uses) and copied into the output.
      edm4hep::TrackCollection candidates;
      std::mutex               candidatesMutex;  // findTracks locks it to add to the collection
      m_ckfRunner->findTracks(*this, measurements, particleSourceLinks, hitContainer,
                              std::vector<Acts::BoundTrackParameters>{*startParams}, magCache, candidates,
                              candidatesMutex, &m_caloMonitor);
      nCandidatesThisParticle = candidates.size();
      nCkfCandidates += candidates.size();
      if (!candidates.empty()) {
        const auto best = std::min_element(candidates.begin(), candidates.end(),
                                           [](const edm4hep::Track& a, const edm4hep::Track& b) {
                                             if (a.getTrackerHits().size() != b.getTrackerHits().size()) {
                                               return a.getTrackerHits().size() > b.getTrackerHits().size();
                                             }
                                             return a.getChi2() < b.getChi2();
                                           });
        // clone() copies the track's data, hits and track states into a new track that belongs to no collection yet
        track = (*best).clone();
      }
    } else {
      std::vector<Acts::SourceLink> groupSourceLinks;
      groupSourceLinks.reserve(group.hits.size());
      for (const ACTSTracking::SeedHit& hit : group.hits) {
        groupSourceLinks.emplace_back(hit.sl);
      }
      track = kfRunner->fit(*this, groupSourceLinks, *startParams, magCache, &m_caloMonitor);
    }

    if (!track) {
      ++nFitFailed;
      debug() << fmt::format("MCParticle {} (PDG {}, pT {:.3f} GeV, {} hits): {}", particle.getObjectID().index,
                             particle.getPDG(), truePt, group.hits.size(),
                             m_useCKF ? "CKF found no track" : "Kalman fit failed")
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
                             "the track, {} holes; chi2/ndf {:.1f}/{}{}",
                             particle.getObjectID().index, particle.getPDG(), truePt, fittedPt, group.hits.size(),
                             track->getTrackerHits().size(), track->getNholes(), track->getChi2(), track->getNdf(),
                             m_useCKF ? fmt::format("; {} CKF candidates", nCandidatesThisParticle) : "")
              << endmsg;
    }
  }

  debug() << fmt::format("{} tracks from {} MCParticles with hits (rejected: {} neutral, {} below MinTruePt, {} below "
                         "MinHitsPerTrack, {} start-parameter failures, {} fit failures); {} CKF candidates",
                         tracks.size(), groups.size(), nNeutral, nLowPt, nTooFewHits, nSeedFailed, nFitFailed,
                         nCkfCandidates)
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
  m_nCkfCandidates += nCkfCandidates;
  m_nTracks += tracks.size();

  return std::make_tuple(std::move(tracks), std::move(trackLinks));
}

StatusCode TruthTrackFinder::finalize() {
  info() << fmt::format("Summary over {} events (FitMode {}): {} hit links ({} to a hit without a source link, {} from "
                        "secondaries skipped, {} with no particle object (BIB)); {} MCParticles with hits -> {} "
                        "tracks (rejected: {} neutral, {} below MinTruePt, {} below MinHitsPerTrack, {} start-parameter "
                        "failures, {} fit failures); {} CKF candidates",
                        m_nEvents.load(), m_fitMode.value(), m_nLinks.load(), m_nHitsNotInFit.load(),
                        m_nSecondaryHits.load(), m_nParticlesNotFound.load(), m_nCandidates.load(), m_nTracks.load(),
                        m_nNeutral.load(), m_nLowPt.load(), m_nTooFewHits.load(), m_nSeedFailed.load(),
                        m_nFitFailed.load(), m_nCkfCandidates.load())
         << endmsg;
  if (m_extrapolateToCalo) {
    info() << m_caloMonitor.summary() << endmsg;
  }
  return StatusCode::SUCCESS;
}

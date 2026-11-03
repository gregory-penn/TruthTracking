#pragma once

// Truth-based track finding: the tracker hits are grouped by the MCParticle that made them (instead of by a pattern
// recognition such as the CKF), and each group is fitted with the ACTS Kalman fitter of k4ActsTracking
// (ACTSTracking::KFRunner). KFRunner uses the same propagator, material, field, measurement calibration, Kalman update
// and edm4hep conversion as the CKF (ACTSTracking::CKFRunner); the difference is that the CKF chooses, surface by
// surface, which hit to add (and rejects hits with chi2 above its cut), while KFRunner fits every hit it is given. The
// fitted momentum therefore still includes the hit resolution, the material and the field, but no wrong, missing,
// duplicate or fake tracks.

#include <k4FWCore/Transformer.h>

#include <k4ActsTracking/IActsGeoSvc.h>
#include <k4ActsTracking/RunnerCommon.hxx>

#include <Acts/Definitions/Units.hpp>

#include <Gaudi/Property.h>
#include <GaudiKernel/SmartIF.h>

#include <edm4hep/MCParticleCollection.h>
#include <edm4hep/TrackCollection.h>
#include <edm4hep/TrackMCParticleLinkCollection.h>
#include <edm4hep/TrackerHitPlaneCollection.h>
#include <edm4hep/TrackerHitSimTrackerHitLinkCollection.h>

#include <atomic>
#include <cstddef>
#include <string>
#include <tuple>

struct TruthTrackFinder final
    : k4FWCore::MultiTransformer<std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection>(
          const edm4hep::TrackerHitPlaneCollection&, const edm4hep::TrackerHitSimTrackerHitLinkCollection&,
          const edm4hep::MCParticleCollection&)> {
  TruthTrackFinder(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize() override;

  /// Called once per event by the framework with the three input collections; the two returned collections are put
  /// into the event store under the output names.
  ///
  /// @param trackerHits  The digitized tracker hits (MergedTrackerHits: VXD, IT and OT hits in one collection).
  /// @param hitLinks     One link per (digitized hit, simulated hit) pair: getFrom() is the digitized hit
  ///                     (edm4hep::TrackerHit, i.e. an element of trackerHits), getTo() is the edm4hep::SimTrackerHit
  ///                     it was made from. The SimTrackerHit gives the MCParticle (getParticle()).
  /// @param mcParticles  The MCParticles that the SimTrackerHits point to. A hit is used only if its particle is in
  ///                     this collection.
  /// @return             The fitted tracks, and one link per track: getFrom() is the track, getTo() the MCParticle.
  std::tuple<edm4hep::TrackCollection, edm4hep::TrackMCParticleLinkCollection> operator()(
      const edm4hep::TrackerHitPlaneCollection& trackerHits,
      const edm4hep::TrackerHitSimTrackerHitLinkCollection& hitLinks,
      const edm4hep::MCParticleCollection& mcParticles) const override;

  /// @name Selection of the particles that get a track
  ///@{
  Gaudi::Property<std::size_t> m_minHitsPerTrk{this, "MinHitsPerTrack", 3,
                                               "Minimum number of hits of one MCParticle for it to be fitted"};
  Gaudi::Property<double>      m_minTruePt{this, "MinTruePt", 0.,
                                      "Minimum true transverse momentum [GeV] of an MCParticle for it to be fitted"};
  Gaudi::Property<bool>        m_skipSecondaryHits{
      this, "SkipSecondaryHits", true,
      "Skip hits whose SimTrackerHit has isProducedBySecondary() set. Geant4 records a hit made by a particle that "
             "is not stored in the MCParticles (e.g. a delta-ray electron) under its nearest stored ancestor, so such a "
             "hit lies near, but not on, the ancestor's path."};
  ///@}

  /// @name Kalman-fit configuration, the same properties and defaults as the CKF (CKFTrackingAlg)
  ///@{
  Gaudi::Property<bool>   m_propagateBackward{this, "PropagateBackward", false, "Extrapolates tracks towards beamline."};
  Gaudi::Property<bool>   m_extrapolateToCalo{
      this, "ExtrapolateToCalo", true,
      "Extrapolate fitted tracks to the calorimeter face and add an AtCalorimeter track state."};
  Gaudi::Property<bool>   m_addEndcapCaloState{
      this, "AddEndcapCaloState", false,
      "Give a track that crosses both calorimeter sections one AtCalorimeter track state per section instead of a "
        "single one. Ignored unless ExtrapolateToCalo is set."};
  Gaudi::Property<double> m_initialTrackError_pos{this, "InitialTrackError_Pos", 10 * Acts::UnitConstants::um,
                                                  "Initial track error for local position."};
  Gaudi::Property<double> m_initialTrackError_phi{this, "InitialTrackError_Phi", 1 * Acts::UnitConstants::degree,
                                                  "Initial track error for phi."};
  Gaudi::Property<double> m_initialTrackError_relP{this, "InitialTrackError_RelP", 0.25,
                                                   "Initial track error for momentum (relative)."};
  Gaudi::Property<double> m_initialTrackError_lambda{this, "InitialTrackError_Lambda", 1 * Acts::UnitConstants::degree,
                                                     "Initial track error for lambda."};
  Gaudi::Property<double> m_initialTrackError_time{this, "InitialTrackError_Time", 100 * Acts::UnitConstants::ns,
                                                   "Initial track error for time."};
  ///@}

private:
  SmartIF<IActsGeoSvc> m_actsGeoSvc{nullptr};

  /// Calorimeter-face extrapolation monitoring, updated by the per-event KFRunner and summarized in finalize().
  ACTSTracking::CaloExtrapMonitor m_caloMonitor{};

  /// Counters summed over all events and printed in finalize(). operator() is const, hence mutable (and atomic,
  /// since Gaudi may run several events at once).
  mutable std::atomic<std::size_t> m_nEvents{0};
  mutable std::atomic<std::size_t> m_nLinks{0};                ///< hit links read
  mutable std::atomic<std::size_t> m_nHitsNotInFit{0};         ///< link to a hit without a source link
  mutable std::atomic<std::size_t> m_nSecondaryHits{0};        ///< hit skipped by SkipSecondaryHits
  mutable std::atomic<std::size_t> m_nParticlesNotFound{0};    ///< hit from a particle not in the input (e.g. BIB)
  mutable std::atomic<std::size_t> m_nCandidates{0};           ///< MCParticles with at least one hit
  mutable std::atomic<std::size_t> m_nNeutral{0};              ///< candidates rejected as neutral
  mutable std::atomic<std::size_t> m_nLowPt{0};                ///< candidates rejected by MinTruePt
  mutable std::atomic<std::size_t> m_nTooFewHits{0};           ///< candidates rejected by MinHitsPerTrack
  mutable std::atomic<std::size_t> m_nSeedFailed{0};           ///< start-parameter estimation failed
  mutable std::atomic<std::size_t> m_nFitFailed{0};            ///< Kalman fit failed
  mutable std::atomic<std::size_t> m_nTracks{0};               ///< tracks written
};

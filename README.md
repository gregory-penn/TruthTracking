# TruthTracking
To perform truth tracking in a muon collider environment.

`TruthTrackFinder` groups the tracker hits by the MCParticle that made them (through the hit → SimTrackerHit links)
and turns each group into a track in one of two ways (property `FitMode`):
- `CKF` (default): the Combinatorial Kalman Filter of k4ActsTracking (`ACTSTracking::CKFRunner`), given only the
  particle's hits, so it applies the same chi2 hit rejection as the standard CKF.
- `KF`: the Kalman fitter of k4ActsTracking (`ACTSTracking::KFRunner::fit`), which fits every hit.

Both use the same propagator, material, field, Kalman update and edm4hep conversion as the standard CKF. The purpose of this is to isolate detector resolutions from reconstruction algorithm performance.

Inputs: `MergedTrackerHits`, `MergedTrackerHitsRelations`. Each hit's MCParticle is read from its SimTrackerHit, so
the output links point to `MCParticles` without overlay and to the overlay's copy (e.g. `NewMCParticles`) with
overlay; with overlay the BIB hits have no particle object behind them and are skipped. Outputs: `TruthTracks`, and `TruthTrackRelations` (one `TrackMCParticleLink` per track,
weight 1).

## Build (inside the mucoll-sim image)
```bash
source /opt/setup_mucoll.sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$PWD/install
cmake --build build -j8
cmake --install build
```
It builds using the k4ActsTracking installed in the image; nothing else needs rebuilding.

## Run (MAIAConfig)
```bash
export LD_LIBRARY_PATH=/path/to/TruthTracking/install/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/path/to/TruthTracking/install/python:$PYTHONPATH
k4run .../digi_reco_steer.py ... --doTruthTracking [--truthTrackingFit CKF|KF|both] [--TruthTrackFinder.OutputLevel 2]
```
`--truthTrackingFit` selects the instances: `CKF` (default) runs `TruthTrackFinder` into `TruthTracks` /
`TruthTrackRelations`; `KF` runs `TruthTrackFinderKF` into `TruthTracksKF` / `TruthTrackRelationsKF`; `both` runs the
two in the same job.
`--doTruthTracking` runs the finder alongside the standard CKF chain (configured in `Tracking/truth_tracking.py`);
`SiTracks` is unchanged. Output level 2 (DEBUG) prints, for each track, the true and fitted pT, the number of hits and
the chi2.

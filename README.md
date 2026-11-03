# TruthTracking
To perform truth tracking in a muon collider environment.

`TruthTrackFinder` groups the tracker hits by the MCParticle that made them (through the hit → SimTrackerHit links)
and fits each group with the ACTS Kalman fitter of k4ActsTracking (`ACTSTracking::KFRunner::fit`). It uses the same
propagator, material, field, Kalman update and edm4hep conversion as the CKF; unlike the CKF, it does not reject hits
by a chi2 cut. The purpose of this is to isolate detector resolutions from reconstruction algorithm performance.

Inputs: `MergedTrackerHits`, `MergedTrackerHitsRelations`, an MCParticle collection (`MCParticles`, or the overlay's
copy, e.g. `NewMCParticles`). Outputs: `TruthTracks`, and `TruthTrackRelations` (one `TrackMCParticleLink` per track,
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
k4run .../digi_reco_steer.py ... --doTruthTracking [--TruthTrackFinder.OutputLevel 2]
```
`--doTruthTracking` runs the finder alongside the standard CKF chain (configured in `Tracking/truth_tracking.py`);
`SiTracks` is unchanged. Output level 2 (DEBUG) prints, for each track, the true and fitted pT, the number of hits and
the chi2.

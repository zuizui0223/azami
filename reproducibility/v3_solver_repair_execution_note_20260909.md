# Chapter 1 v3 solver repair execution note

This note records that the optimizer repair is being validated before any final calibration launch.

- The pinned heterogeneous-slopes failure was reproduced on the same seed and same draw.
- The visible-colour fit reached nearly identical REML solutions from all four starts but SciPy L-BFGS-B returned an abnormal termination flag.
- The repair does not accept the abnormal flag directly. It passes the candidate through the pre-existing score-root stationarity verification and requires stationarity, positive curvature, and non-worsening objective before acceptance.
- No calibration threshold, seed, synthetic truth, spatial grid, or family definition was changed.
- Final calibration remains blocked until repaired smoke, repaired full-family mechanics, and current-head CI all pass.

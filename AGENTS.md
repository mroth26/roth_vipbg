# AGENTS.md

## Project overview

This repository (`twin_sim/`) contains R scripts for a behavior-genetics twin
simulation study. It uses the **OpenMx** structural-equation-modeling package to
simulate MZ/DZ twin data and fit ACE/ACDE + sibling-interaction (univariate) and
ACE + direction-of-causation (bivariate) models, then performs model selection
and AIC model averaging.

Files:
- `twin_sim/univariate_sim.R`, `twin_sim/bivariate_sim.R` — data generators
  (`simulate_twins`). The runners `source()` these from a `models/` subdirectory.
- `twin_sim/univariate_runner.R`, `twin_sim/bivariate_runner.R` — the
  "applications": define the fitting/selection/averaging functions, then run a
  1000-replication Monte-Carlo loop at the bottom that writes CSVs per run.

There is no build system, test suite, linter config, or package manifest — these
are plain `Rscript` programs.

## Cursor Cloud specific instructions

### Required packages
- R 4.3.x with `mvtnorm`, `psych`, and **OpenMx**.
- **OpenMx MUST be the NPSOL-enabled build**, not the CRAN build. Both runners
  hard-code `mxOption(NULL, "Default optimizer", "NPSOL")`, which errors unless
  NPSOL is compiled in. Verify with:
  `Rscript -e 'library(OpenMx); print(mxAvailableOptimizers())'` — the list must
  contain `NPSOL`.

### Rebuilding the NPSOL-enabled OpenMx (only if it is missing)
The startup update script does NOT rebuild OpenMx (it is a source build and is
preserved in the VM snapshot). If `NPSOL` is absent, rebuild once like this — do
NOT just `install.packages("OpenMx")` (that gives the non-NPSOL CRAN build):
1. `mkdir -p /tmp/omxsrc && Rscript -e 'download.packages("OpenMx", destdir="/tmp/omxsrc", repos="https://vipbg.vcu.edu/vipbg/OpenMx2/software/", type="source")'`
2. `cd /tmp/omxsrc && tar xzf OpenMx_*.tar.gz`
3. The repo's `configure` only ships Linux `libnpsol.a` checksums up to gcc12.2,
   while this image has gcc 13.x, so the auto-download fetches an HTML 404 page
   ("file format not recognized" at link time). Pre-place the gcc12.2 lib (same
   ABI, links fine) so `configure` skips the download:
   `curl -fsSL -o /tmp/omxsrc/OpenMx/src/libnpsol.a https://vipbg.vcu.edu/vipbg/OpenMx2/software/packages/npsol/linux/x86_64/a21b02cfe4d55dacd70408182becd0170b060423`
   (sha1 `a21b02cfe4d55dacd70408182becd0170b060423`)
4. `cd /tmp/omxsrc && R CMD INSTALL --no-multiarch OpenMx` (compiles in ~8 min).

### Running the applications
- The runners expect to be launched from a working directory that contains a
  `models/` subfolder with the matching `*_sim.R` file (they call
  `source('models/univariate_sim.R')` / `source('models/bivariate_sim.R')`).
- Running a runner file top-to-bottom executes a **1000-replication loop** —
  impractical for a quick check. For a smoke test, load only the function
  definitions (everything above the `n_reps <- 1000` line) and call
  `run_simulation_and_analysis(seed = 1L)` once. Univariate takes ~2–3 s,
  bivariate ~20 s (it fits ~25 models per replication via `mxTryHard`).
- The bivariate `run_simulation_and_analysis` reads `output_prefix` from the
  enclosing (global) scope rather than from its formals; define
  `output_prefix <- ""` in the calling environment before invoking it (the
  univariate version defaults `output_prefix = ""` internally and needs nothing).
- Both runners emit many OpenMx optimizer warnings during fitting; these are
  normal for this simulation and do not indicate failure.

### Lint / test / build
- No linter, automated test suite, or build step exists. "Running" means
  executing the simulation pipeline as above.

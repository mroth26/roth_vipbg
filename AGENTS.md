# AGENTS.md

## Cursor Cloud specific instructions

### Product overview

This repository is an **R batch research pipeline** for behavioral genetics twin-study Monte Carlo simulations. There are no web servers, databases, or Docker services. End-to-end testing means running the R scripts in `twin_sim/`.

| Pipeline | Entry point | Purpose |
|---|---|---|
| Univariate | `twin_sim/univariate_runner.R` | ACE/ACES twin models with three inference methods |
| Bivariate | `twin_sim/bivariate_runner.R` | Bivariate ACE + DoC models (~25 variants) |

### Required dependencies

- **R** (4.3+): `r-base`, `r-base-dev`
- **OpenMx**: install the Ubuntu binary `r-cran-openmx` (CRAN source builds fail on R 4.3 due to `Rf_isDataFrame`; the VCU NPSOL source build also fails to link `libnpsol.a` on Linux x86_64)
- **CRAN packages** (user library): `psych`, `mvtnorm` — install to `~/R/library` because `/usr/local/lib/R/site-library` is not writable

Runners `source('models/univariate_sim.R')` and `source('models/bivariate_sim.R')`, but sim files live at `twin_sim/*.R`. The update script creates `twin_sim/models/` symlinks; run all scripts with `cd twin_sim`.

### NPSOL optimizer caveat

Both runners call `mxOption(NULL, "Default optimizer", "NPSOL")`. The Ubuntu `r-cran-openmx` package is **SLSQP-only** (no NPSOL on Linux). For smoke tests, substitute `SLSQP` when sourcing, or change the optimizer line locally. Full production runs that require NPSOL need a macOS/Windows OpenMx build or a working NPSOL-linked Linux build.

### Commands

See `twin_sim/` scripts for the application itself. Typical workflow:

```bash
cd /workspace/twin_sim

# Quick module test (no OpenMx fitting)
Rscript -e 'source("models/univariate_sim.R"); str(simulate_twins(n_MZ=50, n_DZ=50, seed=1))'

# Single univariate replication (skip the 1000-rep loop at file bottom)
Rscript -e '
  lines <- readLines("univariate_runner.R")
  cut <- which(grepl("^n_reps <-", lines))[1] - 1
  lines <- gsub("NPSOL", "SLSQP", lines)
  eval(parse(text=lines[1:cut]), envir=.GlobalEnv)
  dir.create("smoke_test", showWarnings=FALSE)
  run_simulation_and_analysis(seed=1L, output_prefix="smoke_test/run_0001_")
'
```

**Lint / test:** not configured in-repo. No `testthat` suite. Use the smoke commands above.

**Full production run:** `Rscript univariate_runner.R` from `twin_sim/` runs **1000 replications** and takes hours. Set `n_reps <- 1` in the runner for local debugging.

### Known repo issues (not env-related)

- `bivariate_runner.R`: `run_simulation_and_analysis()` references `output_prefix` but it is not a function parameter.
- Bottom replication loops expect return fields (`res$data`, `res$averaged_estimates`) that differ from what the function returns (`sim_data`, `method2_metrics`, etc.).

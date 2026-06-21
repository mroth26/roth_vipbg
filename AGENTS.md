# AGENTS.md

## Cursor Cloud specific instructions

### What this repo is
A pure-**R** behavioral-genetics twin-simulation toolkit in `twin_sim/`. There is **no service, server, database, or GUI** — each "product" is a batch script run with `Rscript`. Two products:

- **Univariate** ACE/ACDE + sibling-interaction simulation: `twin_sim/univariate_runner.R` (sources `twin_sim/univariate_sim.R`).
- **Bivariate** Direction-of-Causation + sibling-interaction simulation: `twin_sim/bivariate_runner.R` (sources `twin_sim/bivariate_sim.R`).

Each runner: simulates MZ/DZ twin data → fits OpenMx structural-equation models → does model selection + Akaike multimodel averaging → writes per-replication CSVs. There is **no lint/test/build tooling** (no `testthat`, `lintr`, `Makefile`, CI). "Running" the code is the only verification.

### How to run (non-obvious working-directory + loop caveats)
The runners hardcode `source('models/<name>_sim.R')` **relative to the current working directory**, and write outputs into the CWD. The sim files actually live flat in `twin_sim/`, not in a `models/` subdir. So run from a directory that contains a `models/` folder pointing at the sim file:

```bash
mkdir -p /tmp/uni_run/models
ln -s /workspace/twin_sim/univariate_sim.R /tmp/uni_run/models/univariate_sim.R
cd /tmp/uni_run && Rscript /workspace/twin_sim/univariate_runner.R
```

Each runner loops `n_reps <- 1000` replications in an infinite `repeat{}` retry loop (it advances the seed and retries forever on any per-replication error). For a quick smoke check, let the first few replications finish (each writes `simulation_results_*/run_NNNN_*.csv`) and then stop the process by **PID** (never `pkill`). A single univariate replication takes only ~2–3 s.

### Optimizer / OpenMx + NPSOL gotcha (important)
The runners call `mxOption(NULL, "Default optimizer", "NPSOL")`. **NPSOL is proprietary and is NOT in the CRAN build of OpenMx.** The environment ships an NPSOL-enabled OpenMx build; verify with:

```r
library(OpenMx); imxHasNPSOL()   # must be TRUE
```

If OpenMx/NPSOL ever needs to be rebuilt (e.g. a fresh VM where it is missing), note the official installer `source('https://vipbg.vcu.edu/vipbg/OpenMx2/software/getOpenMx.R')` is **broken on modern gcc**: OpenMx's `configure` builds the path `linux/x86_64/gcc<full-version>` (e.g. `gcc13.3.0`), which has no entry in `inst/npsol/checksum` (it stops at gcc12.2), so it downloads the directory index page as `src/libnpsol.a` and the final link fails with `file format not recognized`. Recovery (what this env uses): build OpenMx **2.21.13** source from the VCU repo, but first drop a real NPSOL archive at `src/libnpsol.a` (the gcc8.3–gcc12.2 archive `a21b02cfe4d55dacd70408182becd0170b060423`, which bundles its own BLAS and is ABI-compatible with gcc13), then `R CMD INSTALL`. R must be linked against reference BLAS (it is by default here); NPSOL misbehaves with OpenBLAS/MKL.

### Known pre-existing bug (NOT an environment issue)
`twin_sim/bivariate_runner.R` fails on **every** replication with `object 'output_prefix' not found`. Its `run_simulation_and_analysis()` signature (around line 658) is missing the `output_prefix = ""` parameter that the body uses (the univariate version has it). Because of the infinite retry loop, the bivariate runner spins forever without producing output. This is a bug in the repo code, not in the environment — the bivariate data generator (`bivariate_sim.R::simulate_twins()`) and OpenMx fitting both work fine. Do not "fix" it unless asked.

# TIEE

**Tail-Informed Extreme Quantile Treatment Effects**

R code accompanying the paper:

> Li, M. & Castro-Camilo, D. (2026). Tail-Informed Extreme Quantile Treatment Effects. arXiv: [2603.23309](https://arxiv.org/abs/2603.23309), under review.

## Repository structure

```
functions/              Core estimation and inference functions
  ├── gpd_cpp.R                  TIEE estimator with influence-function SEs
  ├── tools.R                    Rcpp-compiled C++ utilities (GPD score, density, quantile, KDE path)
  ├── calculate_true_qtes.R      Monte Carlo true QTE computation
  ├── h_fun.R                    Heavy-tailed DGP scenarios (S1–S3)
  ├── h_light.R                  Light-tailed DGP scenarios (S4–S6)
  ├── qte_extrapolation_hill.R   Causal Hill estimator
  ├── qte_extrapolation_pickands Causal Pickands estimator
  ├── qte_extrapolation_bootstrap
  └── qte_firpo_zhang.R          Firpo–Zhang QTE estimator

simulation/             Monte Carlo simulation studies
  ├── simulation.R               Single-scenario driver (edit config block at top)
  ├── run_simu.R                 Batch driver with config list
  ├── analyze.R                  Post-processing: metrics tables + coverage plots
  │                               (replaces the former heavy*/light*/single/summarize/make_table scripts)
  ├── missing.R                  Propensity score misspecification study (Table 3)
  ├── sensitive_u.R / sensitive_u_plot.R   Threshold u sensitivity (Appendix B.2)
  └── sensitive_k.R / k_plot.R             Grid size K sensitivity (Appendix B.3)

application/            River precipitation case study
  ├── application.R              Main analysis: TIEE, Hill, Zhang + LaTeX table
  ├── clim_idx.R                 Download and merge climate indices (NAO, AO, AMO, ENSO, PDO)
  └── data/
      ├── eqte.Rdata             Merged precipitation + covariate dataset
      ├── AMO.txt / ONI.txt / PDO.txt         Monthly climate indices
      ├── norm.daily.*.csv                     Daily NAO/AO indices
      └── SFC_Details.txt                      Station metadata
```

> **Note:** `rr_ens_spread_0.1deg_reg_v31.0e.nc` (1.5 GB) is excluded from this repository due to size.

## How to run

### Simulation (Section 5)

1. Edit the configuration block at the top of `simulation/simulation.R` (DGP, sample size, cores, etc.).
2. Run:
   ```r
   source("simulation/simulation.R")
   ```
   Results are saved under `simulation/scratch/`.

3. Post-process all scenarios (heavy/light × n=1000/5000) in one step:
   ```r
   source("simulation/analyze.R")
   ```
   This produces LaTeX tables (`*_table.tex`) and coverage plots (`*_coverage.pdf`).

### Application (Section 6)

1. Preprocess climate indices (only needed once):
   ```r
   source("application/clim_idx.R")
   ```

2. Run the main analysis:
   ```r
   source("application/application.R")
   ```

## Dependencies

- R (>= 4.0)
- Packages: `Rcpp`, `evd`, `actuar`, `quantreg`, `ismev`, `mev`, `parallel`, `BMisc`, `dplyr`, `tidyr`, `purrr`, `ggplot2`, `kableExtra`, `lubridate`, `readr`, `extRemes`
- A C++ compiler compatible with Rcpp (the core computation functions in `tools.R` are implemented as inline C++ via `cppFunction()`)

## Citation

If you use this code, please cite:

```bibtex
@article{LiCastroCamilo2026,
  title   = {Tail-Informed Extreme Quantile Treatment Effects},
  author  = {Li, Mengran and Castro-Camilo, Daniela},
  journal = {arXiv preprint arXiv:2603.23309},
  year    = {2026}
}
```

## License

Code released under the MIT License. Please cite the paper if you use this code in your research.

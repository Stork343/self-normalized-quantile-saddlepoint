# Self-Normalized Quantile Empirical Saddlepoint Approximation (SNQESA)

R code for the SNQESA method: a density-free approach to accurate quantile
inference without bandwidth selection.

Accompanying paper:

> Hou, J., Meng, T., & Tian, M. (2026). Self-Normalized Quantile Empirical
> Saddlepoint Approximation. *Statistical Papers*, 67(5), 106.
> https://doi.org/10.1007/s00362-026-01882-3

## Repository structure

- `main.R` — core SNQESA implementation (distribution-normalized quantile
  saddlepoint functions).
- `simulation.R` — simulation studies comparing SNQESA with kernel,
  bootstrap, and empirical-likelihood baselines.
- `real_data.R` — real-data application.
- `plot_sim.R` — figures for the simulation results.
- `plot_real.R` — figures for the real-data results.

## Requirements

- R >= 4.0
- Optional packages used in the scripts are listed in the file headers
  (e.g., `pbmcapply`, `dplyr`, `knitr`, `kableExtra`).

## Usage

Each script is self-contained. Source the file of interest from R, for
example:

```r
source("simulation.R")
```

## License

MIT (c) 2025 Hou Jian.
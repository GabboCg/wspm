# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

R implementation of the noise-robust intraday jump test from Christensen, Timmermann & Veliyev (2025), *Journal of Financial Economics* 167, 104010. The methodology uses pre-averaging and subsampling to detect price jumps in tick-by-tick data affected by microstructure noise.

## Running the code

Requires R with `RcppArmadillo` installed. Load everything in one step from the project root (WD is set automatically by `wspm.Rproj` in RStudio):

```r
source("load.R")
```

This compiles the C++ via `Rcpp::sourceCpp("src/subsampler.cpp")` and sources `R/jump_test.R` which defines `intradayJumpTest()`.

Typical call with the bundled sample data:

```r
mQ   <- matrix(c(1, 1, 2, 0), nrow = 2, byrow = TRUE)
pData <- read.csv(file.path("data", "AAPL_20080609_o.csv"))  # regular hours
result <- intradayJumpTest(pData, theta=0.8, p=4, L=10, mQ=mQ, q=3, varpi=0.49)
```

## Project layout

```
wspm/
├── R/
│   └── jump_test.R     # intradayJumpTest() wrapper
├── src/
│   └── subsampler.cpp  # Full C++ implementation (RcppArmadillo)
├── data/
│   ├── AAPL_20080609_e.csv   # AAPL tick data 2008-06-09, extended hours
│   └── AAPL_20080609_o.csv   # AAPL tick data 2008-06-09, regular hours
├── refs/
│   └── 1-s2.0-S0304405X25000182-main.pdf
├── load.R          # Convenience loader (source this to start)
└── wspm.Rproj
```

## Architecture

All computation lives in `src/subsampler.cpp`. `R/jump_test.R` exposes `intradayJumpTest()`, a thin wrapper that coerces the `pData` data frame to a numeric vector before calling the C++ function.

### C++ call hierarchy (`src/subsampler.cpp`)

```
intradayJumpTest_cpp()    ← main entry point (exported, called by R wrapper)
├── omega2_rcpp()         ← long-run noise variance (Jacod, Li, Zheng 2019)
├── f_subsampler_rcpp()   ← subsampled covariance estimator (Sigma*)
│   ├── aux_preavgk_rcpp()← pre-averaged return series (tent kernel)
│   └── f_mu_rcpp()       ← normalisation E[|Z|^r], Z ~ N(0,1)
└── f_psi_rcpp()          ← kernel constants psi1, psi2 for bias correction
```

### Key design choices

- The pre-averaging kernel is hard-coded as `g(x) = min(x, 1-x)` (tent function) everywhere in the C++ code.
- The bandwidth `kn = round(theta * sqrt(n))` is forced even (`kn + kn % 2`).
- Inside `f_subsampler_rcpp`, if the blocking parameter `c = floor(floor(n/(p*K))/L)` equals zero, `p` or `L` is decremented in a loop until `c > 0`. After this automatic adjustment the effective `p` and `L` may differ from the user's inputs.
- Realized variance (`mQ` row with `c(2,0)`) is recomputed **without** truncation at the end of `f_subsampler_rcpp`; bipower variation uses truncated series throughout.
- The long-run noise variance `omega2` uses two window lengths: `ell_n = round(n^(1/5))` for the autocovariance sum and `ell_n = round(n^(1/8))` for the Bartlett kernel weighting.
- The Bonferroni critical value is `qnorm(1 - 0.01/n)` (1% level).
- `intradayJumpTest` returns bias-corrected, annualised volatility measures (`RV`, `BV`) and the raw jump proportion (`JV`). `RV` and `BV` are floored at 0 before `sqrt`.

### Return values

| Element | Description |
|---------|-------------|
| `RV`    | Pre-averaged realized variance, bias-corrected and annualized as `100 * sqrt(252 * RV)` (%) |
| `BV`    | Pre-averaged bipower variation, bias-corrected and annualized as `100 * sqrt(252 * BV)` (%) |
| `JV`    | Jump proportion: `100 * (1 - BV_raw / RV_raw)` (%) |
| `JF`    | Jump indicator: `1` if the test rejects at the 1% Bonferroni level, `0` otherwise |

## Reference

The paper is in `refs/1-s2.0-S0304405X25000182-main.pdf`. Equations referenced in the code:

- Pre-averaged returns: eq. (9)
- `RV*`, `BV*`: eq. (11)
- Subsampled covariance `Sigma*`: eq. (18)
- Test statistic `J_n`: eq. (21)
- Bonferroni rejection rule: eq. (22)

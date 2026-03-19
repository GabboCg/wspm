# Warp Speed Price Moves — Intraday Jump Test

An R implementation of the noise-robust intraday jump test developed in:

> Christensen, K., Timmermann, A., & Veliyev, B. (2025). Warp speed price moves: Jumps after earnings announcements. *Journal of Financial Economics*, 167, 104010. https://doi.org/10.1016/j.jfineco.2025.104010

## Overview

Corporate earnings announcements should, in efficient markets, trigger near-instantaneous jumps in stock prices. Testing this at the tick-by-tick level requires a jump test that is **robust to the high levels of microstructure noise** characteristic of after-hours trading sessions. This package implements the pre-averaging, subsampling-based bipower variation jump test of Christensen et al. (2025), which generalizes the classical Barndorff-Nielsen and Shephard (2006) test.

The key innovations are:

- Pre-averaging observed noisy log-returns using a tent kernel `g(x) = min(x, 1 - x)` to attenuate microstructure noise

- A subsampling estimator for the asymptotic variance matrix, valid under heteroscedasticity and serial dependence in the noise process

- Bias correction for the long-run noise variance following Jacod, Li, and Zheng (2019)

- Truncation of the bipower covariance estimator to achieve robustness to large jumps

## Functions

### `intradayJumpTest(pData, theta, p, L, mQ, q, varpi)`

Main function. Computes the noise-robust jump test statistic for a single intraday price series.

**Arguments:**

| Argument | Description |
|----------|-------------|
| `pData`  | Data frame with a column `PRICE` containing tick-by-tick transaction prices |
| `theta`  | Pre-averaging window parameter. The bandwidth is `kn = round(theta * sqrt(n))` (rounded to even). Typical value: `0.5` |
| `p`      | Block size multiplier for the subsampler (number of pre-averaging windows per block). Typical value: `10` |
| `L`      | Number of subsamples used to estimate the asymptotic variance. Typical value: `10` |
| `mQ`     | A `2 × 2` matrix of power exponents. Each row `(q1, q2)` defines a variation measure. Row 1 should be `c(1, 1)` (bipower variation) and row 2 `c(2, 0)` (realized variance) |
| `q`      | Quantile multiplier for the truncation threshold. Typical value: `5` |
| `varpi`  | Truncation rate exponent. The threshold is `q * sqrt(BV) * n^(0.25 - varpi)`. Typical value: `0.24` |

**Returns:** A named list:

| Element | Description |
|---------|-------------|
| `RV`    | Pre-averaged realized variance, bias-corrected and annualized as `100 * sqrt(252 * RV)` (%) |
| `BV`    | Pre-averaged bipower variation, bias-corrected and annualized as `100 * sqrt(252 * BV)` (%) |
| `JV`    | Jump proportion: `100 * (1 - BV_raw / RV_raw)` (%) |
| `JF`    | Jump indicator: `1` if the test rejects at the 1% Bonferroni level, `0` otherwise |

### Internal functions (C++)

Implemented in `src/subsampler.cpp` via RcppArmadillo:

| Function | Description |
|----------|-------------|
| `aux_preavgk(logp, K)` | Pre-averaged return series (tent kernel hard-coded) |
| `f_mu(vpow)` | Computes `E[|Z|^r]` for `Z ~ N(0,1)` |
| `f_subsampler(logp, K, p, L, mQ, q, varpi)` | Subsampled power variation estimates and covariance matrix `Sigma*` |
| `f_psi(K)` | Kernel constants `psi1`, `psi2` for bias correction (tent kernel) |
| `omega2(logp)` | Long-run noise variance (Jacod, Li & Zheng 2019) |

## The Test Statistic

The jump test statistic (equation (21) in the paper) is:

```
J_n = n^(1/4) * (RV*_n - BV*_n(1,1)) / sqrt(v' * Sigma*_n * v)
```

where `v = [1, -1]'`, `RV*_n` is the pre-averaged realized variance, `BV*_n(1,1)` is the truncated pre-averaged bipower variation, and `Sigma*_n` is the subsampled covariance estimator. Under the null of no jumps, `J_n -> N(0,1)`.

The null is rejected at level `alpha` when `J_n > Phi^{-1}(1 - alpha/n)`, where the Bonferroni correction accounts for the `n` observation-level hypotheses.

## Usage

Load from the project root (compiles C++ and sources `R/jump_test.R`):

```r
source("load.R")
```

Run the test on the bundled AAPL data:

```r
mQ <- matrix(c(1, 2, 1, 0), nrow = 2, ncol = 2)

# Extended trading hours
pData  <- readr::read_csv("data/AAPL_20080609_e.csv")
result <- intradayJumpTest(pData, theta = 0.5, p = 10, L = 10, mQ = mQ, q = 5, varpi = 0.24)

# Regular trading hours
pData  <- readr::read_csv("data/AAPL_20080609_o.csv")
result <- intradayJumpTest(pData, theta = 0.5, p = 10, L = 10, mQ = mQ, q = 5, varpi = 0.24)

result$JF  # 1 = jump detected, 0 = no jump
result$RV  # annualized realized volatility (%)
result$BV  # annualized bipower variation (%)
result$JV  # jump proportion (%)
```

## Notes

- The function expects tick-time sampled data (all price changes retained). Observations with identical timestamps should be aggregated to an average price before calling the function, consistent with Griffin and Oomen (2008).

- The pre-averaging bandwidth `kn` is forced to be even, as required by the theory.

- If the data are too sparse relative to the chosen `p` and `L` (i.e., `c = floor(floor(n / (p * kn)) / L) == 0`), the function automatically reduces `p` or `L` to ensure at least one block is available.

- `RV` and `BV` in the output are bias-corrected for microstructure noise and floored at zero before annualization.

## References

- Barndorff-Nielsen, O.E., & Shephard, N. (2006). Econometrics of testing for jumps in financial economics using bipower variation. *Journal of Financial Econometrics*, 4(1), 1–30.

- Christensen, K., Oomen, R., & Podolskij, M. (2018). Fact or friction: Jumps at ultra high frequency. *Journal of Financial Economics*, 114(3), 576–599.

- Griffin, J.E., & Oomen, R.C.A. (2008). Sampling returns for realized variance calculations: Tick time or transaction time? *Econometric Reviews*, 27(1–3), 230–253.

- Jacod, J., Li, Y., & Zheng, X. (2019). Estimating the integrated volatility with tick observations. *Journal of Econometrics*, 208(1), 80–100.

- Jacod, J., Li, Y., Mykland, P.A., Podolskij, M., & Vetter, M. (2009). Microstructure noise in the continuous case: The pre-averaging approach. *Stochastic Processes and their Applications*, 119(7), 2249–2276.

- Podolskij, M., & Vetter, M. (2009). Estimation of volatility functionals in the simultaneous presence of microstructure noise and jumps. *Bernoulli*, 15(3), 634–658.

# ======================================================== #
#
#         Warp Speed Price Moves — Intraday Jump Test
#                C++ backend (RcppArmadillo) 
#
#                 Gabriel E. Cabrera-Guzmán
#                The University of Manchester
#
#                       Spring, 2026
#
#                https://gcabrerag.rbind.io
#
# ------------------------------ #
# email: gabriel.cabreraguzman@postgrad.manchester.ac.uk
# ======================================================== #

# C++ compilation is handled by load.R (Rcpp::sourceCpp must be called at the
# top level so the loaded DLL is not garbage-collected).  Do not call
# sourceCpp() here.


# -----------------------------------------------------------------------------
# intradayJumpTest()
#
# Noise-robust intraday jump test for a single price series.
#
# Args:
#   pData  data.frame with a column PRICE containing tick-by-tick prices.
#          Observations with identical timestamps should be averaged before
#          calling (Griffin & Oomen, 2008).
#   theta  Pre-averaging window parameter. Bandwidth kn = round(theta*sqrt(n)),
#          forced even. Typical value: 0.8.
#   p      Block size multiplier (number of pre-averaging windows per block).
#          Typical value: 4. Automatically reduced if data are too sparse.
#   L      Number of subsamples for the covariance estimator.
#          Typical value: 10. Automatically reduced if data are too sparse.
#   mQ     2x2 matrix of power exponents. Row 1 must be c(1,1) (bipower
#          variation); row 2 must be c(2,0) (realized variance).
#   q      Quantile multiplier for the truncation threshold. Typical value: 3.
#   varpi  Truncation rate exponent. Threshold = q*sqrt(BV)*n^(0.25-varpi).
#          Typical value: 0.49.
#
# Returns a named list:
#   RV   Bias-corrected, annualised realised volatility: 100*sqrt(252*RV) (%)
#   BV   Bias-corrected, annualised bipower variation:   100*sqrt(252*BV) (%)
#   JV   Jump proportion: 100*(1 - BV_raw/RV_raw) (%)
#   JF   Jump indicator: 1 if H0 rejected at 1% Bonferroni level, 0 otherwise
# -----------------------------------------------------------------------------
intradayJumpTest <- function(pData, theta, p, L, mQ, q, varpi) {
    
    intradayJumpTest_cpp(
        as.numeric(pData$PRICE),
        theta,
        as.integer(p),
        as.integer(L),
        mQ,
        q,
        varpi
    )
    
}

# =============================================================================
# wspm — project loader
#
# Sources all R modules in dependency order (mirrors devtools::load_all()).
# Run once per session from the project root:
#
#   source("load.R")
#
# =============================================================================

# Compile C++ and load exported symbols into the global environment.
# Must be called at top level (not nested inside source()) so the DLL
# is not garbage-collected before use.
Rcpp::sourceCpp(file.path("src", "subsampler.cpp"))

# Import intradayJumpTest() wrapper
source(file.path("R", "jump_test.R"))   

# Powers
mQ <- matrix(c(1, 2, 1, 0), nrow = 2, ncol = 2)

# Extended trading hour
pData <- readr::read_csv("data/AAPL_20080609_e.csv")

# Run intraday jump test (Christensen et al., 2025, JFE)
intradayJumpTest(pData, theta = 0.5, p = 10, L = 10, mQ = mQ, q = 5, varpi = 0.24)

# Regular trading hour
pData <- readr::read_csv("data/AAPL_20080609_o.csv")

# Run intraday jump test (Christensen et al., 2025, JFE)
intradayJumpTest(pData, theta = 0.5, p = 10, L = 10, mQ = mQ, q = 5, varpi = 0.24)

# tvBartSurv
R package implementing TV-BART-Cox from 'Time-Varying Treatments and Bayesian Additive Regression Trees with Cox Regression Model'. Flexible survival analysis with time-varying treatments using log-linear BART and a spline-based baseline via a hybrid Gibbs sampler.

**Dependency note**: This package depends on `loglinearBART`. If it is not on CRAN, the `Remotes:` field in `DESCRIPTION` points to `NingyaW/loglinearBART`. Change it if your repo is different, or remove the `Remotes:` line if `loglinearBART` is on CRAN.

## Install

```r
install.packages(c("remotes","devtools","roxygen2"))
# make sure loglinearBART is available:
# remotes::install_github("NingyaW/loglinearBART")
# or remotes::install_local("/path/to/loglinearBART.tar.gz", upgrade = "never")

devtools::document("tvBartSurv")
devtools::install("tvBartSurv")
```

## Example

```r
library(tvBartSurv)
set.seed(3)
res <- Simulation_Gibbs(N = 100, g = 1, a0 = "Weibull",
                        BARTsetting = 1, Mi = 5, iterations = 200)
str(res)
```

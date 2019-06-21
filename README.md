
# Multiple Discrete-Continuous Extreme Value (MDCEV) Model Estimation and Simulation in R: The rmdcev Package

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/plloydsmith/rmdcev.svg?branch=master)](https://travis-ci.org/plloydsmith/rmdcev)
<!-- badges: end -->

The *rmdcev* R package estimates and simulates multiple
discrete-continuous extreme value (MDCEV) demand models with observed
and unobserved individual heterogneity. Fixed parameter, latent class,
and random parameter models can be estimated. These models are estimated
using maximum likelihood or Bayesian estimation techniques and are
implemented in Stan, which is a C++ package for performing full Bayesian
inference (see <http://mc-stan.org/>). The package also includes
functions for simulating demand and welfare outcomes from policy
scenarios.

## Current Status

Development is in progress. Currently users can estimate the following
models:

1.  Fixed parameter models (maximum likelihood or Bayesian estimation)
2.  Latent class models (maximum likelihood estimation)
3.  Random parameters models

Under development

2.  Phaneuf/von Haefen Kuhn Tucker model

## Installation

I recommend you first install **rstan** and C++ toolchain by following
these steps:

<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>

Once **rstan** is installed, you can install the released version of
rmdcev from GitHub using devtools

``` r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("plloydsmith/rmdcev", build_vignettes = FALSE)
```

Depending on your computer set-up, you may need to adjust your Makevars
file of the .R folder (usually in your computer user folder) to ensure
the first two lines are

``` r
CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y
CXX14FLAGS=-O3 -mtune=native -march=native
```

You can switch build\_vignettes to TRUE but it will take a lot longer to
install (Note: The vignette will be complete soon). If installation
fails, please let me know by [filing an
issue](https://github.com/plloydsmith/rmdcev/issues).

## References

For more details on the model specification and estimation:

Bhat, C.R. (2008) [“The Multiple Discrete-Continuous Extreme Value
(MDCEV) Model: Role of Utility Function Parameters, Identification
Considerations, and Model
Extensions.”](https://www.sciencedirect.com/science/article/pii/S0191261507000677)
Transportation Research Part B, 42(3): 274-303.

For more details on the welfare simulation:

Lloyd-Smith, P (2018). [“A New Approach to Calculating Welfare Measures
in Kuhn-Tucker Demand
Models.”](https://www.sciencedirect.com/science/article/pii/S1755534517300994)
Journal of Choice Modeling, 26: 19-27

## Estimation

As an example, we can simulate some data using the Hybrid specification.
In this example, we are simulating data for 1000 invdividuals and 10
non-numeraire alternatives.

``` r
library(pacman)
p_load(tidyverse, rmdcev)
model <- "hybrid"
nobs <- 1000
ngoods <- 10
sim.data <- GenerateMDCEVData(model = model, nobs = nobs, ngoods = ngoods)
```

Estimate model using MLE

``` r
mdcev_est <- FitMDCEV(psi_formula = ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8-1,
                   data = sim.data$data,
                   model = model,
                   algorithm = "MLE")
#> Checking data...
#> Data is good
#> Using MLE to estimate MDCEV
#> Chain 1: Initial log joint probability = -17885.5
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -14869.2      0.801289       238.643           1           1       24   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -14763.5      0.346585       155.473           1           1       47   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -14752.5     0.0310867       9.44744           1           1       68   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79        -14751     0.0165834       6.64363           1           1       90   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -14750.9     0.0228325       1.88743           1           1      114   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      119      -14750.9    0.00349279       1.57032           1           1      136   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      127      -14750.9   0.000798697      0.365787           1           1      145   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
SummaryMDCEV(mdcev_est)
#> Model run using rmdcev for R, version 0.9.0 
#> Estimation method                : MLE
#> Model type                       : hybrid specification
#> Number of classes                : 1
#> Number of individuals            : 1000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -14750.88
#> AIC                              : 29541.77
#> BIC                              : 29639.92
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.69
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 21.89 15.91  6.16 23.71  7.09 25.39 15.10 32.50 40.47 13.27 
#> 
#> 
#> Psi specification:
#> ~
#> b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 - 1
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -6.022   0.690  -8.72
#> psi_b2     0.552   0.103   5.33
#> psi_b3     2.220   0.148  15.01
#> psi_b4    -1.612   0.122 -13.24
#> psi_b5     3.209   0.206  15.57
#> psi_b6    -2.190   0.150 -14.64
#> psi_b7     1.090   0.116   9.41
#> psi_b8     2.280   0.143  15.90
#> gamma1     1.084   0.143   7.60
#> gamma2     1.000   0.139   7.17
#> gamma3     1.601   0.356   4.50
#> gamma4     2.128   0.353   6.02
#> gamma5     0.784   0.115   6.83
#> gamma6     1.867   0.306   6.09
#> gamma7     2.304   0.349   6.60
#> gamma8     1.839   0.250   7.34
#> gamma9     1.131   0.146   7.73
#> gamma10    1.130   0.207   5.45
#> alpha1     0.454   0.034  13.30
#> scale      1.093   0.065  16.74
#> Note: Alpha parameter is equal for all goods.
```

Compare estimates to true values

``` r
parms_true <- tbl_df(sim.data$parms_true) %>%
    mutate(true = as.numeric(true))

coefs <- mdcev_est$est_pars %>%
    mutate(parms = gsub("\\[|\\]", "", parms)) %>%
    group_by(parms) %>%
    summarise(mean = mean(value),
              sd = sd(value),
              zstat = mean / sd,
              cl_lo = quantile(value, 0.025),
              cl_hi = quantile(value, 0.975)) %>%
    left_join(parms_true, by = "parms") %>%
    print(n=200)
#> # A tibble: 20 x 7
#>    parms     mean     sd  zstat  cl_lo  cl_hi  true
#>    <chr>    <dbl>  <dbl>  <dbl>  <dbl>  <dbl> <dbl>
#>  1 alpha1   0.454 0.0342  13.3   0.386  0.508  0.5 
#>  2 gamma1   1.08  0.143    7.60  0.882  1.38   1.16
#>  3 gamma10  1.13  0.207    5.45  0.785  1.57   1.11
#>  4 gamma2   1.000 0.139    7.17  0.728  1.22   1.04
#>  5 gamma3   1.60  0.356    4.50  1.09   2.34   1.82
#>  6 gamma4   2.13  0.353    6.02  1.50   2.74   1.76
#>  7 gamma5   0.784 0.115    6.83  0.585  0.991  1.12
#>  8 gamma6   1.87  0.306    6.09  1.38   2.41   1.92
#>  9 gamma7   2.30  0.349    6.60  1.82   2.98   1.84
#> 10 gamma8   1.84  0.250    7.34  1.42   2.37   1.64
#> 11 gamma9   1.13  0.146    7.73  0.850  1.38   1.11
#> 12 psi_b1  -6.02  0.690   -8.72 -7.43  -4.86  -5   
#> 13 psi_b2   0.552 0.103    5.33  0.376  0.718  0.5 
#> 14 psi_b3   2.22  0.148   15.0   1.94   2.51   2   
#> 15 psi_b4  -1.61  0.122  -13.2  -1.82  -1.38  -1.5 
#> 16 psi_b5   3.21  0.206   15.6   2.90   3.68   3   
#> 17 psi_b6  -2.19  0.150  -14.6  -2.46  -1.91  -2   
#> 18 psi_b7   1.09  0.116    9.41  0.889  1.27   1   
#> 19 psi_b8   2.28  0.143   15.9   2.03   2.52   2   
#> 20 scale    1.09  0.0653  16.7   0.992  1.21   1
```

Compare outputs using a figure

``` r
coefs %>%
    ggplot(aes(y = mean, x = true))  +
    geom_point(size=2) +
    geom_text(label=coefs$parms) +
    geom_abline(slope = 1) +
    geom_errorbar(aes(ymin=cl_lo,ymax=cl_hi,width=0.2))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

## Welfare simulation

Create policy simulations (these are ‘no change’ policies with no
effects)

``` r
npols <- 2 # Choose number of policies

policies<-  CreateBlankPolicies(npols, ngoods, mdcev_est$stan_data[["dat_psi"]], price_change_only = TRUE)

df_sim <- PrepareSimulationData(mdcev_est, policies)
```

Simulate welfare changes

``` r
wtp <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
                     cond_err = 1, nerrs = 5, sim_type = "welfare")
#> Using hybrid approach to simulation
#> Compiling simulation code
#> Simulating welfare...
#> 
#> 3.00e+05simulations finished in0.06minutes.(77519per second)
SummaryWelfare(wtp)
#> # A tibble: 2 x 5
#>   policy   Mean Std.Dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>   <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1     0       0           0            0
#> 2 policy2     0       0           0            0
```

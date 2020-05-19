
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

1.  Bhat (2008) MDCEV model specifications
2.  Kuhn-Tucker model specification (von Haefen and Phaneuf, 2005)

Models can be estimated using

1.  Fixed parameter models (maximum likelihood or Bayesian estimation)
2.  Latent class models (maximum likelihood estimation)
3.  Random parameters models (Bayesian estimation, except KT model
    specifications)

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
install_github("plloydsmith/rmdcev", build_vignettes = FALSE, INSTALL_opts="--no-multiarch")
```

You can switch build\_vignettes to TRUE but it will take a lot longer to
install. If installation fails, please let me know by [filing an
issue](https://github.com/plloydsmith/rmdcev/issues).

## References

For more details on the model specification and estimation:

Bhat, C.R. (2008) [“The Multiple Discrete-Continuous Extreme Value
(MDCEV) Model: Role of Utility Function Parameters, Identification
Considerations, and Model
Extensions”](https://www.sciencedirect.com/science/article/pii/S0191261507000677)
Transportation Research Part B, 42(3): 274-303.

von Haefen, R. and Phaneuf D. (2005) [“Kuhn-Tucker Demand System
Approaches to Non-Market
Valuation”](https://link.springer.com/chapter/10.1007/1-4020-3684-1_8)
In: Scarpa R., Alberini A. (eds) Applications of Simulation Methods in
Environmental and Resource Economics. The Economics of Non-Market Goods
and Resources, vol 6. Springer, Dordrecht.

For more details on the welfare simulation:

Lloyd-Smith, P (2018). [“A New Approach to Calculating Welfare Measures
in Kuhn-Tucker Demand
Models.”](https://www.sciencedirect.com/science/article/pii/S1755534517300994)
Journal of Choice Modeling, 26: 19-27

## Estimation

As an example, we can simulate some data using Bhat (2008)‘s ’Gamma’
specification. In this example, we are simulating data for 2,000
individuals and 10 non-numeraire alternatives.

``` r
library(pacman)
p_load(tidyverse, rmdcev)
model <- "gamma"
nobs <- 2000
nalts <- 10
sim.data <- GenerateMDCEVData(model = model, nobs = nobs, nalts = nalts)
#> Checking data...
#> Data is good
```

Estimate model using MLE

``` r
mdcev_est <- mdcev(~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8-1,
                   data = sim.data$data,
                   model = model,
                   std_errors = "deltamethod",
                   algorithm = "MLE")
#> Using MLE to estimate MDCEV
#> Chain 1: Initial log joint probability = -51663.6
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -38808.9       1.74731       438.401           1           1       24   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -38330.4       0.27699       72.9864           1           1       49   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -38325.9    0.00620657       9.51927           1           1       74   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       76      -38325.8    0.00207364      0.389301           1           1       94   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.1.2 
#> Estimation method                : MLE
#> Model type                       : gamma specification
#> Number of classes                : 1
#> Number of individuals            : 2000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -38325.79
#> AIC                              : 76691.59
#> BIC                              : 76803.61
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:2.59
#> 
#> Average consumption of non-numeraire alternatives:
#>      1      2      3      4      5      6      7      8      9     10 
#> 169.86  11.77   0.76  40.80   4.12   4.92  52.11   6.63   7.22   1.07 
#> 
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -4.804   0.098 -49.02
#> psi_b2     0.367   0.116   3.16
#> psi_b3     1.899   0.086  22.09
#> psi_b4    -1.494   0.055 -27.16
#> psi_b5     3.010   0.051  59.02
#> psi_b6    -1.992   0.059 -33.76
#> psi_b7     0.994   0.043  23.12
#> psi_b8     1.976   0.045  43.91
#> gamma1     8.858   0.643  13.78
#> gamma2     9.318   0.687  13.56
#> gamma3     2.987   0.467   6.40
#> gamma4     5.570   0.298  18.69
#> gamma5     9.772   1.067   9.16
#> gamma6     1.660   0.108  15.37
#> gamma7     5.289   0.336  15.74
#> gamma8     6.580   0.569  11.56
#> gamma9     4.784   0.357  13.40
#> gamma10    3.378   0.448   7.54
#> alpha1     0.497   0.008  62.16
#> scale      0.969   0.013  74.57
#> Note: All non-numeraire alpha's fixed to 0.
```

Compare estimates to true values

``` r
coefs <- as_tibble(sim.data$parms_true) %>%
    mutate(true = as.numeric(true)) %>%
 cbind(summary(mdcev_est)[["CoefTable"]]) %>%
    mutate(cl_lo = Estimate - 1.96 * Std.err,
           cl_hi = Estimate + 1.96 * Std.err)

head(coefs, 200)
#>           parms      true Estimate Std.err z.stat    cl_lo    cl_hi
#> psi_b1   psi_b1 -5.000000   -4.804   0.098 -49.02 -4.99608 -4.61192
#> psi_b2   psi_b2  0.500000    0.367   0.116   3.16  0.13964  0.59436
#> psi_b3   psi_b3  2.000000    1.899   0.086  22.09  1.73044  2.06756
#> psi_b4   psi_b4 -1.500000   -1.494   0.055 -27.16 -1.60180 -1.38620
#> psi_b5   psi_b5  3.000000    3.010   0.051  59.02  2.91004  3.10996
#> psi_b6   psi_b6 -2.000000   -1.992   0.059 -33.76 -2.10764 -1.87636
#> psi_b7   psi_b7  1.000000    0.994   0.043  23.12  0.90972  1.07828
#> psi_b8   psi_b8  2.000000    1.976   0.045  43.91  1.88780  2.06420
#> gamma1   gamma1  8.390458    8.858   0.643  13.78  7.59772 10.11828
#> gamma2   gamma2  8.583333    9.318   0.687  13.56  7.97148 10.66452
#> gamma3   gamma3  3.290057    2.987   0.467   6.40  2.07168  3.90232
#> gamma4   gamma4  4.803377    5.570   0.298  18.69  4.98592  6.15408
#> gamma5   gamma5  9.502540    9.772   1.067   9.16  7.68068 11.86332
#> gamma6   gamma6  1.576536    1.660   0.108  15.37  1.44832  1.87168
#> gamma7   gamma7  5.370277    5.289   0.336  15.74  4.63044  5.94756
#> gamma8   gamma8  6.842473    6.580   0.569  11.56  5.46476  7.69524
#> gamma9   gamma9  4.657651    4.784   0.357  13.40  4.08428  5.48372
#> gamma10 gamma10  2.686770    3.378   0.448   7.54  2.49992  4.25608
#> alpha1   alpha1  0.500000    0.497   0.008  62.16  0.48132  0.51268
#> scale     scale  1.000000    0.969   0.013  74.57  0.94352  0.99448
```

Compare outputs using a figure

``` r
coefs %>%
    ggplot(aes(y = Estimate, x = true))  +
    geom_point(size=2) +
    geom_text(label=coefs$parms,position=position_jitter(width=.5,height=1)) +
    geom_abline(slope = 1) +
    geom_errorbar(aes(ymin=cl_lo,ymax=cl_hi,width=0.2))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />

## Welfare simulation

Create policy simulations (these are ‘no change’ policies with no
effects)

``` r
npols <- 2 # Choose number of policies

policies<-  CreateBlankPolicies(npols, nalts, mdcev_est$stan_data[["dat_psi"]])

df_sim <- PrepareSimulationData(mdcev_est, policies, nsims = 1)
```

Simulate welfare changes

``` r
wtp <- mdcev.sim(df_sim$df_indiv, 
                 df_common = df_sim$df_common, 
                 sim_options = df_sim$sim_options,
                 cond_err = 1, 
                 nerrs = 15, 
                 sim_type = "welfare")
#> Using general approach to simulation
#> Simulating welfare...
#> 
#> 6.00e+04simulations finished in0.35minutes.(2872per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy      mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#> * <chr>      <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 7.13e-11      NA    7.13e-11     7.13e-11
#> 2 policy2 7.13e-11      NA    7.13e-11     7.13e-11
```

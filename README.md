
<!-- README.md is generated from README.Rmd. Please edit that file -->

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
#> Chain 1: Initial log joint probability = -17583.3
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -14521.4      0.136714        65.618           1           1       25   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39        -14448      0.142581          37.5           1           1       47   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59        -14440     0.0781948       38.7772      0.5104           1       69   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -14438.4     0.0653246       14.4014           1           1       91   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -14438.1     0.0040228       1.51202      0.8957      0.8957      114   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      119      -14438.1    0.00338588       1.02323           1           1      135   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      131      -14438.1   0.000646157      0.096187           1           1      148   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
SummaryMDCEV(mdcev_est)
#> Model run using rmdcev for R, version 0.7.0 
#> Estimation method                : MLE
#> Model type                       : hybrid specification
#> Number of classes                : 1
#> Number of individuals            : 1000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -14438.1
#> AIC                              : 28916.2
#> BIC                              : 29014.35
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.7
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 30.98 15.01  9.68 38.77  6.81 11.71  8.65  7.29 32.65 23.89 
#> 
#> 
#> Psi specification:
#> ~
#> b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 - 1
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -5.042   0.464 -10.87
#> psi_b2     0.448   0.103   4.33
#> psi_b3     1.903   0.142  13.42
#> psi_b4    -1.624   0.108 -14.98
#> psi_b5     3.034   0.150  20.17
#> psi_b6    -1.907   0.119 -16.07
#> psi_b7     1.066   0.104  10.27
#> psi_b8     1.993   0.109  18.23
#> gamma1     1.423   0.180   7.90
#> gamma2     1.454   0.207   7.01
#> gamma3     0.864   0.180   4.80
#> gamma4     2.355   0.305   7.72
#> gamma5     1.283   0.289   4.44
#> gamma6     1.464   0.225   6.50
#> gamma7     1.286   0.228   5.64
#> gamma8     1.839   0.393   4.68
#> gamma9     1.456   0.215   6.76
#> gamma10    2.216   0.293   7.57
#> alpha1     0.490   0.021  23.30
#> scale      1.003   0.042  23.71
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
#>  1 alpha1   0.490 0.0210  23.3   0.451  0.529  0.5 
#>  2 gamma1   1.42  0.180    7.90  1.17   1.82   1.14
#>  3 gamma10  2.22  0.293    7.57  1.70   2.73   1.87
#>  4 gamma2   1.45  0.207    7.01  1.11   1.80   1.18
#>  5 gamma3   0.864 0.180    4.80  0.565  1.19   1.11
#>  6 gamma4   2.36  0.305    7.72  1.92   3.15   1.87
#>  7 gamma5   1.28  0.289    4.44  0.890  1.79   1.32
#>  8 gamma6   1.46  0.225    6.50  1.10   1.90   1.68
#>  9 gamma7   1.29  0.228    5.64  0.913  1.68   1.37
#> 10 gamma8   1.84  0.393    4.68  1.18   2.61   1.17
#> 11 gamma9   1.46  0.215    6.76  1.08   1.87   1.33
#> 12 psi_b1  -5.04  0.464  -10.9  -6.00  -4.13  -5   
#> 13 psi_b2   0.448 0.103    4.33  0.244  0.606  0.5 
#> 14 psi_b3   1.90  0.142   13.4   1.63   2.13   2   
#> 15 psi_b4  -1.62  0.108  -15.0  -1.82  -1.42  -1.5 
#> 16 psi_b5   3.03  0.150   20.2   2.73   3.32   3   
#> 17 psi_b6  -1.91  0.119  -16.1  -2.12  -1.72  -2   
#> 18 psi_b7   1.07  0.104   10.3   0.877  1.24   1   
#> 19 psi_b8   1.99  0.109   18.2   1.82   2.23   2   
#> 20 scale    1.00  0.0423  23.7   0.924  1.07  NA
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
#>  3.00e+05 simulations finished in 0.06 minutes. ( 85470 per second)
SummaryWelfare(wtp)
#> # A tibble: 2 x 5
#>   policy   Mean Std.Dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>   <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1     0       0           0            0
#> 2 policy2     0       0           0            0
```

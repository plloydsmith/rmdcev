
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
3.  Random parameters models (Bayesian estimation)

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
install_github("plloydsmith/rmdcev", build_vignettes = FALSE, INSTALL_opts="--no-multiarch")
```

Depending on your computer set-up, you may need to adjust your Makevars
file of the .R folder (usually in your computer user folder) to ensure
the first two lines are

``` r
CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y
CXX14FLAGS=-O3 -mtune=native -march=native
```

You can switch build\_vignettes to TRUE but it will take a lot longer to
install. If installation fails, please let me know by [filing an
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
                   algorithm = "MLE")
#> Using MLE to estimate MDCEV
#> Chain 1: Initial log joint probability = -18631.4
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -15626.9      0.185145       107.315           1           1       22   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39        -15555       0.14702       84.1025           1           1       43   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -15542.6     0.0775889        40.987           1           1       65   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -15541.4     0.0174648       7.09082           1           1       89   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -15540.9     0.0284867       3.21212           1           1      111   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      119      -15540.9     0.0052752       1.16396           1           1      133   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      128      -15540.9     0.0016659      0.312295           1           1      142   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.1.2 
#> Estimation method                : MLE
#> Model type                       : hybrid specification
#> Number of classes                : 1
#> Number of individuals            : 1000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -15540.88
#> AIC                              : 31121.77
#> BIC                              : 31219.92
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.54
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 25.51 15.15  9.70 23.00 37.43 25.60 12.60 12.74 27.76  9.32 
#> 
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -5.514   0.547 -10.09
#> psi_b2     0.491   0.076   6.45
#> psi_b3     1.986   0.107  18.49
#> psi_b4    -1.433   0.102 -14.10
#> psi_b5     3.190   0.175  18.28
#> psi_b6    -2.126   0.122 -17.50
#> psi_b7     0.988   0.085  11.63
#> psi_b8     2.053   0.124  16.50
#> gamma1     1.323   0.209   6.34
#> gamma2     1.459   0.195   7.50
#> gamma3     1.136   0.219   5.18
#> gamma4     1.522   0.312   4.87
#> gamma5     2.177   0.324   6.71
#> gamma6     1.692   0.212   8.00
#> gamma7     1.260   0.205   6.14
#> gamma8     2.001   0.342   5.85
#> gamma9     1.331   0.193   6.90
#> gamma10    1.232   0.249   4.95
#> alpha1     0.473   0.027  17.73
#> scale      1.038   0.046  22.73
#> Note: Alpha parameter is equal for all alternatives.
```

Compare estimates to true values

``` r
parms_true <- tbl_df(sim.data$parms_true) %>%
    mutate(true = as.numeric(true))

output <- tbl_df(mdcev_est[["stan_fit"]][["theta_tilde"]]) %>%
                dplyr::select(-tidyselect::starts_with("log_like"),         
                              -tidyselect::starts_with("sum_log_lik"))

names(output)[1:mdcev_est$parms_info$n_vars$n_parms_total] <- mdcev_est$parms_info$parm_names$all_names

output<- output %>%
        tibble::rowid_to_column("sim_id") %>%
        tidyr::gather(parms, value, -sim_id)

coefs <- output %>%
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
#>  1 alpha1   0.473 0.0267  17.7   0.427  0.526  0.5 
#>  2 gamma1   1.32  0.209    6.34  0.972  1.71   1.33
#>  3 gamma10  1.23  0.249    4.95  0.880  1.90   1.13
#>  4 gamma2   1.46  0.195    7.50  1.13   1.93   1.72
#>  5 gamma3   1.14  0.219    5.18  0.813  1.57   1.69
#>  6 gamma4   1.52  0.312    4.87  1.09   2.22   1.61
#>  7 gamma5   2.18  0.324    6.71  1.67   2.72   1.75
#>  8 gamma6   1.69  0.212    8.00  1.34   2.14   1.51
#>  9 gamma7   1.26  0.205    6.14  0.911  1.70   1.21
#> 10 gamma8   2.00  0.342    5.85  1.53   2.76   1.91
#> 11 gamma9   1.33  0.193    6.90  0.987  1.64   1.48
#> 12 psi_b1  -5.51  0.547  -10.1  -6.36  -4.39  -5   
#> 13 psi_b2   0.491 0.0761   6.45  0.364  0.628  0.5 
#> 14 psi_b3   1.99  0.107   18.5   1.82   2.19   2   
#> 15 psi_b4  -1.43  0.102  -14.1  -1.62  -1.26  -1.5 
#> 16 psi_b5   3.19  0.175   18.3   2.87   3.54   3   
#> 17 psi_b6  -2.13  0.122  -17.5  -2.32  -1.90  -2   
#> 18 psi_b7   0.988 0.0849  11.6   0.819  1.13   1   
#> 19 psi_b8   2.05  0.124   16.5   1.80   2.25   2   
#> 20 scale    1.04  0.0456  22.7   0.949  1.12   1
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

policies<-  CreateBlankPolicies(npols, nalts, mdcev_est$stan_data[["dat_psi"]], price_change_only = TRUE)

df_sim <- PrepareSimulationData(mdcev_est, policies)
```

Simulate welfare changes

``` r
wtp <- mdcev.sim(df_sim$df_indiv, 
                 df_common = df_sim$df_common, 
                 sim_options = df_sim$sim_options,
                 cond_err = 1, 
                 nerrs = 15, 
                 sim_type = "welfare")
#> Using hybrid approach to simulation
#> Simulating welfare...
#> 
#> 9.00e+05simulations finished in0.18minutes.(82569per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy      mean  std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>      <dbl>    <dbl>       <dbl>        <dbl>
#> 1 policy1 3.84e-12 7.79e-11   -1.16e-10     1.44e-10
#> 2 policy2 3.84e-12 7.79e-11   -1.16e-10     1.44e-10
```

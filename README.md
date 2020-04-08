
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
#> Chain 1: Initial log joint probability = -17556.8
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -14641.9      0.204815       108.929           1           1       24   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -14570.9      0.113784       101.523      0.2645      0.2645       45   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -14559.4     0.0131136       26.8481        0.55        0.55       67   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79        -14558     0.0120502       8.84637           1           1       89   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -14557.1    0.00784037       3.07293           1           1      110   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      115      -14557.1   0.000170927      0.247455      0.3223           1      129   
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
#> LL                               : -14557.05
#> AIC                              : 29154.1
#> BIC                              : 29252.26
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.37
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 22.12  9.23 11.89 13.22 17.99 20.27 12.60 17.36 33.84 27.90 
#> 
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -6.173   0.592 -10.42
#> psi_b2     0.521   0.073   7.13
#> psi_b3     2.042   0.189  10.79
#> psi_b4    -1.485   0.126 -11.79
#> psi_b5     3.325   0.180  18.47
#> psi_b6    -2.203   0.135 -16.30
#> psi_b7     1.218   0.098  12.38
#> psi_b8     2.176   0.142  15.29
#> gamma1     1.596   0.279   5.72
#> gamma2     1.587   0.266   5.96
#> gamma3     1.175   0.199   5.90
#> gamma4     1.989   0.269   7.38
#> gamma5     1.207   0.148   8.14
#> gamma6     1.345   0.186   7.24
#> gamma7     1.252   0.184   6.80
#> gamma8     2.135   0.318   6.71
#> gamma9     1.825   0.247   7.39
#> gamma10    2.076   0.256   8.12
#> alpha1     0.454   0.028  16.24
#> scale      1.077   0.055  19.61
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
#>  1 alpha1   0.454 0.0280  16.2   0.407  0.511  0.5 
#>  2 gamma1   1.60  0.279    5.72  1.14   2.14   1.48
#>  3 gamma10  2.08  0.256    8.12  1.68   2.64   1.71
#>  4 gamma2   1.59  0.266    5.96  1.20   2.15   1.69
#>  5 gamma3   1.18  0.199    5.90  0.842  1.56   1.77
#>  6 gamma4   1.99  0.269    7.38  1.50   2.51   1.87
#>  7 gamma5   1.21  0.148    8.14  0.979  1.49   1.16
#>  8 gamma6   1.34  0.186    7.24  1.10   1.68   1.39
#>  9 gamma7   1.25  0.184    6.80  0.936  1.61   1.40
#> 10 gamma8   2.14  0.318    6.71  1.56   2.76   1.69
#> 11 gamma9   1.82  0.247    7.39  1.41   2.31   1.55
#> 12 psi_b1  -6.17  0.592  -10.4  -7.37  -4.99  -5   
#> 13 psi_b2   0.521 0.0731   7.13  0.371  0.629  0.5 
#> 14 psi_b3   2.04  0.189   10.8   1.70   2.38   2   
#> 15 psi_b4  -1.49  0.126  -11.8  -1.70  -1.24  -1.5 
#> 16 psi_b5   3.32  0.180   18.5   2.95   3.67   3   
#> 17 psi_b6  -2.20  0.135  -16.3  -2.48  -1.98  -2   
#> 18 psi_b7   1.22  0.0983  12.4   1.05   1.38   1   
#> 19 psi_b8   2.18  0.142   15.3   1.89   2.46   2   
#> 20 scale    1.08  0.0549  19.6   0.980  1.17   1
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
#> 9.00e+05simulations finished in0.2minutes.(75949per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy       mean  std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>       <dbl>    <dbl>       <dbl>        <dbl>
#> 1 policy1 -4.41e-12 6.64e-11   -1.14e-10     1.28e-10
#> 2 policy2 -4.41e-12 6.64e-11   -1.14e-10     1.28e-10
```


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
#> Chain 1: Initial log joint probability = -9053.19
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -7551.78      0.145059       58.4705           1           1       25   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -7488.05     0.0472515       33.6947      0.1018           1       50   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59       -7482.5     0.0177925       8.83702           1           1       72   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -7481.67     0.0428037       5.14411           1           1       92   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -7481.63    0.00789223       1.20729           1           1      114   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      119      -7481.62   0.000841962      0.239115           1           1      136   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      128      -7481.62     0.0004145      0.113841           1           1      147   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.1.0 
#> Estimation method                : MLE
#> Model type                       : hybrid specification
#> Number of classes                : 1
#> Number of individuals            : 499
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -7481.62
#> AIC                              : 15003.24
#> BIC                              : 15087.49
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:0.75
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 13.59 12.00  6.94 10.57 12.21 17.57 27.49 32.82 22.63 39.14 
#> 
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -5.209   0.731  -7.12
#> psi_b2     0.403   0.081   5.00
#> psi_b3     2.065   0.206  10.04
#> psi_b4    -1.669   0.160 -10.42
#> psi_b5     3.049   0.210  14.53
#> psi_b6    -1.955   0.155 -12.62
#> psi_b7     1.214   0.128   9.49
#> psi_b8     1.953   0.167  11.70
#> gamma1     2.331   0.541   4.31
#> gamma2     1.779   0.471   3.78
#> gamma3     1.692   0.471   3.59
#> gamma4     1.092   0.269   4.06
#> gamma5     1.371   0.288   4.76
#> gamma6     1.173   0.218   5.38
#> gamma7     1.407   0.272   5.16
#> gamma8     2.470   0.482   5.12
#> gamma9     1.356   0.296   4.58
#> gamma10    1.500   0.278   5.39
#> alpha1     0.483   0.036  13.31
#> scale      0.999   0.065  15.35
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
#>  1 alpha1   0.483 0.0363  13.3   0.427  0.552  0.5 
#>  2 gamma1   2.33  0.541    4.31  1.61   3.54   1.78
#>  3 gamma10  1.50  0.278    5.39  1.01   2.08   1.14
#>  4 gamma2   1.78  0.471    3.78  1.08   2.62   1.49
#>  5 gamma3   1.69  0.471    3.59  0.933  2.72   1.55
#>  6 gamma4   1.09  0.269    4.06  0.693  1.64   1.18
#>  7 gamma5   1.37  0.288    4.76  0.933  2.00   1.29
#>  8 gamma6   1.17  0.218    5.38  0.791  1.55   1.48
#>  9 gamma7   1.41  0.272    5.16  0.982  1.82   1.68
#> 10 gamma8   2.47  0.482    5.12  1.53   3.34   1.67
#> 11 gamma9   1.36  0.296    4.58  0.886  2.01   1.27
#> 12 psi_b1  -5.21  0.731   -7.12 -6.44  -3.70  -5   
#> 13 psi_b2   0.403 0.0806   5.00  0.262  0.530  0.5 
#> 14 psi_b3   2.07  0.206   10.0   1.73   2.46   2   
#> 15 psi_b4  -1.67  0.160  -10.4  -1.93  -1.39  -1.5 
#> 16 psi_b5   3.05  0.210   14.5   2.63   3.41   3   
#> 17 psi_b6  -1.95  0.155  -12.6  -2.18  -1.66  -2   
#> 18 psi_b7   1.21  0.128    9.49  1.01   1.44   1   
#> 19 psi_b8   1.95  0.167   11.7   1.64   2.23   2   
#> 20 scale    0.999 0.0651  15.3   0.883  1.12   1
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
                 nerrs = 5, 
                 sim_type = "welfare")
#> Using hybrid approach to simulation
#> Compiling simulation code
#> Simulating welfare...
#> 
#> 1.50e+05simulations finished in0.04minutes.(70948per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy       mean  std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>       <dbl>    <dbl>       <dbl>        <dbl>
#> 1 policy1 -8.90e-12 7.81e-11   -1.33e-10     1.30e-10
#> 2 policy2 -8.90e-12 7.81e-11   -1.33e-10     1.30e-10
```

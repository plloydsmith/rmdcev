
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
#> Chain 1: Initial log joint probability = -17144.7
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -14148.5      0.664993       279.006           1           1       25   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -14080.6     0.0979007       35.6496           1           1       47   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -14072.3     0.0364687       21.5479      0.2996           1       75   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -14071.7    0.00679215       2.50178      0.2447           1       98   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -14071.6    0.00268428      0.585792           1           1      121   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      111      -14071.6   0.000172627      0.423097      0.5706      0.5706      134   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
SummaryMDCEV(mdcev_est)
#> Model run using rmdcev for R, version 0.8.0 
#> Estimation method                : MLE
#> Model type                       : hybrid specification
#> Number of classes                : 1
#> Number of individuals            : 1000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -14071.61
#> AIC                              : 28183.21
#> BIC                              : 28281.37
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.58
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 25.78 29.63  6.75 31.63 23.53 14.60 13.83  6.79 15.70 22.92 
#> 
#> 
#> Psi specification:
#> ~
#> b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 - 1
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -4.431   0.448  -9.88
#> psi_b2     0.602   0.091   6.62
#> psi_b3     1.769   0.138  12.78
#> psi_b4    -1.321   0.107 -12.35
#> psi_b5     2.838   0.138  20.58
#> psi_b6    -1.845   0.100 -18.45
#> psi_b7     0.914   0.068  13.35
#> psi_b8     1.872   0.092  20.28
#> gamma1     1.753   0.296   5.93
#> gamma2     1.610   0.218   7.38
#> gamma3     0.701   0.161   4.36
#> gamma4     1.340   0.180   7.46
#> gamma5     1.187   0.174   6.80
#> gamma6     1.304   0.235   5.55
#> gamma7     1.743   0.269   6.48
#> gamma8     1.185   0.252   4.71
#> gamma9     1.560   0.257   6.06
#> gamma10    1.654   0.246   6.73
#> alpha1     0.532   0.023  22.71
#> scale      0.945   0.039  24.48
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
#>  1 alpha1   0.532 0.0234  22.7   0.493  0.579  0.5 
#>  2 gamma1   1.75  0.296    5.93  1.36   2.45   1.74
#>  3 gamma10  1.65  0.246    6.73  1.21   2.16   1.62
#>  4 gamma2   1.61  0.218    7.38  1.27   2.00   1.93
#>  5 gamma3   0.701 0.161    4.36  0.420  1.02   1.06
#>  6 gamma4   1.34  0.180    7.46  1.07   1.69   1.10
#>  7 gamma5   1.19  0.174    6.80  0.872  1.51   1.24
#>  8 gamma6   1.30  0.235    5.55  0.956  1.70   1.28
#>  9 gamma7   1.74  0.269    6.48  1.27   2.24   1.89
#> 10 gamma8   1.18  0.252    4.71  0.846  1.85   1.26
#> 11 gamma9   1.56  0.257    6.06  1.14   2.10   1.97
#> 12 psi_b1  -4.43  0.448   -9.88 -5.17  -3.73  -5   
#> 13 psi_b2   0.602 0.0909   6.62  0.460  0.799  0.5 
#> 14 psi_b3   1.77  0.138   12.8   1.48   1.98   2   
#> 15 psi_b4  -1.32  0.107  -12.3  -1.52  -1.11  -1.5 
#> 16 psi_b5   2.84  0.138   20.6   2.63   3.13   3   
#> 17 psi_b6  -1.85  0.1000 -18.5  -2.03  -1.63  -2   
#> 18 psi_b7   0.914 0.0685  13.3   0.808  1.04   1   
#> 19 psi_b8   1.87  0.0923  20.3   1.71   2.07   2   
#> 20 scale    0.945 0.0386  24.5   0.874  1.01   1
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
#>  3.00e+05 simulations finished in 0.06 minutes. ( 86957 per second)
SummaryWelfare(wtp)
#> # A tibble: 2 x 5
#>   policy   Mean Std.Dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>   <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1     0       0           0            0
#> 2 policy2     0       0           0            0
```

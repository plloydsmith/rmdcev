
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
#> Chain 1: Initial log joint probability = -16971.5
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -13890.4      0.492817       79.2636           1           1       28   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -13834.7      0.223538       27.9628           1           1       50   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -13828.1      0.020021       21.5442      0.4965      0.4965       75   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -13827.4    0.00219684        1.5256      0.6438      0.6438       99   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       90      -13827.4    0.00013726      0.270762      0.3971           1      113   
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
#> LL                               : -13827.44
#> AIC                              : 27694.87
#> BIC                              : 27793.03
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.42
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 20.85 38.47 26.24 14.57 12.88 15.79  8.47 13.66 10.58 30.18 
#> 
#> 
#> Psi specification:
#> ~
#> b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 - 1
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -5.032   0.527  -9.55
#> psi_b2     0.642   0.103   6.25
#> psi_b3     2.004   0.108  18.55
#> psi_b4    -1.524   0.093 -16.41
#> psi_b5     3.023   0.148  20.44
#> psi_b6    -2.031   0.111 -18.34
#> psi_b7     1.106   0.093  11.95
#> psi_b8     2.025   0.106  19.16
#> gamma1     1.094   0.184   5.93
#> gamma2     1.321   0.184   7.18
#> gamma3     1.133   0.167   6.77
#> gamma4     1.662   0.293   5.68
#> gamma5     1.337   0.182   7.33
#> gamma6     1.308   0.239   5.47
#> gamma7     1.261   0.313   4.03
#> gamma8     2.554   0.463   5.52
#> gamma9     1.527   0.287   5.32
#> gamma10    2.083   0.295   7.06
#> alpha1     0.515   0.022  23.57
#> scale      0.997   0.042  23.86
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
#>  1 alpha1   0.515 0.0218  23.6   0.471  0.547  0.5 
#>  2 gamma1   1.09  0.184    5.93  0.864  1.48   1.04
#>  3 gamma10  2.08  0.295    7.06  1.71   2.85   1.96
#>  4 gamma2   1.32  0.184    7.18  0.989  1.69   1.57
#>  5 gamma3   1.13  0.167    6.77  0.854  1.49   1.34
#>  6 gamma4   1.66  0.293    5.68  1.22   2.15   1.84
#>  7 gamma5   1.34  0.182    7.33  1.04   1.66   1.39
#>  8 gamma6   1.31  0.239    5.47  0.853  1.75   1.28
#>  9 gamma7   1.26  0.313    4.03  0.813  1.95   1.62
#> 10 gamma8   2.55  0.463    5.52  1.88   3.50   1.79
#> 11 gamma9   1.53  0.287    5.32  1.12   2.22   1.66
#> 12 psi_b1  -5.03  0.527   -9.55 -6.08  -4.23  -5   
#> 13 psi_b2   0.642 0.103    6.25  0.460  0.864  0.5 
#> 14 psi_b3   2.00  0.108   18.5   1.82   2.21   2   
#> 15 psi_b4  -1.52  0.0928 -16.4  -1.68  -1.32  -1.5 
#> 16 psi_b5   3.02  0.148   20.4   2.80   3.33   3   
#> 17 psi_b6  -2.03  0.111  -18.3  -2.24  -1.85  -2   
#> 18 psi_b7   1.11  0.0926  11.9   0.946  1.25   1   
#> 19 psi_b8   2.02  0.106   19.2   1.85   2.24   2   
#> 20 scale    0.997 0.0418  23.9   0.937  1.08   1
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
#> 3.00e+05simulations finished in0.06minutes.(85227per second)
SummaryWelfare(wtp)
#> # A tibble: 2 x 5
#>   policy   Mean Std.Dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>   <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1     0       0           0            0
#> 2 policy2     0       0           0            0
```

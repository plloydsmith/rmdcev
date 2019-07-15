
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
#> Chain 1: Initial log joint probability = -16735.6
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -13725.1      0.278647       194.904           1           1       25   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -13620.1      0.447664       41.4844           1           1       46   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -13608.7      0.100036       46.5401           1           1       70   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -13607.6     0.0674252       4.29565           1           1       93   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -13607.4    0.00798931       1.31137           1           1      115   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      112      -13607.4   0.000482846      0.279241      0.5146      0.5146      131   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
SummaryMDCEV(mdcev_est)
#> Model run using rmdcev for R, version 1.0.0 
#> Estimation method                : MLE
#> Model type                       : hybrid specification
#> Number of classes                : 1
#> Number of individuals            : 1000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -13607.42
#> AIC                              : 27254.84
#> BIC                              : 27352.99
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.3
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 49.74 12.70 18.87 16.29  8.87 12.86 24.48 20.67  9.70 18.39 
#> 
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -5.265   0.547  -9.63
#> psi_b2     0.575   0.098   5.85
#> psi_b3     1.973   0.120  16.42
#> psi_b4    -1.558   0.124 -12.61
#> psi_b5     3.128   0.183  17.07
#> psi_b6    -2.050   0.151 -13.62
#> psi_b7     0.965   0.093  10.32
#> psi_b8     2.200   0.119  18.46
#> gamma1     1.675   0.275   6.09
#> gamma2     1.743   0.292   5.98
#> gamma3     0.943   0.144   6.57
#> gamma4     1.021   0.167   6.10
#> gamma5     1.306   0.242   5.41
#> gamma6     1.321   0.255   5.19
#> gamma7     1.558   0.233   6.70
#> gamma8     2.047   0.332   6.17
#> gamma9     0.888   0.165   5.40
#> gamma10    1.844   0.347   5.32
#> alpha1     0.502   0.028  17.70
#> scale      1.026   0.047  21.64
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
#>  1 alpha1   0.502 0.0284  17.7   0.451  0.556  0.5 
#>  2 gamma1   1.68  0.275    6.09  1.22   2.22   1.62
#>  3 gamma10  1.84  0.347    5.32  1.25   2.53   1.58
#>  4 gamma2   1.74  0.292    5.98  1.33   2.26   1.62
#>  5 gamma3   0.943 0.144    6.57  0.698  1.22   1.25
#>  6 gamma4   1.02  0.167    6.10  0.755  1.36   1.05
#>  7 gamma5   1.31  0.242    5.41  0.906  1.83   1.25
#>  8 gamma6   1.32  0.255    5.19  0.845  1.80   1.38
#>  9 gamma7   1.56  0.233    6.70  1.18   2.05   1.46
#> 10 gamma8   2.05  0.332    6.17  1.50   2.79   1.67
#> 11 gamma9   0.888 0.165    5.40  0.589  1.21   1.10
#> 12 psi_b1  -5.26  0.547   -9.63 -6.27  -4.32  -5   
#> 13 psi_b2   0.575 0.0982   5.85  0.397  0.754  0.5 
#> 14 psi_b3   1.97  0.120   16.4   1.74   2.18   2   
#> 15 psi_b4  -1.56  0.124  -12.6  -1.81  -1.29  -1.5 
#> 16 psi_b5   3.13  0.183   17.1   2.84   3.46   3   
#> 17 psi_b6  -2.05  0.151  -13.6  -2.27  -1.72  -2   
#> 18 psi_b7   0.965 0.0935  10.3   0.805  1.10   1   
#> 19 psi_b8   2.20  0.119   18.5   2.02   2.47   2   
#> 20 scale    1.03  0.0474  21.6   0.948  1.11   1
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
#> 3.00e+05simulations finished in0.06minutes.(81081per second)
SummaryWelfare(wtp)
#> # A tibble: 2 x 5
#>   policy   Mean Std.Dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>   <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1     0       0           0            0
#> 2 policy2     0       0           0            0
```

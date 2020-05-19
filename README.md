
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
#> Chain 1: Initial log joint probability = -62156.4
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -52830.5      0.991886       1907.31           1           1       23   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -52332.6     0.0546787       99.2564           1           1       45   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -52321.3   0.000845003       2.59614           1           1       67   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       71      -52321.3   0.000377997      0.663885           1           1       80   
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
#> LL                               : -52321.29
#> AIC                              : 104682.6
#> BIC                              : 104794.6
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.9
#> 
#> Average consumption of non-numeraire alternatives:
#>      1      2      3      4      5      6      7      8      9     10 
#>  11.68  32.93  27.80  25.61   7.96   4.66 132.25  86.05   5.34  19.15 
#> 
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -5.012   0.098 -51.14
#> psi_b2     0.466   0.049   9.50
#> psi_b3     1.991   0.060  33.19
#> psi_b4    -1.510   0.055 -27.45
#> psi_b5     3.039   0.048  63.30
#> psi_b6    -1.952   0.058 -33.65
#> psi_b7     0.982   0.041  23.95
#> psi_b8     1.990   0.042  47.38
#> gamma1     5.686   0.386  14.73
#> gamma2     6.886   0.373  18.46
#> gamma3     5.822   0.329  17.70
#> gamma4     7.217   0.431  16.75
#> gamma5     5.453   0.407  13.40
#> gamma6     3.648   0.284  12.85
#> gamma7     5.182   0.339  15.29
#> gamma8     8.426   0.442  19.06
#> gamma9     2.030   0.132  15.38
#> gamma10    4.289   0.241  17.80
#> alpha1     0.502   0.007  71.75
#> scale      0.988   0.012  82.36
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
#> psi_b1   psi_b1 -5.000000   -5.012   0.098 -51.14 -5.20408 -4.81992
#> psi_b2   psi_b2  0.500000    0.466   0.049   9.50  0.36996  0.56204
#> psi_b3   psi_b3  2.000000    1.991   0.060  33.19  1.87340  2.10860
#> psi_b4   psi_b4 -1.500000   -1.510   0.055 -27.45 -1.61780 -1.40220
#> psi_b5   psi_b5  3.000000    3.039   0.048  63.30  2.94492  3.13308
#> psi_b6   psi_b6 -2.000000   -1.952   0.058 -33.65 -2.06568 -1.83832
#> psi_b7   psi_b7  1.000000    0.982   0.041  23.95  0.90164  1.06236
#> psi_b8   psi_b8  2.000000    1.990   0.042  47.38  1.90768  2.07232
#> gamma1   gamma1  6.240237    5.686   0.386  14.73  4.92944  6.44256
#> gamma2   gamma2  7.291805    6.886   0.373  18.46  6.15492  7.61708
#> gamma3   gamma3  5.849051    5.822   0.329  17.70  5.17716  6.46684
#> gamma4   gamma4  7.488369    7.217   0.431  16.75  6.37224  8.06176
#> gamma5   gamma5  4.711628    5.453   0.407  13.40  4.65528  6.25072
#> gamma6   gamma6  3.709536    3.648   0.284  12.85  3.09136  4.20464
#> gamma7   gamma7  5.188255    5.182   0.339  15.29  4.51756  5.84644
#> gamma8   gamma8  7.724965    8.426   0.442  19.06  7.55968  9.29232
#> gamma9   gamma9  1.789074    2.030   0.132  15.38  1.77128  2.28872
#> gamma10 gamma10  4.284462    4.289   0.241  17.80  3.81664  4.76136
#> alpha1   alpha1  0.500000    0.502   0.007  71.75  0.48828  0.51572
#> scale     scale  1.000000    0.988   0.012  82.36  0.96448  1.01152
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
#> 6.00e+04simulations finished in0.44minutes.(2294per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy       mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#> * <chr>       <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 -5.76e-12      NA   -5.76e-12    -5.76e-12
#> 2 policy2 -5.76e-12      NA   -5.76e-12    -5.76e-12
```

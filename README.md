
# Kuhn-Tucker and Multiple Discrete-Continuous Extreme Value (MDCEV) Model Estimation and Simulation in R: The rmdcev Package

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/plloydsmith/rmdcev.svg?branch=master)](https://travis-ci.org/plloydsmith/rmdcev)
<!-- badges: end -->

The **rmdcev** R package estimate and simulates Kuhn-Tucker demand
models with individual heterogeneity. The models supported by **rmdcev**
are the multiple-discrete continuous extreme value (MDCEV) model and
Kuhn-Tucker specification common in the environmental economics
literature on recreation demand. Latent class and random parameters
specifications can be implemented and the models are fit using maximum
likelihood estimation or Bayesian estimation. All models are implemented
in Stan, which is a C++ package for performing full Bayesian inference
(see <http://mc-stan.org>). The **rmdcev** package also implements
demand forecasting and welfare calculation for policy simulation.

## Current Status

Development is in progress. Currently users can estimate the following
models:

1.  Bhat (2008) MDCEV model specifications
2.  Kuhn-Tucker model specification in environmental economics (von
    Haefen and Phaneuf, 2005)

Models can be estimated using

1.  Fixed parameter models (maximum likelihood or Bayesian estimation)
2.  Latent class models (maximum likelihood estimation)
3.  Random parameters models (Bayesian estimation)

## Installation

I recommend you first install **rstan** by following these steps:

<https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started>

Once **rstan** is installed, you can install the released version of
**rmdcev** from GitHub using devtools

``` r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("plloydsmith/rmdcev", build_vignettes = FALSE, INSTALL_opts="--no-multiarch")
```

If you have any issues with installation or use of the package, please
let me know by [filing an
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

For more details on the demand and welfare simulation:

Pinjari, A.R. and Bhat , C.R. (2011) [“Computationally Efficient
Forecasting Procedures for Kuhn-Tucker Consumer Demand Model Systems:
Application to Residential Energy Consumption
Analysis.”](https://www.caee.utexas.edu/prof/bhat/MDCEV_Forecasting.html)
Technical paper, Department of Civil & Environmental Engineering,
University of South Florida.

Lloyd-Smith, P (2018). [“A New Approach to Calculating Welfare Measures
in Kuhn-Tucker Demand
Models.”](https://www.sciencedirect.com/science/article/pii/S1755534517300994)
Journal of Choice Modeling, 26: 19-27

## Estimation

As an example, we can simulate some data using Bhat (2008)‘s ’Gamma’
specification. In this example, we are simulating data for 2,000
individuals and 10 non-numeraire alternatives. We will randomly generate
the parameter values to simulate the data and then check these values to
our estimation results.

``` r
library(pacman)
p_load(tidyverse, rmdcev)
model <- "gamma"
nobs <- 2000
nalts <- 10
sim.data <- GenerateMDCEVData(model = model, nobs = nobs, nalts = nalts)
#> Sorting data by id.var then alt...
#> Checking data...
#> Data is good
```

Estimate model using MLE (note that we set “psi\_ascs = 0” to omit any
alternative-specific constants)

``` r
mdcev_est <- mdcev(~ b1 + b2 + b3 + b4 + b5 + b6,
                   data = sim.data$data,
                   psi_ascs = 0,
                   model = model,
                   algorithm = "MLE")
#> Using MLE to estimate KT model
#> Chain 1: Initial log joint probability = -91770.7
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1: Error evaluating model log probability: Non-finite gradient.
#> Error evaluating model log probability: Non-finite gradient.
#> 
#> Chain 1:       19      -36100.6      0.232987       233.736           1           1       30   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -35881.7      0.132855       181.545           1           1       52   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -35866.5     0.0184806       4.88176           1           1       75   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -35866.4     0.0010381       1.19302           1           1      101   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       81      -35866.4   0.000987269      0.677771           1           1      103   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.2.0 
#> Estimation method                : MLE
#> Model type                       : gamma specification
#> Number of classes                : 1
#> Number of individuals            : 2000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 18
#> LL                               : -35866.41
#> AIC                              : 71768.83
#> BIC                              : 71869.64
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:2.75
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 21.05 13.39  1.45 16.38 21.13  3.74 14.54 11.85 16.69  2.54 
#> 
#> Parameter estimates --------------------------------  
#>           Estimate Std.err z.stat
#> psi_b1      -5.082   0.122 -41.66
#> psi_b2       0.659   0.074   8.91
#> psi_b3       2.042   0.065  31.41
#> psi_b4      -1.438   0.053 -27.13
#> psi_b5       1.975   0.047  42.02
#> psi_b6      -0.897   0.049 -18.31
#> gamma_1      3.681   0.216  17.04
#> gamma_2      3.976   0.267  14.89
#> gamma_3      6.467   1.206   5.36
#> gamma_4      9.246   0.652  14.18
#> gamma_5      2.008   0.117  17.16
#> gamma_6      8.791   1.061   8.29
#> gamma_7      8.215   0.571  14.39
#> gamma_8      2.875   0.181  15.89
#> gamma_9      2.963   0.169  17.53
#> gamma_10     8.937   1.319   6.78
#> alpha_num    0.518   0.008  64.74
#> scale        0.986   0.015  65.75
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
#>      parms      true Estimate Std.err z.stat    cl_lo    cl_hi
#> 1   psi_b1 -5.000000   -5.082   0.122 -41.66 -5.32112 -4.84288
#> 2   psi_b2  0.500000    0.659   0.074   8.91  0.51396  0.80404
#> 3   psi_b3  2.000000    2.042   0.065  31.41  1.91460  2.16940
#> 4   psi_b4 -1.500000   -1.438   0.053 -27.13 -1.54188 -1.33412
#> 5   psi_b5  2.000000    1.975   0.047  42.02  1.88288  2.06712
#> 6   psi_b6 -1.000000   -0.897   0.049 -18.31 -0.99304 -0.80096
#> 7   gamma1  3.823750    3.681   0.216  17.04  3.25764  4.10436
#> 8   gamma2  3.524959    3.976   0.267  14.89  3.45268  4.49932
#> 9   gamma3  6.022275    6.467   1.206   5.36  4.10324  8.83076
#> 10  gamma4  9.488652    9.246   0.652  14.18  7.96808 10.52392
#> 11  gamma5  1.878443    2.008   0.117  17.16  1.77868  2.23732
#> 12  gamma6  8.347107    8.791   1.061   8.29  6.71144 10.87056
#> 13  gamma7  9.217354    8.215   0.571  14.39  7.09584  9.33416
#> 14  gamma8  2.597435    2.875   0.181  15.89  2.52024  3.22976
#> 15  gamma9  2.613302    2.963   0.169  17.53  2.63176  3.29424
#> 16 gamma10  8.459189    8.937   1.319   6.78  6.35176 11.52224
#> 17  alpha1  0.500000    0.518   0.008  64.74  0.50232  0.53368
#> 18   scale  1.000000    0.986   0.015  65.75  0.95660  1.01540
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

policies <- CreateBlankPolicies(npols, mdcev_est)

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
#> Using general approach in simulation...
#> 
#> 6.00e+04simulations finished in0.41minutes.(2414per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy      mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>      <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 4.07e-11      NA    4.07e-11     4.07e-11
#> 2 policy2 4.07e-11      NA    4.07e-11     4.07e-11
```

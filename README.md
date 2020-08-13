
# Kuhn-Tucker and Multiple Discrete-Continuous Extreme Value (MDCEV) Model Estimation and Simulation in R: The rmdcev Package

<!-- badges: start -->

[![R build
status](https://github.com/plloydsmith/rmdcev/workflows/R-CMD-check/badge.svg)](https://github.com/plloydsmith/rmdcev/actions)
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
#> Chain 1: Initial log joint probability = -86260
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1: Error evaluating model log probability: Non-finite gradient.
#> Error evaluating model log probability: Non-finite gradient.
#> 
#> Chain 1:       19      -34176.3      0.928606       434.214           1           1       31   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -33933.1       0.16448       51.5043           1           1       51   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -33925.2     0.0267053       15.2788           1           1       77   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -33924.9    0.00121969      0.941475           1           1       99   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       92      -33924.9   0.000774317      0.507544      0.9556      0.9556      113   
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
#> LL                               : -33924.91
#> AIC                              : 67885.82
#> BIC                              : 67986.64
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:2.81
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#>  0.53 16.33 14.19  2.26 69.68  0.58  1.36 15.44 13.93 16.86 
#> 
#> Parameter estimates --------------------------------  
#>           Estimate Std.err z.stat
#> psi_b1      -5.382   0.168 -32.03
#> psi_b2       0.364   0.099   3.68
#> psi_b3       2.088   0.060  34.80
#> psi_b4      -1.553   0.059 -26.32
#> psi_b5       2.123   0.052  40.82
#> psi_b6      -1.074   0.055 -19.53
#> gamma_1      5.993   1.399   4.28
#> gamma_2      6.580   0.460  14.31
#> gamma_3      3.814   0.259  14.73
#> gamma_4      7.681   1.052   7.30
#> gamma_5      4.556   0.274  16.63
#> gamma_6      2.493   0.425   5.87
#> gamma_7      7.059   1.173   6.02
#> gamma_8      3.265   0.208  15.70
#> gamma_9      3.952   0.260  15.20
#> gamma_10     1.767   0.102  17.32
#> alpha_num    0.489   0.011  44.43
#> scale        1.033   0.016  64.59
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
#> 1   psi_b1 -5.000000   -5.382   0.168 -32.03 -5.71128 -5.05272
#> 2   psi_b2  0.500000    0.364   0.099   3.68  0.16996  0.55804
#> 3   psi_b3  2.000000    2.088   0.060  34.80  1.97040  2.20560
#> 4   psi_b4 -1.500000   -1.553   0.059 -26.32 -1.66864 -1.43736
#> 5   psi_b5  2.000000    2.123   0.052  40.82  2.02108  2.22492
#> 6   psi_b6 -1.000000   -1.074   0.055 -19.53 -1.18180 -0.96620
#> 7   gamma1  7.478434    5.993   1.399   4.28  3.25096  8.73504
#> 8   gamma2  7.460036    6.580   0.460  14.31  5.67840  7.48160
#> 9   gamma3  4.381761    3.814   0.259  14.73  3.30636  4.32164
#> 10  gamma4  7.720269    7.681   1.052   7.30  5.61908  9.74292
#> 11  gamma5  4.819395    4.556   0.274  16.63  4.01896  5.09304
#> 12  gamma6  2.431377    2.493   0.425   5.87  1.66000  3.32600
#> 13  gamma7  8.083688    7.059   1.173   6.02  4.75992  9.35808
#> 14  gamma8  3.240990    3.265   0.208  15.70  2.85732  3.67268
#> 15  gamma9  4.032100    3.952   0.260  15.20  3.44240  4.46160
#> 16 gamma10  1.962492    1.767   0.102  17.32  1.56708  1.96692
#> 17  alpha1  0.500000    0.489   0.011  44.43  0.46744  0.51056
#> 18   scale  1.000000    1.033   0.016  64.59  1.00164  1.06436
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
#> 6.00e+04simulations finished in0.39minutes.(2562per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy       mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>       <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 -7.91e-11      NA   -7.91e-11    -7.91e-11
#> 2 policy2 -7.91e-11      NA   -7.91e-11    -7.91e-11
```

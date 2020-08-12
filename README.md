
# Kuhn-Tucker and Multiple Discrete-Continuous Extreme Value (MDCEV) Model Estimation and Simulation in R: The rmdcev Package

<!-- badges: start -->

[![Build
Status](https://travis-ci.org/plloydsmith/rmdcev.svg?branch=master)](https://travis-ci.org/plloydsmith/rmdcev)
<!-- badges: end -->

The *rmdcev* R package estimate and simulates Kuhn-Tucker demand models
with individual heterogeneity. The models supported by *rmdcev* are the
multiple-discrete continuous extreme value (MDCEV) model and Kuhn-Tucker
specification common in the environmental economics literature on
recreation demand. Latent class and random parameters specifications can
be implemented and the models are fit using maximum likelihood
estimation or Bayesian estimation. All models are implemented in Stan,
which is a C++ package for performing full Bayesian inference (see
<http://mc-stan.org/>). The *rmdcev* package also implements demand
forecasting and welfare calculation for policy simulation.

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
rmdcev from GitHub using devtools

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
#> Chain 1: Initial log joint probability = -72987.3
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1: Error evaluating model log probability: Non-finite gradient.
#> 
#> Chain 1:       19      -23345.5      0.156605       222.227           1           1       31   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -23251.6     0.0940246       80.7835           1           1       55   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -23238.4     0.0436308       64.8183      0.2972           1       77   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -23234.8      0.014661       4.36058           1           1       98   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -23234.6    0.00315425       2.96504      0.8754      0.8754      119   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      119      -23234.5    0.00094661      0.475567           1           1      142   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      128      -23234.5   0.000268274      0.287345           1           1      151   
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
#> LL                               : -23234.55
#> AIC                              : 46505.1
#> BIC                              : 46605.91
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:3.44
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#>  7.23  1.73  0.12 58.74  1.63 16.13  6.38  1.94  0.10  3.93 
#> 
#> Parameter estimates --------------------------------  
#>           Estimate Std.err z.stat
#> psi_b1      -4.770   0.110 -43.37
#> psi_b2       0.600   0.063   9.52
#> psi_b3       2.127   0.089  23.90
#> psi_b4      -1.487   0.060 -24.78
#> psi_b5       1.850   0.053  34.90
#> psi_b6      -0.951   0.054 -17.61
#> gamma_1     12.042   1.299   9.27
#> gamma_2      2.905   0.327   8.88
#> gamma_3      4.154   1.375   3.02
#> gamma_4      4.411   0.305  14.46
#> gamma_5      3.700   0.444   8.33
#> gamma_6      2.107   0.139  15.16
#> gamma_7      4.290   0.363  11.82
#> gamma_8      3.953   0.494   8.00
#> gamma_9      1.881   0.530   3.55
#> gamma_10     7.642   0.915   8.35
#> alpha_num    0.504   0.009  55.99
#> scale        0.970   0.018  53.90
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
#> 1   psi_b1 -5.000000   -4.770   0.110 -43.37 -4.98560 -4.55440
#> 2   psi_b2  0.500000    0.600   0.063   9.52  0.47652  0.72348
#> 3   psi_b3  2.000000    2.127   0.089  23.90  1.95256  2.30144
#> 4   psi_b4 -1.500000   -1.487   0.060 -24.78 -1.60460 -1.36940
#> 5   psi_b5  2.000000    1.850   0.053  34.90  1.74612  1.95388
#> 6   psi_b6 -1.000000   -0.951   0.054 -17.61 -1.05684 -0.84516
#> 7   gamma1  9.571871   12.042   1.299   9.27  9.49596 14.58804
#> 8   gamma2  2.662405    2.905   0.327   8.88  2.26408  3.54592
#> 9   gamma3  5.278969    4.154   1.375   3.02  1.45900  6.84900
#> 10  gamma4  4.149304    4.411   0.305  14.46  3.81320  5.00880
#> 11  gamma5  3.880138    3.700   0.444   8.33  2.82976  4.57024
#> 12  gamma6  2.032531    2.107   0.139  15.16  1.83456  2.37944
#> 13  gamma7  4.464965    4.290   0.363  11.82  3.57852  5.00148
#> 14  gamma8  4.318229    3.953   0.494   8.00  2.98476  4.92124
#> 15  gamma9  1.439274    1.881   0.530   3.55  0.84220  2.91980
#> 16 gamma10  7.136159    7.642   0.915   8.35  5.84860  9.43540
#> 17  alpha1  0.500000    0.504   0.009  55.99  0.48636  0.52164
#> 18   scale  1.000000    0.970   0.018  53.90  0.93472  1.00528
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
#> 6.00e+04simulations finished in0.34minutes.(2970per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy       mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>       <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 -1.65e-10      NA   -1.65e-10    -1.65e-10
#> 2 policy2 -1.65e-10      NA   -1.65e-10    -1.65e-10
```

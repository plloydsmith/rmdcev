

# Kuhn-Tucker and Multiple Discrete-Continuous Extreme Value (MDCEV) Model Estimation and Simulation in R: The rmdcev Package

<!-- badges: start -->

[![R build
status](https://github.com/plloydsmith/rmdcev/workflows/R-CMD-check/badge.svg)](https://github.com/plloydsmith/rmdcev/actions)
[![R-CMD-check](https://github.com/plloydsmith/rmdcev/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/plloydsmith/rmdcev/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The **rmdcev** R package estimate and simulates Kuhn-Tucker demand
models with individual heterogeneity. The models supported by **rmdcev**
are the multiple-discrete continuous extreme value (MDCEV) model and
Kuhn-Tucker specification common in the environmental economics
literature on recreation demand. Latent class and random parameters
specifications can be implemented and the models are fit using maximum
likelihood estimation or Bayesian estimation. All models are implemented
in Stan, which is a C++ package for performing full Bayesian inference
(see https://mc-stan.org/). The **rmdcev** package also implements
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

https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

Once **rstan** is installed, you can install **rmdcev** from CRAN using

``` r
install.packages("rmdcev")
```

Or install the latest version of **rmdcev** from GitHub using devtools

``` r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("plloydsmith/rmdcev", dependencies = TRUE, INSTALL_opts="--no-multiarch")
```

If you have any issues with installation or use of the package, please
let me know by [filing an
issue](https://github.com/plloydsmith/rmdcev/issues).

## References

Background on the models, estimation, and simulation as well as a walk
through of the package is provided in

Lloyd-Smith, P (2021). [“Kuhn-Tucker and Multiple Discrete-Continuous
Extreme Value Model Estimation and Simulation in R: The rmdcev
Package”](https://doi.org/10.32614/RJ-2021-015) The R Journal, 12(2):
251-265.

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
set.seed(12345)
model <- "gamma"
nobs <- 2000
nalts <- 10
sim.data <- GenerateMDCEVData(model = model, nobs = nobs, nalts = nalts)
#> Sorting data by id.var then alt...
#> Checking data...
#> Data is good
```

Estimate model using MLE (note that we set “psi_ascs = 0” to omit any
alternative-specific constants)

``` r
mdcev_est <- mdcev(~ b1 + b2 + b3 + b4 + b5 + b6,
                   data = sim.data$data,
                   psi_ascs = 0,
                   model = model,
                   algorithm = "MLE")
#> Using MLE to estimate KT model
#> Initial log joint probability = -47202 
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>       99      -35723.4    0.00529771       3.53273      0.3392      0.9583      118    
#>     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes  
#>      112      -35723.4     0.0004941       0.47622           1           1      132    
#> Optimization terminated normally:  
#>   Convergence detected: relative gradient magnitude is below tolerance 
#> Finished in  0.9 seconds.
```

Summarize results

``` r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.3.0 
#> Estimation method                : MLE
#> Model type                       : gamma specification
#> Number of classes                : 1
#> Number of individuals            : 2000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 18
#> LL                               : -35723.35
#> AIC                              : 71482.7
#> BIC                              : 71583.52
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:2.57
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 59.90 10.43  0.84 71.64  4.89  1.46 10.40 13.49 21.86  0.37 
#> 
#> Parameter estimates --------------------------------  
#>           Estimate Std.err z.stat
#> psi_b1      -4.894   0.115 -42.55
#> psi_b2       0.555   0.091   6.09
#> psi_b3       2.010   0.062  32.42
#> psi_b4      -1.500   0.057 -26.32
#> psi_b5       2.078   0.046  45.18
#> psi_b6      -1.088   0.055 -19.78
#> gamma_1      6.982   0.411  16.99
#> gamma_2      8.449   0.741  11.40
#> gamma_3      7.349   1.521   4.83
#> gamma_4      8.743   0.535  16.34
#> gamma_5      4.881   0.425  11.48
#> gamma_6      2.147   0.234   9.17
#> gamma_7      3.449   0.232  14.87
#> gamma_8      5.589   0.385  14.52
#> gamma_9      7.676   0.509  15.08
#> gamma_10     7.797   2.748   2.84
#> alpha_num    0.503   0.008  62.86
#> scale        1.000   0.015  66.67
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
#>             parms      true Estimate Std.err z.stat    cl_lo    cl_hi
#> psi_b1     psi_b1 -5.000000   -4.894   0.115 -42.55 -5.11940 -4.66860
#> psi_b2     psi_b2  0.500000    0.555   0.091   6.09  0.37664  0.73336
#> psi_b3     psi_b3  2.000000    2.010   0.062  32.42  1.88848  2.13152
#> psi_b4     psi_b4 -1.500000   -1.500   0.057 -26.32 -1.61172 -1.38828
#> psi_b5     psi_b5  2.000000    2.078   0.046  45.18  1.98784  2.16816
#> psi_b6     psi_b6 -1.000000   -1.088   0.055 -19.78 -1.19580 -0.98020
#> gamma_1    gamma1  7.488135    6.982   0.411  16.99  6.17644  7.78756
#> gamma_2    gamma2  8.881959    8.449   0.741  11.40  6.99664  9.90136
#> gamma_3    gamma3  7.848841    7.349   1.521   4.83  4.36784 10.33016
#> gamma_4    gamma4  8.975121    8.743   0.535  16.34  7.69440  9.79160
#> gamma_5    gamma5  5.108329    4.881   0.425  11.48  4.04800  5.71400
#> gamma_6    gamma6  2.497346    2.147   0.234   9.17  1.68836  2.60564
#> gamma_7    gamma7  3.925858    3.449   0.232  14.87  2.99428  3.90372
#> gamma_8    gamma8  5.583019    5.589   0.385  14.52  4.83440  6.34360
#> gamma_9    gamma9  7.549347    7.676   0.509  15.08  6.67836  8.67364
#> gamma_10  gamma10  9.907632    7.797   2.748   2.84  2.41092 13.18308
#> alpha_num  alpha1  0.500000    0.503   0.008  62.86  0.48732  0.51868
#> scale       scale  1.000000    1.000   0.015  66.67  0.97060  1.02940
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

<img src="man/figures/README-unnamed-chunk-6-1.png"
style="width:100.0%" />

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
#> 6.00e+04simulations finished in0.19minutes.(5155per second)
summary(wtp)
#> # A tibble: 2 × 5
#>   policy      mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>      <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 1.51e-10      NA    1.51e-10     1.51e-10
#> 2 policy2 1.51e-10      NA    1.51e-10     1.51e-10
```

## Thanks

This package was not developed in isolation and I gratefully acknowledge
Joshua Abbott, Allen Klaiber, Lusi Xie, the [apollo
team](http://www.apollochoicemodelling.com/), and the Stan team, whose
codes or suggestions were helpful in putting this package together.

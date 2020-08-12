
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
3.  Random parameters models (Bayesian estimation)

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
#> Chain 1: Initial log joint probability = -77850.8
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1: Error evaluating model log probability: Non-finite gradient.
#> Error evaluating model log probability: Non-finite gradient.
#> 
#> Chain 1:       19      -29151.4      0.387104       268.866           1           1       30   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -29074.5     0.0456729       33.6302           1           1       54   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -29064.8      0.017239        26.829      0.8147      0.8147       75   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -29064.3     0.0067155       2.44497           1           1       97   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       92      -29064.3   0.000473185      0.375743           1           1      114   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.1.3 
#> Estimation method                : MLE
#> Model type                       : gamma specification
#> Number of classes                : 1
#> Number of individuals            : 2000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 18
#> LL                               : -29064.31
#> AIC                              : 58164.62
#> BIC                              : 58265.44
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:2.7
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 22.36 17.37  6.77  1.94  6.16 20.07  0.26 19.86  1.57  2.35 
#> 
#> Parameter estimates --------------------------------  
#>           Estimate Std.err z.stat
#> psi_b1      -4.884   0.109 -44.81
#> psi_b2       0.413   0.070   5.90
#> psi_b3       1.760   0.093  18.92
#> psi_b4      -1.493   0.058 -25.75
#> psi_b5       2.039   0.054  37.75
#> psi_b6      -0.985   0.053 -18.58
#> gamma_1      7.574   0.511  14.82
#> gamma_2      2.192   0.139  15.77
#> gamma_3      4.694   0.391  12.00
#> gamma_4      3.531   0.390   9.05
#> gamma_5      3.584   0.284  12.62
#> gamma_6      2.493   0.170  14.67
#> gamma_7      2.986   0.813   3.67
#> gamma_8      4.449   0.301  14.78
#> gamma_9      4.241   0.564   7.52
#> gamma_10     7.839   1.107   7.08
#> alpha_num    0.491   0.009  54.53
#> scale        1.003   0.017  58.98
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
#> 1   psi_b1 -5.000000   -4.884   0.109 -44.81 -5.09764 -4.67036
#> 2   psi_b2  0.500000    0.413   0.070   5.90  0.27580  0.55020
#> 3   psi_b3  2.000000    1.760   0.093  18.92  1.57772  1.94228
#> 4   psi_b4 -1.500000   -1.493   0.058 -25.75 -1.60668 -1.37932
#> 5   psi_b5  2.000000    2.039   0.054  37.75  1.93316  2.14484
#> 6   psi_b6 -1.000000   -0.985   0.053 -18.58 -1.08888 -0.88112
#> 7   gamma1  7.314567    7.574   0.511  14.82  6.57244  8.57556
#> 8   gamma2  2.097729    2.192   0.139  15.77  1.91956  2.46444
#> 9   gamma3  5.392488    4.694   0.391  12.00  3.92764  5.46036
#> 10  gamma4  3.970054    3.531   0.390   9.05  2.76660  4.29540
#> 11  gamma5  3.523847    3.584   0.284  12.62  3.02736  4.14064
#> 12  gamma6  2.684004    2.493   0.170  14.67  2.15980  2.82620
#> 13  gamma7  2.787718    2.986   0.813   3.67  1.39252  4.57948
#> 14  gamma8  4.648740    4.449   0.301  14.78  3.85904  5.03896
#> 15  gamma9  5.393200    4.241   0.564   7.52  3.13556  5.34644
#> 16 gamma10  8.095358    7.839   1.107   7.08  5.66928 10.00872
#> 17  alpha1  0.500000    0.491   0.009  54.53  0.47336  0.50864
#> 18   scale  1.000000    1.003   0.017  58.98  0.96968  1.03632
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
#> 6.00e+04simulations finished in0.36minutes.(2783per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy       mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>       <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 -3.86e-12      NA   -3.86e-12    -3.86e-12
#> 2 policy2 -3.86e-12      NA   -3.86e-12    -3.86e-12
```

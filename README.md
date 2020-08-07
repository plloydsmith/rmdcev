
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
individuals and 10 non-numeraire alternatives. We will use randomly
generate the parameter values to estimate and then check our estimation
results.

``` r
library(pacman)
p_load(tidyverse, rmdcev)
model <- "gamma"
nobs <- 10000
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
#> Chain 1: Initial log joint probability = -470237
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1: Error evaluating model log probability: Non-finite gradient.
#> Error evaluating model log probability: Non-finite gradient.
#> 
#> Chain 1:       19       -155997     0.0626126       2703.73       0.248       0.248       32   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39       -154683      0.155554        465.91           1           1       54   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59       -154652      0.020183       48.4093      0.9566      0.9566       75   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79       -154651    0.00202988       8.74596           1           1       99   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       93       -154651   0.000157087      0.664501      0.2613      0.9589      118   
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
#> Number of individuals            : 10000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 18
#> LL                               : -154651.2
#> AIC                              : 309338.4
#> BIC                              : 309468.2
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:17.19
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#>  0.54  1.98 21.62  4.26 44.29  1.48 25.87  9.11 32.02  0.60 
#> 
#> Parameter estimates --------------------------------  
#>          Estimate Std.err z.stat
#> psi_b1     -4.923   0.059 -83.43
#> psi_b2      0.509   0.030  16.97
#> psi_b3      2.144   0.054  39.71
#> psi_b4     -1.508   0.026 -57.99
#> psi_b5      1.966   0.024  81.93
#> psi_b6     -1.012   0.024 -42.18
#> gamma_1     3.623   0.318  11.39
#> gamma_2     7.407   0.532  13.92
#> gamma_3     1.637   0.044  37.21
#> gamma_4     8.749   0.450  19.44
#> gamma_5     5.938   0.166  35.77
#> gamma_6     4.988   0.348  14.33
#> gamma_7     9.460   0.280  33.79
#> gamma_8     5.674   0.200  28.37
#> gamma_9     8.487   0.245  34.64
#> gamma_10    2.707   0.206  13.14
#> alpha_1     0.506   0.005 101.25
#> scale       1.000   0.007 142.90
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
#> 1   psi_b1 -5.000000   -4.923   0.059 -83.43 -5.03864 -4.80736
#> 2   psi_b2  0.500000    0.509   0.030  16.97  0.45020  0.56780
#> 3   psi_b3  2.000000    2.144   0.054  39.71  2.03816  2.24984
#> 4   psi_b4 -1.500000   -1.508   0.026 -57.99 -1.55896 -1.45704
#> 5   psi_b5  2.000000    1.966   0.024  81.93  1.91896  2.01304
#> 6   psi_b6 -1.000000   -1.012   0.024 -42.18 -1.05904 -0.96496
#> 7   gamma1  3.790078    3.623   0.318  11.39  2.99972  4.24628
#> 8   gamma2  7.156351    7.407   0.532  13.92  6.36428  8.44972
#> 9   gamma3  1.650154    1.637   0.044  37.21  1.55076  1.72324
#> 10  gamma4  9.089130    8.749   0.450  19.44  7.86700  9.63100
#> 11  gamma5  5.855912    5.938   0.166  35.77  5.61264  6.26336
#> 12  gamma6  4.372138    4.988   0.348  14.33  4.30592  5.67008
#> 13  gamma7  9.572931    9.460   0.280  33.79  8.91120 10.00880
#> 14  gamma8  5.798455    5.674   0.200  28.37  5.28200  6.06600
#> 15  gamma9  8.031057    8.487   0.245  34.64  8.00680  8.96720
#> 16 gamma10  2.088844    2.707   0.206  13.14  2.30324  3.11076
#> 17  alpha1  0.500000    0.506   0.005 101.25  0.49620  0.51580
#> 18   scale  1.000000    1.000   0.007 142.90  0.98628  1.01372
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
#> Using general approach in simulation...
#> 
#> 3.00e+05simulations finished in1.49minutes.(3366per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy      mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>      <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 9.22e-12      NA    9.22e-12     9.22e-12
#> 2 policy2 9.22e-12      NA    9.22e-12     9.22e-12
```

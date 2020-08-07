
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
#> Chain 1: Initial log joint probability = -60303.6
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1: Error evaluating model log probability: Non-finite gradient.
#> Error evaluating model log probability: Non-finite gradient.
#> 
#> Chain 1:       19      -22587.6      0.147465        239.08           1           1       29   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -22417.9      0.105684       88.3198           1           1       52   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -22414.6    0.00557064        6.4863      0.4163           1       75   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       77      -22414.6   0.000197231      0.596081       0.341      0.9852      100   
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
#> Estimated parameters             : 18
#> LL                               : -22414.6
#> AIC                              : 44865.21
#> BIC                              : 44966.02
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:2.58
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#>  1.11  1.65  9.71  1.06 14.70  5.16  1.31  5.32  5.64  9.06 
#> 
#> Parameter estimates --------------------------------  
#>          Estimate Std.err z.stat
#> psi_b1     -5.017   0.129 -38.90
#> psi_b2      0.585   0.080   7.31
#> psi_b3      2.067   0.071  29.11
#> psi_b4     -1.496   0.059 -25.36
#> psi_b5      1.921   0.059  32.56
#> psi_b6     -0.988   0.056 -17.64
#> gamma_1     1.890   0.225   8.40
#> gamma_2     3.831   0.586   6.54
#> gamma_3     1.074   0.074  14.51
#> gamma_4     1.506   0.168   8.96
#> gamma_5     3.137   0.218  14.39
#> gamma_6     7.937   0.838   9.47
#> gamma_7     7.234   1.453   4.98
#> gamma_8     1.314   0.093  14.13
#> gamma_9     1.676   0.123  13.62
#> gamma_10   10.196   0.995  10.25
#> alpha_1     0.499   0.010  49.89
#> scale       0.994   0.019  52.33
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
#> 1   psi_b1 -5.000000   -5.017   0.129 -38.90 -5.26984 -4.76416
#> 2   psi_b2  0.500000    0.585   0.080   7.31  0.42820  0.74180
#> 3   psi_b3  2.000000    2.067   0.071  29.11  1.92784  2.20616
#> 4   psi_b4 -1.500000   -1.496   0.059 -25.36 -1.61164 -1.38036
#> 5   psi_b5  2.000000    1.921   0.059  32.56  1.80536  2.03664
#> 6   psi_b6 -1.000000   -0.988   0.056 -17.64 -1.09776 -0.87824
#> 7   gamma1  1.672031    1.890   0.225   8.40  1.44900  2.33100
#> 8   gamma2  4.454715    3.831   0.586   6.54  2.68244  4.97956
#> 9   gamma3  1.121848    1.074   0.074  14.51  0.92896  1.21904
#> 10  gamma4  1.779145    1.506   0.168   8.96  1.17672  1.83528
#> 11  gamma5  2.666243    3.137   0.218  14.39  2.70972  3.56428
#> 12  gamma6  7.543130    7.937   0.838   9.47  6.29452  9.57948
#> 13  gamma7  8.978440    7.234   1.453   4.98  4.38612 10.08188
#> 14  gamma8  1.198326    1.314   0.093  14.13  1.13172  1.49628
#> 15  gamma9  1.641405    1.676   0.123  13.62  1.43492  1.91708
#> 16 gamma10  9.692787   10.196   0.995  10.25  8.24580 12.14620
#> 17  alpha1  0.500000    0.499   0.010  49.89  0.47940  0.51860
#> 18   scale  1.000000    0.994   0.019  52.33  0.95676  1.03124
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
#> 6.00e+04simulations finished in0.28minutes.(3636per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy       mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>       <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 -1.77e-10      NA   -1.77e-10    -1.77e-10
#> 2 policy2 -1.77e-10      NA   -1.77e-10    -1.77e-10
```

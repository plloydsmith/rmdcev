
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
(see <http://mc-stan.org/>). The **rmdcev** package also implements
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
#> Chain 1: Initial log joint probability = -102091
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1: Error evaluating model log probability: Non-finite gradient.
#> Error evaluating model log probability: Non-finite gradient.
#> 
#> Chain 1:       19      -36191.9      0.238697        170.87           1           1       33   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -36099.7       0.09316       104.396           1           1       54   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -36075.1     0.0173156       53.9469      0.3078      0.3078       80   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -36072.3     0.0278811       12.2672           1           1      103   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -36072.2    0.00180446      0.959796           1           1      124   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      107      -36072.2   0.000747187      0.522778      0.9515      0.9515      135   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

``` r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.2.1 
#> Estimation method                : MLE
#> Model type                       : gamma specification
#> Number of classes                : 1
#> Number of individuals            : 2000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 18
#> LL                               : -36072.19
#> AIC                              : 72180.38
#> BIC                              : 72281.2
#> Standard errors calculated using : Delta method
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:2.85
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#>  7.18  0.95 16.85 11.57  8.81 58.15 18.51 17.74  0.30  6.56 
#> 
#> Parameter estimates --------------------------------  
#>           Estimate Std.err z.stat
#> psi_b1      -5.012   0.126 -39.78
#> psi_b2       0.484   0.082   5.91
#> psi_b3       2.107   0.074  28.47
#> psi_b4      -1.513   0.055 -27.52
#> psi_b5       2.015   0.047  42.86
#> psi_b6      -1.027   0.053 -19.38
#> gamma_1      5.673   0.480  11.82
#> gamma_2      1.799   0.207   8.69
#> gamma_3      2.237   0.140  15.98
#> gamma_4      2.440   0.147  16.60
#> gamma_5      6.482   0.527  12.30
#> gamma_6      8.840   0.512  17.27
#> gamma_7      2.684   0.156  17.20
#> gamma_8      8.299   0.573  14.48
#> gamma_9      6.507   1.892   3.44
#> gamma_10     3.867   0.313  12.36
#> alpha_num    0.501   0.008  62.67
#> scale        1.006   0.015  67.04
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
#> 1   psi_b1 -5.000000   -5.012   0.126 -39.78 -5.25896 -4.76504
#> 2   psi_b2  0.500000    0.484   0.082   5.91  0.32328  0.64472
#> 3   psi_b3  2.000000    2.107   0.074  28.47  1.96196  2.25204
#> 4   psi_b4 -1.500000   -1.513   0.055 -27.52 -1.62080 -1.40520
#> 5   psi_b5  2.000000    2.015   0.047  42.86  1.92288  2.10712
#> 6   psi_b6 -1.000000   -1.027   0.053 -19.38 -1.13088 -0.92312
#> 7   gamma1  5.880124    5.673   0.480  11.82  4.73220  6.61380
#> 8   gamma2  2.037718    1.799   0.207   8.69  1.39328  2.20472
#> 9   gamma3  2.326599    2.237   0.140  15.98  1.96260  2.51140
#> 10  gamma4  2.378944    2.440   0.147  16.60  2.15188  2.72812
#> 11  gamma5  6.188133    6.482   0.527  12.30  5.44908  7.51492
#> 12  gamma6  9.298275    8.840   0.512  17.27  7.83648  9.84352
#> 13  gamma7  2.845774    2.684   0.156  17.20  2.37824  2.98976
#> 14  gamma8  7.553436    8.299   0.573  14.48  7.17592  9.42208
#> 15  gamma9  8.992650    6.507   1.892   3.44  2.79868 10.21532
#> 16 gamma10  3.494510    3.867   0.313  12.36  3.25352  4.48048
#> 17  alpha1  0.500000    0.501   0.008  62.67  0.48532  0.51668
#> 18   scale  1.000000    1.006   0.015  67.04  0.97660  1.03540
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
#> 6.00e+04simulations finished in0.4minutes.(2477per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy      mean std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>      <dbl>   <dbl>       <dbl>        <dbl>
#> 1 policy1 4.19e-11      NA    4.19e-11     4.19e-11
#> 2 policy2 4.19e-11      NA    4.19e-11     4.19e-11
```

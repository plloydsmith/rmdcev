---
output: 
  html_document: 
    keep_md: yes
---



# Multiple Discrete-Continuous Extreme Value (MDCEV) Model Estimation and Simulation in R: The rmdcev Package 

<!-- badges: start -->
[![Build Status](https://travis-ci.org/plloydsmith/rmdcev.svg?branch=master)](https://travis-ci.org/plloydsmith/rmdcev)
<!-- badges: end -->

The _rmdcev_ R package estimates and simulates multiple discrete-continuous extreme value (MDCEV) demand models with observed and unobserved individual heterogneity. Fixed parameter, latent class, and random parameter models can be estimated. These models are estimated using maximum likelihood or Bayesian estimation techniques and are implemented in Stan, which is a C++ package for performing full Bayesian inference (see http://mc-stan.org/). The package also includes functions for simulating demand and welfare outcomes from policy scenarios.

## Current Status

Development is in progress. Currently users can estimate the following models:

1. Fixed parameter models (maximum likelihood or Bayesian estimation)
2. Latent class models (maximum likelihood estimation)
3. Random parameters models (Bayesian estimation)

Under development

2. Phaneuf/von Haefen Kuhn Tucker model

## Installation

I recommend you first install **rstan** and C++ toolchain by following these steps:

https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

Once **rstan** is installed, you can install the released version of rmdcev from GitHub using devtools

``` r
if (!require(devtools)) {
  install.packages("devtools")
  library(devtools)
}
install_github("plloydsmith/rmdcev", build_vignettes = FALSE)
```

Depending on your computer set-up, you may need to adjust your Makevars file of the .R folder (usually in your computer user folder) to ensure the first two lines are

``` r
CXX14 = $(BINPREF)g++ -m$(WIN) -std=c++1y
CXX14FLAGS=-O3 -mtune=native -march=native
```

You can switch build_vignettes to TRUE but it will take a lot longer to install. If installation fails, please let me know by [filing an issue](https://github.com/plloydsmith/rmdcev/issues).

## References

For more details on the model specification and estimation:

Bhat, C.R. (2008) ["The Multiple Discrete-Continuous Extreme Value (MDCEV) Model: Role of Utility Function Parameters, Identification Considerations, and Model Extensions."](https://www.sciencedirect.com/science/article/pii/S0191261507000677) Transportation Research Part B, 42(3): 274-303.

For more details on the welfare simulation:

Lloyd-Smith, P (2018). ["A New Approach to Calculating Welfare Measures in Kuhn-Tucker Demand Models."](https://www.sciencedirect.com/science/article/pii/S1755534517300994) Journal of Choice Modeling, 26: 19-27

## Estimation

As an example, we can simulate some data using the Hybrid specification. In this example, we are simulating data for 1000 invdividuals and 10 non-numeraire alternatives. 


```r
library(pacman)
p_load(tidyverse, rmdcev)
model <- "hybrid"
nobs <- 1000
nalts <- 10
sim.data <- GenerateMDCEVData(model = model, nobs = nobs, nalts = nalts)
#> Checking data...
#> Data is good
```

Estimate model using MLE

```r
mdcev_est <- mdcev(~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8-1,
				   data = sim.data$data,
				   model = model,
				   algorithm = "MLE")
#> Using MLE to estimate MDCEV
#> Chain 1: Initial log joint probability = -16648.7
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       19      -13713.2      0.568511       266.997           1           1       30   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       39      -13605.2     0.0658926       56.7751      0.6705      0.6705       51   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       59      -13593.6     0.0337973        19.519           1           1       73   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       79      -13591.6     0.0423431        19.982      0.1854           1       95   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:       99      -13590.9     0.0138761        1.7051           1           1      120   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      119      -13590.9    0.00115508      0.330859       0.508       0.508      148   
#> Chain 1:     Iter      log prob        ||dx||      ||grad||       alpha      alpha0  # evals  Notes 
#> Chain 1:      120      -13590.9    0.00164954       0.14598           1           1      149   
#> Chain 1: Optimization terminated normally: 
#> Chain 1:   Convergence detected: relative gradient magnitude is below tolerance
```

Summarize results

```r
summary(mdcev_est)
#> Model run using rmdcev for R, version 1.1.2 
#> Estimation method                : MLE
#> Model type                       : hybrid specification
#> Number of classes                : 1
#> Number of individuals            : 1000
#> Number of non-numeraire alts     : 10
#> Estimated parameters             : 20
#> LL                               : -13590.87
#> AIC                              : 27221.74
#> BIC                              : 27319.9
#> Standard errors calculated using : 50 MVN draws
#> Exit of MLE                      : successful convergence
#> Time taken (hh:mm:ss)            : 00:00:1.42
#> 
#> Average consumption of non-numeraire alternatives:
#>     1     2     3     4     5     6     7     8     9    10 
#> 10.87  5.83  5.94 37.14  5.83 49.17  4.65 18.19 20.62 24.88 
#> 
#> Parameter estimates --------------------------------  
#>         Estimate Std.err z.stat
#> psi_b1    -5.743   0.713  -8.05
#> psi_b2     0.871   0.101   8.59
#> psi_b3     2.055   0.180  11.43
#> psi_b4    -1.610   0.120 -13.39
#> psi_b5     3.058   0.220  13.92
#> psi_b6    -2.049   0.150 -13.68
#> psi_b7     1.047   0.113   9.29
#> psi_b8     2.136   0.162  13.18
#> gamma1     1.533   0.321   4.77
#> gamma2     1.389   0.305   4.56
#> gamma3     1.971   0.378   5.22
#> gamma4     1.677   0.247   6.78
#> gamma5     2.014   0.358   5.62
#> gamma6     1.338   0.168   7.95
#> gamma7     1.039   0.242   4.29
#> gamma8     1.888   0.286   6.60
#> gamma9     1.840   0.364   5.06
#> gamma10    1.368   0.172   7.97
#> alpha1     0.466   0.032  14.73
#> scale      1.073   0.062  17.40
#> Note: Alpha parameter is equal for all alternatives.
```


Compare estimates to true values

```r
parms_true <- tbl_df(sim.data$parms_true) %>%
	mutate(true = as.numeric(true))

output <- tbl_df(mdcev_est[["stan_fit"]][["theta_tilde"]]) %>%
				dplyr::select(-tidyselect::starts_with("log_like"),         
							  -tidyselect::starts_with("sum_log_lik"))

names(output)[1:mdcev_est$parms_info$n_vars$n_parms_total] <- mdcev_est$parms_info$parm_names$all_names

output<- output %>%
		tibble::rowid_to_column("sim_id") %>%
		tidyr::gather(parms, value, -sim_id)

coefs <- output %>%
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
#>  1 alpha1   0.466 0.0317  14.7   0.406  0.513  0.5 
#>  2 gamma1   1.53  0.321    4.77  1.13   2.24   1.36
#>  3 gamma10  1.37  0.172    7.97  1.07   1.69   1.37
#>  4 gamma2   1.39  0.305    4.56  0.948  1.99   1.35
#>  5 gamma3   1.97  0.378    5.22  1.40   2.75   1.57
#>  6 gamma4   1.68  0.247    6.78  1.36   2.11   1.90
#>  7 gamma5   2.01  0.358    5.62  1.52   2.95   1.90
#>  8 gamma6   1.34  0.168    7.95  1.06   1.70   1.40
#>  9 gamma7   1.04  0.242    4.29  0.705  1.59   1.13
#> 10 gamma8   1.89  0.286    6.60  1.42   2.37   1.85
#> 11 gamma9   1.84  0.364    5.06  1.32   2.61   1.73
#> 12 psi_b1  -5.74  0.713   -8.05 -6.91  -4.78  -5   
#> 13 psi_b2   0.871 0.101    8.59  0.689  1.06   0.5 
#> 14 psi_b3   2.06  0.180   11.4   1.74   2.32   2   
#> 15 psi_b4  -1.61  0.120  -13.4  -1.81  -1.37  -1.5 
#> 16 psi_b5   3.06  0.220   13.9   2.68   3.39   3   
#> 17 psi_b6  -2.05  0.150  -13.7  -2.32  -1.74  -2   
#> 18 psi_b7   1.05  0.113    9.29  0.808  1.23   1   
#> 19 psi_b8   2.14  0.162   13.2   1.89   2.40   2   
#> 20 scale    1.07  0.0616  17.4   0.966  1.17   1
```

Compare outputs using a figure

```r
coefs %>%
	ggplot(aes(y = mean, x = true))  +
	geom_point(size=2) +
	geom_text(label=coefs$parms) +
	geom_abline(slope = 1) +
	geom_errorbar(aes(ymin=cl_lo,ymax=cl_hi,width=0.2))
```

<img src="man/figures/README-unnamed-chunk-6-1.png" width="100%" />


## Welfare simulation

Create policy simulations (these are 'no change' policies with no effects)

```r
npols <- 2 # Choose number of policies

policies<-	CreateBlankPolicies(npols, nalts, mdcev_est$stan_data[["dat_psi"]], price_change_only = TRUE)

df_sim <- PrepareSimulationData(mdcev_est, policies)
```

Simulate welfare changes

```r
wtp <- mdcev.sim(df_sim$df_indiv, 
				 df_common = df_sim$df_common, 
				 sim_options = df_sim$sim_options,
				 cond_err = 1, 
				 nerrs = 15, 
				 sim_type = "welfare")
#> Using hybrid approach to simulation
#> Simulating welfare...
#> 
#> 9.00e+05simulations finished in0.15minutes.(98793per second)
summary(wtp)
#> # A tibble: 2 x 5
#>   policy      mean  std.dev `ci_lo2.5%` `ci_hi97.5%`
#>   <chr>      <dbl>    <dbl>       <dbl>        <dbl>
#> 1 policy1 1.23e-11 8.50e-11   -1.43e-10     1.44e-10
#> 2 policy2 1.23e-11 8.50e-11   -1.43e-10     1.44e-10
```

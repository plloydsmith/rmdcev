[![Build Status](https://travis-ci.org/plloydsmith/rmdcev.svg?branch=master)](https://travis-ci.org/plloydsmith/rmdcev)

# R code for multiple discrete-continuous extreme value (MDCEV) model 

The _rmdcev_ R package estimates and simulates multiple discrete-continuous extreme value (MDCEV) demand models with observed and unobserved individual heterogneity. Fixed parameter, latent class, and random parameter models can be estimated. These models are estimated using maximum likelihood or Bayesian estimation techniques and are implemented in Stan, which is a C++ package for performing full Bayesian inference (see http://mc-stan.org/). The package also supports Phaneuf and von Haefen's Kuhn-Tucker model specification.

## Current Status

Development is in progress. Currently users can estimate the following models:

1. Fixed parameter models (maximum likelihood or Bayesian)
2. Latent class models (maximum likelihood)

Under development
1. Random parameters models
2. Phaneuf/von Haefen Kuhn Tucker model

For more details on the model specification and estimation:

Bhat, C.R. (2008) ["The Multiple Discrete-Continuous Extreme Value (MDCEV) Model: Role of Utility Function Parameters, Identification Considerations, and Model Extensions."](https://www.sciencedirect.com/science/article/pii/S0191261507000677) Transportation Research Part B, 42(3): 274-303. 

For more details on the welfare simulation:

Lloyd-Smith, P (2018). ["A New Approach to Calculating Welfare Measures in Kuhn-Tucker Demand Models."](https://www.sciencedirect.com/science/article/pii/S1755534517300994) Journal of Choice Modeling, 26: 19-27


## Installation

I recommend you first install rstan following these steps:

https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

Once rstan is installed, you can install the released version of rmdcev from GitHub using devtools

``` r
library(devtools)
install_github("plloydsmith/rmdcev")
```

## Estimation

As an example, we can simulate some data using the Gamma specification.

```{r} 

library(pacman)
p_load(tidyverse, rmdcev)
model <- "gamma"
sim.data <- GenerateMDCEVData(model = model, nobs = 1000, ngoods = 10)
```

Estimate model using MLE
``` {r}
stan_est <- FitMDCEV(psi_formula = ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8-1,
				   data = sim.data$data,
				   model = model,
				   algorithm = "MLE")
```

Compare estimates to true values
``` {r}
parms_true <- tbl_df(sim.data$parms_true) %>%
	mutate(true = as.numeric(true))

coefs <- stan_est$est_pars %>%
	mutate(parms = gsub("\\[|\\]", "", parms)) %>%
	group_by(parms) %>%
	summarise(mean = mean(value),
			  sd = sd(value),
			  zstat = mean / sd,
			  cl_lo = quantile(value, 0.025),
			  cl_hi = quantile(value, 0.975)) %>%
	left_join(parms_true, by = "parms") %>%
	print(n=200)
```

Compare outputs using a figure
```{r} 
coefs %>%
	ggplot(aes(y = mean, x = true))  +
	geom_point(size=2) +
	geom_text(label=coefs$parms) +
	geom_abline(slope = 1) +
	geom_errorbar(aes(ymin=cl_lo,ymax=cl_hi,width=0.2))
```


## Welfare simulation

Create policy simulations (these are blank policies with no effects)
```{r}
###########################################################################
# Create policy scenarios (can affect price only at this point)
###########################################################################
npols <- 2 # Choose number of policies
policies<-	CreateBlankPolicies(npols, stan_est$stan_data[["J"]], stan_est$stan_data[["dat_psi"]])

# This modifies the policies to lose each alternative at once
# price_p <- cbind(0,diag(stan_est$stan_data[["J"]]))
# policies$price_p <- CreateListsRow(price_p)
```

Simulate welfare changes
```{r}
df_sim <- PrepareSimulationData(stan_est, policies)

wtp <- SimulateMDCEV(df_sim$df_indiv, df_sim$df_common, df_sim$sim_options,
		sim_type = "welfare")

SummaryWelfare(wtp)
```


# R code for multiple discrete-continuous extreme value (MDCEV) model 

The _rmdcev_ R package estimates and simulates multiple discrete-continuous extreme value (MDCEV) demand models with observed and unobserved individual heterogneity. Fixed parameter, latent class, and random parameter models can be estimated. These models are estimated using maximum likelihood or Bayesian estimation techniques and are implemented in Stan, which is a C++ package for performing full Bayesian inference (see http://mc-stan.org/). The package also supports Phaneuf and von Haefen's Kuhn-Tucker model specification.

## Installation

You can install the released version of rmdcev from GitHub using devtools

``` r
library(devtools)
install_github("plloydsmith/rmdcev")
```

Development is in progress. Currently users can estimate the following models:

1. Fixed parameter models
2. Latent class models

Under development
1. Random parameters models
2. Phaneuf/von Haefen Kuhn Tucker model

For more details on the model specification and estimation:

Bhat, C.R. (2008) ["The Multiple Discrete-Continuous Extreme Value (MDCEV) Model: Role of Utility Function Parameters, Identification Considerations, and Model Extensions."](https://www.sciencedirect.com/science/article/pii/S0191261507000677) Transportation Research Part B, 42(3): 274-303. 

For more details on the welfare simulation:

Lloyd-Smith, P (2018). ["A New Approach to Calculating Welfare Measures in Kuhn-Tucker Demand Models."](https://www.sciencedirect.com/science/article/pii/S1755534517300994) Journal of Choice Modeling, 26: 19-27


## Estimation

As an example, let's simulate some data using Gamma specification.

```r 
# MDCEV Model with fixed parameters estimated using MLE or Stan

# Load Packages ------------------------------------#

rm(list=ls(all=TRUE))
ls()

library(pacman)

p_load(tidyverse, rstan, rmdcev)


#---------------------------------------------------------------------------------------
# user input
nobs <- 1000 # number of individuals
inc <- 100000 + runif(nobs, 0, 50000) # budget
ngoods <- 10 # number of goods
price <- 5 + matrix(runif(nobs*ngoods, 0, 100), nobs, ngoods)  # price of non-numeraire good

model <- "gamma0"
beta <- c(-5, 0.5, 2, -1.5, 3, -2, 1, 2)

gamma <- 1 + runif(ngoods, 0, 2)
scale <- 1
fixed_scale <- 0

if (model == "les"){
	model_num <- 1
	alpha <- c(0.8, rep(0, ngoods))
	parms_true <- c(beta, gamma, alpha[1], scale)
} else if (model == "alpha"){
	model_num <- 2
	alpha <- rep(0.7,ngoods+1)#0 + runif(ngoods+1, 0.01, .98)
	gamma <- rep(1, ngoods)
	parms_true <- c(beta, alpha, scale)
} else if (model == "gamma"){
	model_num <- 3
	alpha <- rep(0.5, ngoods+1)
	parms_true <- c(beta, gamma, alpha[1], scale)
} else if (model == "gamma0"){
	model_num <- 4
	alpha <- rep(0, ngoods+1)
	parms_true <- c(beta, gamma, scale)
} else
	stop("No model specificied. Choose a model")

parms_true <- tbl_df(parms_true) %>%
	rename(parms_true = value)

# Create psi variables that vary over alternatives
b1 <- rep(1,ngoods)
b2 <- rbinom(ngoods, 1, 0.5)
b3 <- runif(ngoods, 0, 1)

psi_j_temp <- list(b1 = b1,
				   b2 = b2,
				   b3 = b3)

# Create psi variables that vary by person
psi_socio = 2*matrix(runif(nobs * (length(beta)-3)), nobs,length(beta)-3)

psi_i_temp <- list(b4 = psi_socio[,1],
				   b5 = psi_socio[,2],
				   b5 = psi_socio[,3],
				   b5 = psi_socio[,4],
				   b5 = psi_socio[,5])
```

Collect data into form for simulating quant
```r 
# Create full set of base psi variables
psi_j_temp <- map(psi_j_temp, function(x) {rep(x, each=nobs)})
psi_i_temp <- map(psi_i_temp, function(x) {rep(x, times= ngoods)})

dat_psi = c(psi_j_temp, psi_i_temp)
dat_psi = matrix(unlist(dat_psi), ncol = length(beta))

psi_sims <- matrix(dat_psi %*% beta, ncol = ngoods, byrow = TRUE)
psi_sims <- CreateListsRow(psi_sims)
psi_sims <- list(psi_sims )
names(psi_sims) <- "psi_sims"


inc_list <- list(as.list(inc))
names(inc_list) <- "inc" # price normalized MU at zero

price_list <- cbind(1, price) #add numeraire price to price matrix (<-1)
price_list <- list(CreateListsRow(price_list))
names(price_list) <- "price" # price normalized MU at zero

###########################################################################
# Pull individual level data into one list
###########################################################################

df_temp <- c(inc_list, price_list, psi_sims)
```

Simulate quantites for non-numeraire alternative
```r
expose_stan_functions(stanmodels$SimulationFunctions)

quant <- pmap(df_temp, CalcmdemandOne_rng,
			  gamma_sim=gamma,
			  alpha_sim=alpha,
			  scale_sim=scale,
			  nerrs=nerrs,algo_gen = algo_gen)
```

Convert simulated data into estimation data
```r
quant <- matrix(unlist(quant), nrow = nobs, byrow = TRUE)

quant <- quant[,2:(ncol(quant))]

stan.dat <- list(quant = quant,
			price = price,
			inc = as.vector(inc),
			dat_psi = dat_psi)
```

Estimate model using MLE
``` r
stan_est <- FitMDCEV(stan.dat,
				   model = model,
				   n_classes = n_classes,
				   algorithm = "MLE",
				   hessian = TRUE,
				   n_draws = 30)
```

Compare estimates to true values
``` r
coefs <- stan_est$est_pars %>%
	group_by(parms) %>%
	summarise(mean = mean(value),
			  sd = sd(value),
			  zstat = mean / sd,
			  cl_lo = quantile(value, 0.025),
			  cl_hi = quantile(value, 0.975)) %>%
	bind_cols(parms_true) %>%
	print(n=200)
coefs
```

Compare outputs using a figure
```r 
coefs %>%
	ggplot(aes(y = mean, x = parms_true))  +
	geom_point(size=2) +
	geom_text(label=coefs$parms) +
	geom_abline(slope = 1) +
	geom_errorbar(aes(ymin=cl_lo,ymax=cl_hi,width=0.2))
```


## Welfare simulation

Create policy simulations (these are blank)
```r 
###########################################################################
# Create policy scenarios (can affect price only at this point)
###########################################################################
npols <- 2 # Choose number of policies
policies<-	CreateBlankPolicies(npols, stan_est$stan_data[["J"]], stan_est$stan_data[["dat_psi"]])

# This modifiesthe policies to lose each alternative at once
# price_p <- cbind(0,diag(stan_est$stan_data[["J"]]))
# policies$price_p <- CreateListsRow(price_p)
```

Simulate welfare changes
```r
wtp <- SimulateWTP(stan_est, policies, nsims = 30, nerrs = 30)
```

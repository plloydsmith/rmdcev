context("Test wtp")
#Sys.setenv("R_TESTS" = "")
#data(eggs, package = "rmdcev")
library(pacman)

p_load(rstan, tidyverse, rmdcev)

expose_stan_functions(rmdcev:::stanmodels$SimulationFunctions)

#library(tidyverse, rstan, rmdcev)
nerrs <- 3
n_draws <- 5
nsims <- n_draws
nobs <- 1000 # number of individuals
ngoods <- 4 # number of goods

# default is one
n_classes = 2

n_chains = 1
n_cores = 1
n_iterations = 10
hb_random_parameters = "fixed"

CheckWTP <- function(model, algorithm,
					 cond_error = 1,
					 fixed_scale = 0,
					 trunc_data = 0){

	inc <- 100000 + runif(nobs, 0, 50000) # budget
	price <- 5 + matrix(runif(nobs*ngoods, 0, 100), nobs, ngoods)  # price of non-numeraire good
	beta <- c(-5, 0.5, 2, -1.5, 3, -2, 1, 2)
	gamma <- 1 + runif(ngoods, 0, 2)
	scale <- 1

	if (model == "les"){
		model_num <- 1
		alpha <- c(0.8, rep(0, ngoods))
	} else if (model == "alpha"){
		model_num <- 2
		alpha <- 0 + runif(ngoods+1, 0.01, .98)
		gamma <- rep(1, ngoods)
	} else if (model == "gamma"){
		model_num <- 3
		alpha <- rep(0.5, ngoods+1)
	} else if (model == "gamma0"){
		model_num <- 4
		alpha <- rep(1e-6, ngoods+1)
	} else
		stop("No model specificied. Choose a model")

	# Create psi variables that vary over alternatives
	b1 <- rep(1,ngoods)
	b2 <- c(rep(1,ngoods/2), rep(0, ngoods / 2))
	b3 <- rep(c(0,1), ngoods / 2)

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

	df_temp <- c(inc_list, price_list, psi_sims)

	quant <- pmap(df_temp, CalcmdemandOne_rng,
				  gamma_sim=gamma,
				  alpha_sim=alpha,
				  scale_sim=scale,
				  nerrs=nerrs,algo_gen = 1)

	# Convert simulated data into estimation data
	quant <- matrix(unlist(quant), nrow = nobs, byrow = TRUE)

	quant <- quant[,2:(ncol(quant))]

	stan.dat <- list(quant = quant,
					 price = price,
					 inc = as.vector(inc),
					 dat_psi = dat_psi)

	weights <-  rep(1, length(inc))

	stan_est <- FitMDCEV(stan.dat,
						 model = model,
						 n_classes = n_classes,
						 fixed_scale = fixed_scale,
						 trunc_data = trunc_data,
						 seed = 123,
						 algorithm = algorithm,
						 n_iterations = n_iterations, n_chains = n_chains,
						 hb_random_parameters = hb_random_parameters,
						 print_ll = 0,
						 hessian = TRUE,
						 n_draws = n_draws,
						 keep_loglik = 0)

	npols <- 2
	policies<-	CreateBlankPolicies(npols, stan_est$stan_data[["J"]], stan_est$stan_data[["dat_psi"]])

	wtp <- SimulateWTP(stan_est, policies, nsims = nsims, nerrs = nerrs, cond_error = cond_error)
	wtp_sum <- SummaryWelfare(wtp)
	out <- round(sum(wtp_sum$mean),0)
	return(out)
}

test_that("basic test", {
    # If the number below needs to be increased due to additional outputs,
    # ensure that the output size does not get too big when there are multiple
    # classes and many iterations.
    expect_true(1 < 1245000)
    # Add the following option to expect_true to print out the size in Travis
    # info = print(as.numeric(object.size(result))))
})

test_that("test gamma model", {
out <- CheckWTP(model = "gamma0", algorithm = "MLE")
expect_true(out == 0)
})

#test_that("test gamma0 model", {
#	out2 <- CheckWTP(model = "gamma", algorithm = "MLE")
#	expect_true(out2 == 0)
#})


#out2 <- CheckWTP("gamma0", algorithm = "MLE")
#expect_true(out2 == 0)
#out3 <- CheckWTP("alpha", algorithm = "MLE")
#expect_true(out3 == 0)
#out4 <- CheckWTP("les", algorithm = "MLE")
#expect_true(out4 == 0)
#expect_true(CheckWTP("alpha", algorithm = "MLE") == 0)
#expect_true(CheckWTP("les", algorithm = "MLE") == 0)


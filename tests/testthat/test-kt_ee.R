context("Test KT-EE version")

# Create log(alpha version)
# ensure Rp works
# Need to update SimulationFunctions file....not the case right now
# Email Ben about

tol <- 0.01
#library(tidyverse)
#library(rmdcev)
data(data_rec, package = "rmdcev")

data_rec$beach = ifelse(data_rec$alt == "beach", 1, 0)

data_rec <- mdcev.data(data_rec, subset = id < 100,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

# initial starting values
init = list(psi = array(rep(0,2), dim=c(1,2)),
			phi = array(rep(1,1), dim=c(1,1)),
			scale = array(1, dim=c(1)),
			alpha = array(1-exp(-1), dim=c(1,1)),
			gamma = array(rep(1,1), dim=c(1,1)))

# initial starting values
init = list(psi = rep(0,2),
			phi = rep(1,1),
			scale = 1,
			alpha = (1-exp(-1)),
			gamma = rep(1,1))


# initial starting values
init = list(psi = matrix(rep(0,1), nrow = 1, ncol = 1),
			phi = array(rep(1,1), dim=c(1,1)),
			scale = array(1, dim=c(1)),
			alpha = array(1-exp(-1), dim=c(1,1)),
			gamma = array(rep(1,1), dim=c(1,1)))

test_that("kt_ee model estimation", {

	output <- mdcev(formula = ~ ageindex| 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = 0,
					algorithm = "MLE",
					initial.parameters = init,
					print_iterations = F)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 5)
	expect_equal(output$model, "kt_ee")
	print(output$log.likelihood , digits =10)

	expect_true(abs(output$log.likelihood - (-2770.00194)) < tol)
	expect_true(abs(output$bic - 5562.979479) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["scale"]] - 0.8917970787) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["psi"]][[1]] - -0.2772774938) < tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 1)
})

test_that("kt_ee model estimation using num grad", {

	output <- mdcev(formula = ~ ageindex| 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = 0,
					algorithm = "MLE",
					initial.parameters = init,
					jacobian_analytical_grad = 0,
					print_iterations = F)
	expect_true(abs(output$log.likelihood - (-2770.00194)) < tol)
})

test_that("kt_ee model estimation using trunc_data", {

	output <- mdcev(formula = ~ ageindex| 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = 0,
					algorithm = "MLE",
					trunc_data = 1,
					print_iterations = F)
	expect_true(abs(output$log.likelihood - (-2768.959809)) < tol)
})

test_that("Conditional error draw", {

output <- mdcev(formula = ~ ageindex| 0 | beach,
				data = data_rec,
				model = "kt_ee",
				gamma_ascs = 0,
				algorithm = "MLE",
				print_iterations = F)
nalts <- output$stan_data[["J"]]
model_num <- 5
npols <- 2
algo_gen = 1

policies<-	CreateBlankPolicies(npols = 2, output, price_change_only = TRUE)

df_sim <- PrepareSimulationData(output, policies, nsims = 1)

income <- df_sim[["df_indiv"]][["income"]][[2]]
quant_j <- df_sim[["df_indiv"]][["quant_j"]][[2]]
price <- df_sim[["df_indiv"]][["price"]][[2]]

quant_num <- income - quant_j  %*% price[-1]
quant <- c(quant_num, quant_j)
psi_j <- c(output[["stan_data"]][["dat_psi"]][1:nalts,] %*% t(output[["stan_fit"]][["par"]][["psi"]]))
phi_j <- exp(c(output[["stan_data"]][["dat_phi"]][1:nalts,] %*% t(output[["stan_fit"]][["par"]][["phi"]])))
gamma_j <- rep(output[["stan_fit"]][["par"]][["gamma"]], nalts)
gamma <- c(1, gamma_j)
alpha <- c(output[["stan_fit"]][["par"]][["alpha"]], rep(0, nalts))
#alpha <- rep(output[["stan_fit"]][["par"]][["alpha"]], nalts+1)
#alpha <- c(output[["stan_fit"]][["par"]][["alpha"]], rep(0,nalts))
scale <- output[["stan_fit"]][["par"]][["scale"]]


#		library(rstan)
#		expose_stan_functions("inst/stan/SimulationFunctions.stan")

	tol_e <- 1e-20
	tol_l <- 1e-20
	max_loop = 500

	PRNG <-rstan::get_rng(seed = 3)
	o <- rstan::get_stream() # Need for Expecting an external pointer error

	error <- DrawError_rng(quant_num, quant_j, price[-1],
						   psi_j, phi_j, gamma_j, alpha, scale, model_num = 5,
						   nalts = nalts, nerrs = 1, cond_error = 1, draw_mlhs = 1,
						   PRNG, o)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err*c(1, phi_j) / (price * gamma)
	#	 Test general algo
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen, model_num, tol_e = tol_e, max_loop = max_loop, o)
	expect_true(sum(abs(mdemand - quant)) < tol)

	mdemand <- c(quant_num, quant_j)

	util <- ComputeUtilJ(income, mdemand[-1], price[-1],
						 psi_b_err[-1], phi_j, gamma[-1], alpha,
						 nalts, model_num, o)
	util
#	expect_true(abs(util - 253.1646) < tol)

	price_p <- price + c(0,rep(0,nalts))
	MUzero_p <- psi_b_err*c(1, phi_j) / (price_p * gamma)

	hdemand <- HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
								nalts, algo_gen, model_num, tol_l = tol_l, max_loop = max_loop, o)
	mdemand
	hdemand # Will not equal mdemand here because of errors
	wtp_err <- income - t(price_p) %*% hdemand
	wtp_err
	expect_true(abs(wtp_err) < tol)

	price_p <- price + c(0,rep(1000000,nalts))
	MUzero_p <- psi_b_err*c(1, phi_j) / (price_p * gamma)

	hdemand <- HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
							  nalts, algo_gen, model_num, tol_l = tol_l, max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
#	expect_true(abs(wtp_err - (-375.0025)) < tol)
})

test_that("Test full simulation function", {

	data_rec <- mdcev.data(data_rec, subset = id < 100,
						   id.var = "id",
						   alt.var = "alt",
						   choice = "quant")

	output <- mdcev(formula = ~ 0 | 0 | 0,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = 0,
					algorithm = "MLE",
					print_iterations = F)
	summary(output)
	policies<-	CreateBlankPolicies(npols = 2, output, price_change_only = TRUE)

	df_sim <- PrepareSimulationData(output, policies, nsims = 1)
	# Test conditional errors
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
					 cond_err =1, nerrs = 3, sim_type = "welfare")
	sum_wtp <- summary(wtp)
	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < .01)

	# Test unconditional errors (currently returns -Inf for gamam0)
	#	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
	#						 cond_err = 0, nerrs = 3, sim_type = "welfare")
	#	sum_wtp <- summary(wtp)
	#	expect_true(sum(abs(sum_wtp$Mean)) < .01)

})


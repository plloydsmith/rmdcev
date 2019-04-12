context("Test Simulations")
#Sys.setenv("R_TESTS" = "")

library(pacman)
p_load(tidyverse, rstan, rmdcev)

tol <- 0.00001
data(recreation, package = "rmdcev")

result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
				   data = data,
				   model = "gamma0",
				   algorithm = "MLE",
				   n_draws = 30)
npols <- 2 #stan_est$stan_data[["J"]]
policies<-	CreateBlankPolicies(npols, result$stan_data[["J"]], result$stan_data[["dat_psi"]])

df_wtp <- PrepareSimulationData(result, policies, nsims = 3)
inc <- df_wtp[["df_indiv"]][["inc"]][[1]]
quant_j <- df_wtp[["df_indiv"]][["quant_j"]][[1]]
price <- df_wtp[["df_indiv"]][["price"]][[1]]
psi_sim <- df_wtp[["df_indiv"]][["psi_sim"]][[1]]
psi_p_sim <- df_wtp[["df_indiv"]][["psi_p_sim"]][[1]]
gamma_sim <- df_wtp[["df_common"]][["gamma_sim_list"]]
alpha_sim <- df_wtp[["df_common"]][["alpha_sim_list"]]
scale_sim <- df_wtp[["df_common"]][["scale_sim"]]

sim.data <- SimulateMdcevData(model = "gamma0",  nobs = 5)

price_j <- price[-1]
quant_num <- inc - sum(quant_j *price_j)
quant <- c(quant_num, quant_j)
psi_j <- psi_sim[1,]
gamma_j <- gamma_sim[[1]]
gamma <- c(1, gamma_j)
alpha <- alpha_sim[[1]]
scale <- scale_sim[1]
algo_gen <- 0
ngoods <- length(price_j)

test_that("Conditional error draw", {
error <- DrawError_rng(quant_num, quant_j, price_j,
				  psi_j, gamma_j, alpha, scale,
				  ngoods = ngoods, nerrs = 10, cond_error = 1)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err / price
	# Test hybrid algo
	mdemand <- MarshallianDemand(inc, price, MUzero_b, gamma, alpha, ngoods, algo_gen = 0);
	expect_true(sum(abs(mdemand - quant)) < tol)
	# Test hybrid algo
	mdemand <- MarshallianDemand(inc, price, MUzero_b, gamma, alpha, ngoods, algo_gen = 1);
	expect_true(sum(abs(mdemand - quant)) < tol)

#	expect_true(abs(sum(unlist(error)) - sum(-log(-log(runif(3000*18, 0, 1))) * scale)) < tol)

#	tt <- -log(-log(runif(300*18, 0, 1))) * scale
})

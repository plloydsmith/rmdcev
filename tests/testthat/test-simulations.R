context("Test Simulations")
#Sys.setenv("R_TESTS" = "")

library(pacman)
p_load(tidyverse, rstan, rmdcev)

tol <- 0.00001
data(recreation, package = "rmdcev")

result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
				   data = data,
				   model = "gamma0",
				   algorithm = "MLE")
npols <- 2 #stan_est$stan_data[["J"]]
policies<-	CreateBlankPolicies(npols, result$stan_data[["J"]], result$stan_data[["dat_psi"]])

df_wtp <- PrepareSimulationData(result, policies, nsims = 3)
inc <- df_wtp[["df_indiv"]][["inc"]][[1]]
quant_j <- df_wtp[["df_indiv"]][["quant_j"]][[1]]
price <- df_wtp[["df_indiv"]][["price"]][[1]]

#sim.data <- SimulateMdcevData(model = "gamma0",  nobs = 5)
expose_stan_functions(stanmodels$SimulationFunctions)

price_j <- price[-1]
ngoods <- length(price_j)
model_num <- 4

quant_num <- inc - quant_j  %*% price[-1]
quant <- c(quant_num, quant_j)
psi_j <- result[["stan_data"]][["dat_psi"]][1:ngoods,] %*% result[["stan_fit"]][["par"]][["psi"]]
gamma_j <- result[["stan_fit"]][["par"]][["gamma"]]
gamma <- c(1, gamma_j)
alpha <- rep(1e-06, ngoods+1)# result[["stan_fit"]][["par"]][["alpha"]]
scale <- result[["stan_fit"]][["par"]][["scale"]]

test_that("Conditional error draw", {

tol_e <- 1e-20
tol_l <- 1e-20
max_loop = 999

error <- DrawError_rng(quant_num, quant_j, price[-1],
				  psi_j, gamma_j, alpha, scale,
				  ngoods = ngoods, nerrs = 10, cond_error = 1)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err / price
	# Test hybrid algo
	mdemand <- MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
								 ngoods, algo_gen = 0, tol_e = tol_e, max_loop = max_loop);
	expect_true(sum(abs(mdemand - quant)) < tol)
	# Test general algo
	mdemand <- MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
								 ngoods, algo_gen = 1, tol_e = tol_e, max_loop = max_loop);
	expect_true(sum(abs(mdemand - quant)) < tol)

error <- c(0.0000000000, -0.3971093521, -0.5374972429,  1.3553488637, -0.7734435227, -0.3326890061, 0.0008661508, -0.3271314442,
  -0.8966159933, -0.1718724179,  0.7474894252, 1.4910795290, -0.3873935380,  0.0875434277,  0.1386813642, -0.0456836959,  0.1453761206,  0.4672311416)

psi_b_err <- exp(c(0, psi_j) + error)
MUzero_b <- psi_b_err / price
mdemand <- MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
							 ngoods, algo_gen = 0, tol_e = tol_e, max_loop = max_loop);

	util <- ComputeUtilJ(inc, mdemand[-1], price[-1],
							 psi_b_err[-1], gamma[-1], alpha,
							 ngoods, model_num)
	expect_true(abs(util - 1000011.39019661) < tol)

	price_p <- price + c(0,rep(1,ngoods))
 	MUzero_p <- psi_b_err / price_p

	hdemand <- HicksianDemand(util, price_p, MUzero_p, gamma, alpha,
			ngoods, algo_gen = 0, model_num, tol_l = tol_l, max_loop = max_loop)
	wtp_err <- inc - t(price_p) %*% hdemand
	expect_true(abs(wtp_err - (-60.5094477323)) < tol)

	hdemand <- HicksianDemand(util, price_p, MUzero_p, gamma, alpha,
							  ngoods, algo_gen = 1, model_num, tol_l = tol_l, max_loop = max_loop)
	wtp_err <- inc - t(price_p) %*% hdemand
	expect_true(abs(wtp_err - (-60.5094477323)) < tol)

})

test_that("Test full simulation function", {

	# Test conditional errors
	wtp <- SimulateWTP(df_wtp$df_indiv, df_common = df_wtp$df_common, sim_options = df_wtp$sim_options, cond_err = 1, nerrs = 3)
	sum_wtp <- SummaryWelfare(wtp)
	expect_true(sum(abs(sum_wtp$mean)) < tol)

	# Test unconditional errors
	wtp <- SimulateWTP(df_wtp$df_indiv, df_common = df_wtp$df_common, sim_options = df_wtp$sim_options, cond_err = 0, nerrs = 3)
	sum_wtp <- SummaryWelfare(wtp)
	expect_true(sum(abs(sum_wtp$mean)) < tol)

})


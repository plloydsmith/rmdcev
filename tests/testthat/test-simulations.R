context("Test Simulations")

tol <- 0.00001
data(data_rec, package = "rmdcev")
data_rec
result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
				   data = subset(data_rec, id < 100),
				   model = "gamma0",
				   algorithm = "MLE",
				   print_iterations = FALSE)
npols <- 2
policies<-	CreateBlankPolicies(npols, result$stan_data[["J"]], result$stan_data[["dat_psi"]])

df_sim <- PrepareSimulationData(result, policies, nsims = 3)
inc <- df_sim[["df_indiv"]][["inc"]][[1]]
quant_j <- df_sim[["df_indiv"]][["quant_j"]][[1]]
price <- df_sim[["df_indiv"]][["price"]][[1]]

#sim.data <- GenerateMDCEVData(model = "gamma0",  nobs = 5)
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
				  ngoods = ngoods, nerrs = 2, cond_error = 1, draw_mlhs = 1)

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

error <- c(0.0000000000, -0.77612399, -0.55169780, -0.02143232, -0.81241994, -0.19768961,
		 -0.85621517, -0.82511555, -0.60241553, -1.14157948, -0.62773715, 0.69217576,
		0.61801650, 0.26789713, -0.66958419, -0.35102966, -0.90048354, 0.32263623)

psi_b_err <- exp(c(0, psi_j) + error)
MUzero_b <- psi_b_err / price
mdemand <- MarshallianDemand(inc, price, MUzero_b, gamma, alpha,
							 ngoods, algo_gen = 0, tol_e = tol_e, max_loop = max_loop);

	util <- ComputeUtilJ(inc, mdemand[-1], price[-1],
							 psi_b_err[-1], gamma[-1], alpha,
							 ngoods, model_num)
	expect_true(abs(util - 1000011.04297481) < tol)

	price_p <- price + c(.001,rep(1,ngoods))
 	MUzero_p <- psi_b_err / price_p

	hdemand <- HicksianDemand(util, price_p, MUzero_p, gamma, alpha,
			ngoods, algo_gen = 0, model_num, tol_l = tol_l, max_loop = max_loop)
	wtp_err <- inc - t(price_p) %*% hdemand
	expect_true(abs(wtp_err - (-62.4994989953)) < tol)

	hdemand <- HicksianDemand(util, price_p, MUzero_p, gamma, alpha,
							  ngoods, algo_gen = 1, model_num, tol_l = tol_l, max_loop = max_loop)
	wtp_err <- inc - t(price_p) %*% hdemand
	expect_true(abs(wtp_err - (-62.4994989953)) < tol)

})

test_that("Test full simulation function", {

	# Test conditional errors
	wtp <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 3, sim_type = "welfare")
	sum_wtp <- SummaryWelfare(wtp)
	expect_true(sum(abs(sum_wtp$mean)) < .01)

	# Test unconditional errors
	wtp <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 0, nerrs = 3, sim_type = "welfare")
	sum_wtp <- SummaryWelfare(wtp)
	expect_true(sum(abs(sum_wtp$mean)) < .01)

})


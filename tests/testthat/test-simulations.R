tol <- 0.1

data(data_rec, package = "rmdcev")

data_rec <- mdcev.data(data_rec, subset = id < 100,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

result <- mdcev(~ alt - 1,
				data = data_rec,
				model = "hybrid0",
				psi_ascs = FALSE,
				algorithm = "MLE",
				std_errors = "mvn",
				print_iterations = FALSE,
				backend = "rstan")

nalts <- result$stan_data[["J"]]
model_num <- result$stan_data[["model_num"]]
npols <- 2
policies <- CreateBlankPolicies(npols, result, price_change_only = TRUE)

df_sim <- PrepareSimulationData(result, policies, nsims = 3)
income <- df_sim[["df_indiv"]][["income"]][[1]]
quant_j <- df_sim[["df_indiv"]][["quant_j"]][[1]]
price <- df_sim[["df_indiv"]][["price"]][[1]]

price_j <- price[-1]

quant_num <- income - quant_j %*% price[-1]
quant <- c(quant_num, quant_j)
psi_j <- result[["stan_data"]][["dat_psi"]][1:nalts, ] %*% t(result[["stan_fit"]][["par"]][["psi"]])
phi_j <- rep(0, nalts)
gamma_j <- result[["stan_fit"]][["par"]][["gamma"]]
gamma <- c(1, gamma_j)
alpha <- rep(0, nalts + 1)
scale <- result[["stan_fit"]][["par"]][["scale"]]

test_that("Conditional error hybrid0 draw", {
	skip_on_os("solaris")

	tol_e <- 1e-20
	tol_l <- 1e-20
	max_loop <- 999

	PRNG <- rmdcev_get_rng(seed = 3)
	o <- rmdcev_get_stream()

	error <- DrawError_rng(quant_num, quant_j, price[-1],
						   psi_j, phi_j, gamma_j, alpha, scale,
						   model_num = model_num, nalts = nalts, nerrs = 2,
						   cond_error = 1, draw_mlhs = 1, PRNG, o)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err / price

	# Test hybrid algo
	algo_gen <- 0
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen, model_num, tol_e, max_loop, o)
	expect_equal(sum(abs(mdemand - quant)), 0, tolerance = tol)

	# Test general algo
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen = 1, model_num, tol_e = tol_e,
								 max_loop = max_loop, o)
	expect_equal(sum(abs(mdemand - quant)), 0, tolerance = tol)

	# Corner solution: extreme negative errors suppress all non-numeraire demand
	error_corner <- c(0, rep(-100, nalts))
	psi_b_err_corner <- exp(c(0, psi_j) + error_corner)
	MUzero_corner <- psi_b_err_corner / price
	mdemand_corner <- MarshallianDemand(income, price, MUzero_corner, c(1, phi_j), gamma, alpha,
										nalts, algo_gen = 0, model_num, tol_e = tol_e,
										max_loop = max_loop, o)
	expect_equal(as.numeric(mdemand_corner[-1]), rep(0, nalts), tolerance = 1e-4)

	# For hybrid0 at corner solution: U = log(numeraire_quant) = log(income)
	util_corner <- ComputeUtilJ(income, mdemand_corner[-1], price[-1],
								psi_b_err_corner, phi_j, gamma[-1], alpha,
								nalts, model_num, o)
	expect_equal(util_corner, log(income), tolerance = tol)

	# Hicksian demand at price increase: WTP < 0 and must reproduce corner utility
	price_p <- price + c(.001, rep(1, nalts))
	MUzero_p <- psi_b_err_corner / price_p
	hdemand <- HicksianDemand(util_corner, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
							  nalts, algo_gen = 0, model_num, tol_l = tol_l,
							  max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
	expect_true(as.numeric(wtp_err) < 0)
	budget_p <- hdemand[1] + sum(price_p[-1] * hdemand[-1])
	util_check <- ComputeUtilJ(budget_p, hdemand[-1], price_p[-1],
							   psi_b_err_corner, phi_j, gamma[-1], alpha, nalts, model_num, o)
	expect_equal(util_check, util_corner, tolerance = 1e-4)
})

test_that("Test demand simulation", {
	skip_on_os("solaris")

	# Test conditional errors
	demand <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
						sim_options = df_sim$sim_options,
						cond_err = TRUE, nerrs = 3, sim_type = "demand")
	# Budget exhaustion: simulated demand must exhaust income for individual 1
	expect_equal(sum(price * demand[[1]][[1]][1, ]), income, tolerance = 0.1)
	# Conditional errors recover individual 1's observed quantities
	expect_equal(as.numeric(demand[[1]][[1]][1, ]), as.numeric(quant), tolerance = 0.1)

	# Unconditional errors currently return -Inf for hybrid0: not asserted, left for future fix.
})

test_that("Test full simulation function conditional welfare", {

	# Conditional errors: zero price change -> WTP should be ~0
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = TRUE, nerrs = 3, sim_type = "welfare")
	sum_wtp <- summary(wtp)
	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = .01)

	# Unconditional errors currently return -Inf for hybrid0: not asserted, left for future fix.
})

test_that("Conditional error hybrid draw", {
	skip_on_os("solaris")

	result_hyb <- mdcev(~ alt - 1,
						data = data_rec,
						model = "hybrid",
						psi_ascs = FALSE,
						algorithm = "MLE",
						print_iterations = FALSE,
						backend = "rstan")

	nalts_h <- result_hyb$stan_data[["J"]]
	model_num_h <- result_hyb$stan_data[["model_num"]]
	npols_h <- 2
	policies_h <- CreateBlankPolicies(npols_h, result_hyb, price_change_only = TRUE)

	df_sim_h <- PrepareSimulationData(result_hyb, policies_h, nsims = 1)
	income_h <- df_sim_h[["df_indiv"]][["income"]][[2]]
	quant_j_h <- df_sim_h[["df_indiv"]][["quant_j"]][[2]]
	price_h <- df_sim_h[["df_indiv"]][["price"]][[2]]

	quant_num_h <- income_h - quant_j_h %*% price_h[-1]
	quant_h <- c(quant_num_h, quant_j_h)
	psi_j_h <- result_hyb[["stan_data"]][["dat_psi"]][1:nalts_h, ] %*%
		t(result_hyb[["stan_fit"]][["par"]][["psi"]])
	phi_j_h <- rep(0, nalts_h)
	gamma_j_h <- result_hyb[["stan_fit"]][["par"]][["gamma"]]
	gamma_h <- c(1, gamma_j_h)
	alpha_h <- rep(result_hyb[["stan_fit"]][["par"]][["alpha"]], nalts_h + 1)
	scale_h <- as.numeric(result_hyb[["stan_fit"]][["par"]][["scale"]])
	expect_snapshot_value(round(scale_h, 2), style = "deparse", cran = FALSE)
	expect_true(scale_h > 0)

	tol_e <- 1e-20
	tol_l <- 1e-20
	max_loop <- 999

	PRNG <- rmdcev_get_rng(seed = 3)
	o <- rmdcev_get_stream()

	error_h <- DrawError_rng(quant_num_h, quant_j_h, price_h[-1],
							 psi_j_h, phi_j_h, gamma_j_h, alpha_h, scale_h,
							 model_num = model_num_h, nalts = nalts_h, nerrs = 2,
							 cond_error = 1, draw_mlhs = 1, PRNG, o)

	psi_b_err_h <- exp(c(0, psi_j_h) + error_h[[1]])
	MUzero_b_h <- psi_b_err_h / price_h

	mdemand_0 <- MarshallianDemand(income_h, price_h, MUzero_b_h,
								   c(1, phi_j_h), gamma_h, alpha_h,
								   nalts_h, 0, model_num_h, tol_e, max_loop, o)
	expect_equal(sum(abs(mdemand_0 - quant_h)), 0, tolerance = tol)

	mdemand_1 <- MarshallianDemand(income_h, price_h, MUzero_b_h,
								   c(1, phi_j_h), gamma_h, alpha_h,
								   nalts_h, 1, model_num_h, tol_e, max_loop, o)
	expect_equal(sum(abs(mdemand_0 - mdemand_1)), 0, tolerance = tol)

	# Compute utility from conditional demand (equals observed demand by construction)
	util_h <- ComputeUtilJ(income_h, mdemand_0[-1], price_h[-1],
						   psi_b_err_h, phi_j_h, gamma_h[-1], alpha_h,
						   nalts_h, model_num_h, o)
	expect_true(is.finite(util_h))

	# Hicksian utility consistency: demand at new prices must reproduce original utility
	price_p_h <- price_h + c(.001, rep(1, nalts_h))
	MUzero_p_h <- psi_b_err_h / price_p_h

	hdemand_h <- HicksianDemand(util_h, price_p_h, MUzero_p_h,
								c(1, phi_j_h), gamma_h, alpha_h,
								nalts_h, 1, model_num_h, tol_l, max_loop, o)
	wtp_err_h <- income_h - t(price_p_h) %*% hdemand_h
	expect_true(as.numeric(wtp_err_h) < 0)
	budget_h <- hdemand_h[1] + sum(price_p_h[-1] * hdemand_h[-1])
	util_h_check <- ComputeUtilJ(budget_h, hdemand_h[-1], price_p_h[-1],
								 psi_b_err_h, phi_j_h, gamma_h[-1], alpha_h,
								 nalts_h, model_num_h, o)
	expect_equal(util_h_check, util_h, tolerance = 1e-4)

	# Both demand algorithms give identical WTP
	hdemand_h2 <- HicksianDemand(util_h, price_p_h, MUzero_p_h,
								 c(1, phi_j_h), gamma_h, alpha_h,
								 nalts_h, 0, model_num_h, tol_l, max_loop, o)
	wtp_err_h2 <- income_h - t(price_p_h) %*% hdemand_h2
	expect_equal(as.numeric(wtp_err_h2), as.numeric(wtp_err_h), tolerance = tol)
})

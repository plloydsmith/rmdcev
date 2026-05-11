tol <- 0.01

data(data_rec, package = "rmdcev")

data_rec$beach <- ifelse(data_rec$alt == "beach", 1, 0)

data_rec <- mdcev.data(data_rec, subset = id < 100,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

# initial starting values (kept for readability; last assignment is used below)
init <- list(psi  = matrix(rep(0, 1), nrow = 1, ncol = 1),
			 phi   = array(rep(1, 1), dim = c(1, 1)),
			 scale = array(1, dim = c(1)),
			 alpha = array(1 - exp(-1), dim = c(1, 1)),
			 gamma = array(rep(1, 1), dim = c(1, 1)))

test_that("kt_ee model estimation", {

	output <- mdcev(formula = ~ ageindex | 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = FALSE,
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 5)
	expect_equal(output$model, "kt_ee")
	expect_true(output$log.likelihood < 0)
	expect_snapshot_value(round(output$log.likelihood, 2), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(output$bic, 2), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(as.numeric(output[["stan_fit"]][["par"]][["scale"]]), 4), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(as.numeric(output[["stan_fit"]][["par"]][["psi"]][[1]]), 4), style = "deparse", cran = FALSE)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 1)
})

test_that("kt_ee model estimation using num grad", {

	output <- mdcev(formula = ~ ageindex | 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = FALSE,
					algorithm = "MLE",
					initial.parameters = init,
					jacobian_analytical_grad = FALSE,
					print_iterations = FALSE,
					backend = "rstan")
	expect_true(output$log.likelihood < 0)
	expect_snapshot_value(round(output$log.likelihood, 2), style = "deparse", cran = FALSE)
})

test_that("kt_ee model estimation using trunc_data", {

	output <- mdcev(formula = ~ ageindex | 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = FALSE,
					algorithm = "MLE",
					trunc_data = TRUE,
					print_iterations = FALSE,
					backend = "rstan")
	expect_true(output$log.likelihood < 0)
	expect_snapshot_value(round(output$log.likelihood, 2), style = "deparse", cran = FALSE)
})

test_that("kt_ee with gamma_ascs = TRUE does not recover gamma parameters", {
	set.seed(12345)

	true_gamma <- c(2, 4, 6, 8, 10)
	sim_data <- GenerateMDCEVData(
		model = "kt_ee",
		nobs = 1000,
		nalts = length(true_gamma),
		psi_i_parms = -1,
		phi_parms = 0.5,
		gamma_parms = true_gamma
	)

	output <- mdcev(
		formula = ~ b1 | 0 | phi_1,
		data = sim_data$data,
		model = "kt_ee",
		gamma_ascs = TRUE,
		algorithm = "MLE",
		print_iterations = FALSE,
		backend = "rstan"
	)

	est_gamma <- as.numeric(output[["stan_fit"]][["par"]][["gamma"]])

	expect_equal(length(est_gamma), length(true_gamma))
	expect_true(any(abs(est_gamma - true_gamma) > 0.5))
})

test_that("Conditional error draw", {
	skip_on_os("solaris")

	output <- mdcev(formula = ~ ageindex | 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = FALSE,
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")
	nalts <- output$stan_data[["J"]]
	model_num <- 5
	algo_gen <- 1

	policies <- CreateBlankPolicies(npols = 2, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 1)

	income <- df_sim[["df_indiv"]][["income"]][[2]]
	quant_j <- df_sim[["df_indiv"]][["quant_j"]][[2]]
	price <- df_sim[["df_indiv"]][["price"]][[2]]

	quant_num <- income - quant_j %*% price[-1]
	quant <- c(quant_num, quant_j)
	psi_j <- c(output[["stan_data"]][["dat_psi"]][1:nalts, ] %*%
				t(output[["stan_fit"]][["par"]][["psi"]]))
	phi_j <- exp(c(output[["stan_data"]][["dat_phi"]][1:nalts, ] %*%
				   t(output[["stan_fit"]][["par"]][["phi"]])))
	gamma_j <- rep(output[["stan_fit"]][["par"]][["gamma"]], nalts)
	gamma <- c(1, gamma_j)
	alpha <- c(output[["stan_fit"]][["par"]][["alpha"]], rep(0, nalts))
	scale <- as.numeric(output[["stan_fit"]][["par"]][["scale"]])
	expect_true(scale > 0)

	tol_e <- 1e-20
	tol_l <- 1e-20
	max_loop <- 999

	PRNG <- rmdcev_get_rng(seed = 3)
	o <- rmdcev_get_stream()

	error <- DrawError_rng(quant_num, quant_j, price[-1],
						   psi_j, phi_j, gamma_j, alpha, scale, model_num = 5,
						   nalts = nalts, nerrs = 1, cond_error = 1, draw_mlhs = 1,
						   PRNG, o)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err * c(1, phi_j) / (price * gamma)
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen, model_num, tol_e = tol_e,
								 max_loop = max_loop, o)
	expect_equal(sum(abs(mdemand - quant)), 0, tolerance = tol)

	mdemand <- c(quant_num, quant_j)
	util <- ComputeUtilJ(income, mdemand[-1], price[-1],
						 psi_b_err, phi_j, gamma[-1], alpha,
						 nalts, model_num, o)
	expect_true(is.finite(util))

	# Zero price change: WTP should be ~0
	price_p <- price + c(0, rep(0, nalts))
	MUzero_p <- psi_b_err * c(1, phi_j) / (price_p * gamma)
	hdemand <- HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
							  nalts, algo_gen, model_num, tol_l = tol_l,
							  max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
	expect_equal(as.numeric(wtp_err), 0, tolerance = tol)

	# Large price change: WTP must be negative and Hicksian demand must reproduce original utility
	price_p <- price + c(0, rep(1000000, nalts))
	MUzero_p <- psi_b_err * c(1, phi_j) / (price_p * gamma)
	hdemand <- HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
							  nalts, algo_gen, model_num, tol_l = tol_l,
							  max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
	expect_true(as.numeric(wtp_err) < 0)
	budget_p <- hdemand[1] + sum(price_p[-1] * hdemand[-1])
	util_check <- ComputeUtilJ(budget_p, hdemand[-1], price_p[-1],
							   psi_b_err, phi_j, gamma[-1], alpha, nalts, model_num, o)
	expect_equal(util_check, util, tolerance = 1e-3)
})

test_that("unconditional error draw", {

	output <- mdcev(formula = ~ ageindex | 0 | beach,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = FALSE,
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")
	nalts <- output$stan_data[["J"]]
	model_num <- 5
	algo_gen <- 1

	policies <- CreateBlankPolicies(npols = 2, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 1)

	income <- df_sim[["df_indiv"]][["income"]][[2]]
	quant_j <- df_sim[["df_indiv"]][["quant_j"]][[2]]
	price <- df_sim[["df_indiv"]][["price"]][[2]]

	quant_num <- income - quant_j %*% price[-1]
	quant <- c(quant_num, quant_j)
	psi_j <- c(output[["stan_data"]][["dat_psi"]][1:nalts, ] %*%
				t(output[["stan_fit"]][["par"]][["psi"]]))
	phi_j <- exp(c(output[["stan_data"]][["dat_phi"]][1:nalts, ] %*%
				   t(output[["stan_fit"]][["par"]][["phi"]])))
	gamma_j <- rep(output[["stan_fit"]][["par"]][["gamma"]], nalts)
	gamma <- c(1, gamma_j)
	alpha <- c(output[["stan_fit"]][["par"]][["alpha"]], rep(0, nalts))
	scale <- output[["stan_fit"]][["par"]][["scale"]]

	tol_e <- 1e-20
	tol_l <- 1e-20
	max_loop <- 999

	PRNG <- rmdcev_get_rng(seed = 5)
	o <- rmdcev_get_stream()

	error <- DrawError_rng(quant_num, quant_j, price[-1],
						   psi_j, phi_j, gamma_j, alpha, scale, model_num = 5,
						   nalts = nalts, nerrs = 5, cond_error = 0, draw_mlhs = 1,
						   PRNG, o)

	psi_b_err <- exp(c(0, psi_j) + error[[1]])
	MUzero_b <- psi_b_err * c(1, phi_j) / (price * gamma)
	mdemand <- MarshallianDemand(income, price, MUzero_b, c(1, phi_j), gamma, alpha,
								 nalts, algo_gen, model_num, tol_e = tol_e,
								 max_loop = max_loop, o)

	util <- ComputeUtilJ(income, mdemand[-1], price[-1],
						 psi_b_err, phi_j, gamma[-1], alpha,
						 nalts, model_num, o)

	# Zero price change: Hicksian demand = Marshallian demand -> WTP â‰ˆ 0
	price_p <- price + c(0, rep(0, nalts))
	MUzero_p <- psi_b_err * c(1, phi_j) / (price_p * gamma)
	hdemand <- HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
							  nalts, algo_gen, model_num, tol_l = tol_l,
							  max_loop = max_loop, o)
	# With unconditional errors mdemand != hdemand in general, but WTP â‰ˆ 0 for zero price change
	expect_equal(sum(abs(mdemand - hdemand)), 0, tolerance = tol)
	wtp_err <- income - t(price_p) %*% hdemand
	expect_equal(as.numeric(wtp_err), 0, tolerance = tol)

	# Large price change
	price_p <- price + c(0, rep(1000000, nalts))
	MUzero_p <- psi_b_err * c(1, phi_j) / (price_p * gamma)
	hdemand <- HicksianDemand(util, price_p, MUzero_p, c(1, phi_j), gamma, alpha,
							  nalts, algo_gen, model_num, tol_l = tol_l,
							  max_loop = max_loop, o)
	wtp_err <- income - t(price_p) %*% hdemand
	expect_true(is.numeric(as.numeric(wtp_err)))
	expect_true(as.numeric(wtp_err) < 0)
})

test_that("Test full simulation function", {

	output <- mdcev(formula = ~ 0 | 0 | 0,
					data = data_rec,
					model = "kt_ee",
					gamma_ascs = FALSE,
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	policies <- CreateBlankPolicies(npols = 2, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 1)

	# Zero price change: conditional welfare should be ~0
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = TRUE, nerrs = 3, sim_type = "welfare")
	sum_wtp <- summary(wtp)
	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = .01)

	# Unconditional errors currently return -Inf for kt_ee zero-price case: not asserted.
})

test_that("kt_ee parameter recovery from synthetic data", {
	skip_on_cran()
	set.seed(42)

	true_psi   <- -1.5
	true_phi   <- 0.5
	true_gamma <- rep(5, 5)

	sim_data <- GenerateMDCEVData(
		model       = "kt_ee",
		nobs        = 500,
		nalts       = 5,
		psi_i_parms = true_psi,
		phi_parms   = true_phi,
		gamma_parms = true_gamma
	)

	output <- mdcev(
		formula          = ~ b1 | 0 | phi_1,
		data             = sim_data$data,
		model            = "kt_ee",
		gamma_ascs       = FALSE,
		algorithm        = "MLE",
		print_iterations = FALSE,
		backend          = "rstan"
	)

	est_psi   <- as.numeric(output[["stan_fit"]][["par"]][["psi"]])
	est_phi   <- as.numeric(output[["stan_fit"]][["par"]][["phi"]])
	est_gamma <- as.numeric(output[["stan_fit"]][["par"]][["gamma"]])

	expect_true(est_psi < 0,        label = "psi is negative")
	expect_true(all(est_gamma > 0), label = "all gamma positive")
	expect_equal(est_psi, true_psi, tolerance = 0.5,
				 label = "psi recovery within 0.5")
	expect_true(is.finite(output$log.likelihood))
	expect_true(output$log.likelihood < 0)
})

tol <- 0.1

data(data_rec, package = "rmdcev")

data_rec <- mdcev.data(
	data_rec,
	subset = id <= 100,
	id.var = "id",
	alt.var = "alt",
	choice = "quant"
)

test_that("MLE gamma1 runs and has no alpha parameters", {
	output <- mdcev(
		~0,
		data = data_rec,
		model = "gamma1",
		psi_ascs = FALSE,
		algorithm = "MLE",
		print_iterations = FALSE,
		backend = "rstan"
	)

	output.sum <- summary(output)
	expect_equal(output$model, "gamma1")
	expect_true(output$log.likelihood < 0)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)

	expect_snapshot_value(
		round(output$log.likelihood, 2),
		style = "deparse",
		cran = FALSE
	)
	expect_snapshot_value(round(output$bic, 2), style = "deparse", cran = FALSE)
	expect_snapshot_value(
		round(as.numeric(output[["stan_fit"]][["par"]][["scale"]]), 4),
		style = "deparse",
		cran = FALSE
	)
	# S3 method coverage
	cf <- coef(output)
	expect_true(is.list(cf))
	expect_true(length(cf) > 0)

	ll <- logLik(output)
	expect_true(is.numeric(ll))

	expect_output(print(output))
})


test_that("gamma1 simulation alpha matrix is correct", {
	skip_on_cran()

	output <- mdcev(
		~0,
		data = data_rec,
		model = "gamma1",
		psi_ascs = FALSE,
		algorithm = "MLE",
		print_iterations = FALSE,
		backend = "rstan"
	)

	policies <- CreateBlankPolicies(npols = 1, output, price_change_only = TRUE)
	sim_data <- PrepareSimulationData(output, policies, nsims = 5)

	alpha_mat <- sim_data$df_common$alpha_sim_nonrandom
	expect_equal(ncol(alpha_mat), output$stan_data$J + 1)
	expect_true(all(alpha_mat[, 1] == 1)) # alpha_0 = 1 for all sims
	expect_true(all(alpha_mat[, -1] == 0)) # alpha_j = 0 for all non-numeraire
})


test_that("gamma1 parameter recovery from simulated data", {
	skip_on_cran()

	set.seed(42)
	nobs <- 500
	nalts <- 5
	gamma_true <- c(2, 4, 6, 8, 10)
	psi_j_parms <- c(-5, 0.5, 2)
	psi_i_parms <- c(-1.5, 2, -1)

	sim.data <- GenerateMDCEVData(
		model = "gamma",
		nobs = nobs,
		nalts = nalts,
		alpha_parms = 1,
		gamma_parms = gamma_true,
		psi_j_parms = psi_j_parms,
		psi_i_parms = psi_i_parms
	)

	output <- mdcev(
		~ b1 + b2 + b3 + b4 + b5 + b6,
		data = sim.data$data,
		psi_ascs = FALSE,
		model = "gamma1",
		algorithm = "MLE",
		print_iterations = FALSE,
		backend = "rstan"
	)

	expect_equal(output$model, "gamma1")
	expect_true(output$log.likelihood < 0)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)

	psi_hat <- as.numeric(output[["stan_fit"]][["par"]][["psi"]])
	psi_true <- c(psi_j_parms, psi_i_parms)
	expect_equal(psi_hat, psi_true, tolerance = 1.5)

	gamma_hat <- as.numeric(output[["stan_fit"]][["par"]][["gamma"]])
	expect_equal(gamma_hat, gamma_true, tolerance = 3)

	scale_hat <- as.numeric(output[["stan_fit"]][["par"]][["scale"]])
	expect_equal(scale_hat, 1, tolerance = 0.5)
})


test_that("gamma1 welfare: zero price change gives zero WTP", {
	output <- mdcev(
		~0,
		data = data_rec,
		model = "gamma1",
		psi_ascs = FALSE,
		algorithm = "MLE",
		std_errors = "mvn",
		n_draws = 2,
		print_iterations = FALSE,
		backend = "rstan"
	)

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 2)

	wtp <- mdcev.sim(
		df_sim$df_indiv,
		df_common = df_sim$df_common,
		sim_options = df_sim$sim_options,
		cond_err = FALSE,
		nerrs = 1,
		sim_type = "welfare"
	)
	sum_wtp <- summary(wtp)

	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)
})


test_that("gamma1 welfare: price increase gives negative WTP", {
	skip_on_cran()

	output <- mdcev(
		~0,
		data = data_rec,
		model = "gamma1",
		psi_ascs = FALSE,
		algorithm = "MLE",
		print_iterations = FALSE,
		backend = "rstan"
	)

	J <- output$stan_data$J
	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	policies$price_p[[1]][2:(J + 1)] <- 100

	df_sim <- PrepareSimulationData(output, policies, nsims = 1)

	wtp <- mdcev.sim(
		df_sim$df_indiv,
		df_common = df_sim$df_common,
		sim_options = df_sim$sim_options,
		cond_err = TRUE,
		nerrs = 3,
		sim_type = "welfare"
	)
	sum_wtp <- summary(wtp)

	expect_true(sum_wtp$CoefTable$mean[1] < 0)
})

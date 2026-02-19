# Tests for the default cmdstanr backend.
# All tests are marked skip_on_cran() because they require cmdstan to be installed
# and compiled, which is not available in CRAN's check environment.

tol <- 0.5   # wider tolerance: cmdstanr and rstan optimizers may converge to slightly different points

data(data_rec, package = "rmdcev")

data_rec_small <- mdcev.data(data_rec, subset = id <= 100,
							 id.var = "id",
							 alt.var = "alt",
							 choice = "quant")

test_that("cmdstanr MLE hybrid0 smoke test", {
	skip_on_cran()

	output <- mdcev(~ alt - 1,
					data = data_rec_small,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					print_iterations = FALSE)
	# No backend argument -> default is cmdstanr
	expect_equal(output$algorithm, "MLE")
	expect_equal(output$model, "hybrid0")
	expect_true(is.numeric(output$log.likelihood))
	expect_true(output$log.likelihood < 0)

	# LL should be in the same ballpark as the rstan result
	expect_equal(output$log.likelihood, -2653.237031, tolerance = tol)
})

test_that("cmdstanr MLE gamma smoke test", {
	skip_on_cran()

	output <- mdcev(~ 0,
					data = data_rec_small,
					model = "gamma",
					algorithm = "MLE",
					print_iterations = FALSE)

	expect_equal(output$model, "gamma")
	expect_true(output$log.likelihood < 0)
	expect_equal(output$log.likelihood, -2628.01, tolerance = tol)
})

test_that("cmdstanr MLE hybrid0 coef and logLik", {
	skip_on_cran()

	output <- mdcev(~ alt - 1,
					data = data_rec_small,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					print_iterations = FALSE)

	# coef() returns a list of parameter arrays
	cf <- coef(output)
	expect_true(is.list(cf))
	expect_true(length(cf) > 0)

	ll <- logLik(output)
	expect_true(is.numeric(ll))
})

test_that("cmdstanr MLE hybrid0 + welfare simulation", {
	skip_on_cran()

	output <- mdcev(~ alt - 1,
					data = data_rec_small,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					std_errors = "mvn",
					n_draws = 5,
					print_iterations = FALSE)

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)

	# Zero price change: welfare should be ~0
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = 0, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)
	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = 0.1)
})

test_that("cmdstanr Bayes hybrid0 fixed smoke test", {
	skip_on_cran()

	data_rec_bayes <- mdcev.data(data_rec, subset = id <= 200,
								 id.var = "id",
								 alt.var = "alt",
								 choice = "quant")

	output <- mdcev(~ 0,
					data = data_rec_bayes,
					model = "hybrid0",
					algorithm = "Bayes",
					random_parameters = "fixed",
					fixed_scale1 = 0,
					print_iterations = FALSE,
					n_cores = 1,
					n_chains = 1,
					n_iterations = 10,
					show_stan_warnings = FALSE)
	# No backend argument -> default is cmdstanr
	expect_equal(output$algorithm, "Bayes")
	expect_equal(output$parms_info$n_vars$n_parms_total, 34)

	# PrepareSimulationData + welfare should work with cmdstanr draws
	policies <- CreateBlankPolicies(2, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)
	expect_true(is.list(df_sim))

	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = 1, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)
	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = 0.1)
})

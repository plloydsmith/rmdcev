tol <- 0.01

data(data_rec, package = "rmdcev")

data_rec <- mdcev.data(data_rec, subset = id <= 200,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

# Fit LC model once; shared by several tests in this file.
result_test <- mdcev(~ alt | university,
					 data = data_rec,
					 model = "hybrid0",
					 psi_ascs = FALSE,
					 algorithm = "MLE",
					 n_classes = 2,
					 print_iterations = FALSE,
					 backend = "rstan")

test_that("LC 2-classes", {
	expect_true(result_test$log.likelihood < 0)
	expect_snapshot_value(round(result_test$log.likelihood, 1), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(result_test[["stan_fit"]][["par"]][["scale"]][[1]], 2), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(result_test[["stan_fit"]][["par"]][["delta"]][[1, 2]], 1), style = "deparse", cran = FALSE)
})

test_that("LC 2-classes with starting values", {

	result_start <- mdcev(~ 0 | university,
						  data = data_rec,
						  model = "hybrid0",
						  algorithm = "MLE",
						  n_classes = 2,
						  print_iterations = FALSE,
						  backend = "rstan")

	result_final <- mdcev(~ 0 | university,
						  data = data_rec,
						  model = "hybrid0",
						  initial.parameters = result_start$stan_fit$par,
						  algorithm = "MLE",
						  n_classes = 2,
						  print_iterations = FALSE,
						  backend = "rstan")

	expect_true(result_final[["stan_fit"]][["par"]][["scale"]][[1]] > 0)
	expect_snapshot_value(round(result_final[["stan_fit"]][["par"]][["scale"]][[1]], 3), style = "deparse", cran = FALSE)
})

test_that("Test LC simulations", {
	npols <- 2
	policies <- CreateBlankPolicies(npols, result_test, price_change_only = TRUE)

	df_sim_c1 <- PrepareSimulationData(result_test, policies, nsims = 1, class = "class1")
	df_sim    <- PrepareSimulationData(result_test, policies, nsims = 1, class = "class2")

	# Test welfare (zero price change -> WTP â‰ˆ 0)
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = TRUE, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)
	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)

	demand <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
						sim_options = df_sim$sim_options,
						cond_err = TRUE, nerrs = 1, sim_type = "demand")

	expect_equal(sum(demand[[5]][[1]][1, -1]), sum(result_test$stan_data$quant_j[5, ]))
})

context("Test LC")
#library(rmdcev)
tol <- 0.01
data(data_rec, package = "rmdcev")
data_rec

data_rec <- mdcev.data(data_rec, subset = id <= 200,
				   id.var = "id",
				   alt.var = "alt",
				   choice = "quant")

result_test <- mdcev( ~ alt | university,
				   data = data_rec,
				   model = "hybrid0",
				   psi_ascs = 0,
				   algorithm = "MLE",
				   n_classes = 2,
				   print_iterations = FALSE,
				   backend = "rstan")

test_that("LC 2-classes", {
	print(result_test$log.likelihood, digits =10)
	expect_true(abs(result_test$log.likelihood - (-5073.717925)) < tol)
	expect_true(abs(result_test[["stan_fit"]][["par"]][["scale"]][[1]] - 0.8072385) < tol)
	expect_true(abs(result_test[["stan_fit"]][["par"]][["delta"]][[1,2]] - -0.9649197) < tol)
})


test_that("LC 2-classes with starting values", {

	result_test <- mdcev( ~ 0| university,
					 data = data_rec,
					 model = "hybrid0",
					 algorithm = "MLE",
					 n_classes = 2,
					 print_iterations = FALSE,
					 backend = "rstan")

	result_test <- mdcev( ~ 0 | university,
					 data = data_rec,
					 model = "hybrid0",
					 initial.parameters = result_test$stan_fit$par,
					 algorithm = "MLE",
					 n_classes = 2,
					 print_iterations = FALSE,
					 backend = "rstan")

	expect_true(abs(result_test[["stan_fit"]][["par"]][["scale"]][[1]] - 0.7593204) < tol)

})

context("Test LC simulations")

test_that("Test LC simulations", {
	npols <- 2
	policies <- CreateBlankPolicies(npols, result_test,
									price_change_only = TRUE)

	df_sim <- PrepareSimulationData(result_test, policies, nsims = 1, class = "class1")
	df_sim <- PrepareSimulationData(result_test, policies, nsims = 1, class = "class2")

	# Test welfare
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < tol)

	demand <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 1, sim_type = "demand")

	print(demand[[5]][[1]][1,-1], digits =10)
	expect_equal(sum(demand[[5]][[1]][1,-1]), sum(result_test$stan_data$quant_j[5,]))
})

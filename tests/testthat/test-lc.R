context("Test LC")
#library(rmdcev)
tol <- 0.01
data(data_rec, package = "rmdcev")
data_rec

data_rec <- mdcev.data(data_rec, subset = id < 501,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

result <- mdcev( ~ 1 | university + ageindex,
				   data = data_rec,
				   model = "hybrid0",
				   algorithm = "MLE",
				   std_error = "deltamethod",
				   n_classes = 2,
				   print_iterations = FALSE)

test_that("LC 2-classes", {
#	expect_true(abs(result$log.likelihood - (-12612.51)) < tol)

	print(result$log.likelihood, digits =10)
#	expect_true(abs(result$bic - 25479.74) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]][[1]] - 0.682) < tol)
#	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][2,1] - -7.493115) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["delta"]][[1,2]] - .431) < tol)
})

context("Test LC simulations")

test_that("Test LC simulations", {
	npols <- 2
	policies <- CreateBlankPolicies(npols, result$stan_data[["J"]], result$stan_data[["dat_psi"]], price_change_only = TRUE)
	df_sim <- PrepareSimulationData(result, policies, nsims = 3, class = "class1")

	# Test welfare
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 3, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < tol)

	demand <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 3, sim_type = "demand")

#	print(demand[["class2"]][[5]][[2]][1,-1], digits =10)
	expect_equal(sum(demand[[5]][[2]][1,-1]), sum(result$stan_data$quant_j[5,]))
})


context("Test LC")
#library(rmdcev)
tol <- 0.01
data(data_rec, package = "rmdcev")
data_rec

result <- FitMDCEV(psi_formula = ~ 1,
				   lc_formula = ~ university + ageindex,
				   data = subset(data_rec, id < 500),
				   model = "hybrid0",
				   algorithm = "MLE",
				   std_error = "deltamethod",
				   n_classes = 2,
				   print_iterations = FALSE)

output.sum <- SummaryMDCEV(result)

test_that("LC 2-classes", {
#	expect_true(abs(result$log.likelihood - (-12612.51)) < tol)

	print(result$log.likelihood, digits =10)
#	expect_true(abs(result$bic - 25479.74) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]][[1]] - 0.6832412) < tol)
#	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][2,1] - -7.493115) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["delta"]][[1,2]] - .4318349) < tol)
})


 #test_that("Test LC simulations", {
#	npols <- 2
#	policies <-	CreateBlankPolicies(npols, result$stan_data[["J"]], result$stan_data[["dat_psi"]])
#	df_sim <- PrepareSimulationData(result, policies, nsims = 3)

	# Test welfare
#	wtp <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
#						 cond_err = 1, nerrs = 3, sim_type = "welfare")
#	sum_wtp <- purrr:::map(wtp, SummaryWelfare)
#	print(sum_wtp[["class1"]][["Mean"]], digits =10)

#	expect_true(sum(abs(sum_wtp[["class1"]][["Mean"]]), abs(sum_wtp[["class2"]][["Mean"]])) < tol)

#	demand <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
#						 cond_err = 1, nerrs = 3, sim_type = "demand")
#	print(demand[["class2"]][[5]][[2]][1,-1], digits =10)
#	expect_equal(sum(demand[["class2"]][[5]][[2]][1,-1]), sum(result$stan_data$j_quant[5,]))
#})


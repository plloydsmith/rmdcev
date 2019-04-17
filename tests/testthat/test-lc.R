context("Test LC")
library(rmdcev)

tol <- 0.00001
data(data_rec, package = "rmdcev")
data_rec

result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
				   lc_formula = ~ university + ageindex,
				   data = subset(data_rec, id < 500),
				   model = "gamma0",
				   algorithm = "MLE",
					n_classes = 2,
					print_iterations = FALSE)

test_that("LC 2-classes", {
	expect_true(abs(result$log.likelihood - (-12032.0652195997)) < tol)
	expect_true(abs(result$bic - 24517.6506841892) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]][[1]] - 0.8112873253) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][2,2] - -8.1932172861) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["beta_m"]][[1,2]] - -1.3655889455) < tol)
})


test_that("Test LC simulations", {
	npols <- 2
	policies <-	CreateBlankPolicies(npols, result$stan_data[["J"]], result$stan_data[["dat_psi"]])
	df_sim <- PrepareSimulationData(result, policies, nsims = 3)

	# Test welfare
	wtp <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 3, sim_type = "welfare")
	sum_wtp <- purrr:::map(wtp, SummaryWelfare)
	expect_true(sum(abs(sum_wtp[["class1"]][["Mean"]]), abs(sum_wtp[["class2"]][["Mean"]])) < 1e-3)

	demand <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 3, sim_type = "demand")
	expect_equal(sum(demand[["class2"]][[5]][[2]][1,-1]), sum(result$stan_data$j_quant[5,]))
})


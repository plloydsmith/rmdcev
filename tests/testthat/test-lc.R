context("Test LC")

library(pacman)

p_load(tidyverse, rmdcev, rstan)

tol <- 0.00001
data(data_rec, package = "rmdcev")
data_rec
result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
				   lc_formula = ~ university,
				   data = data_rec,
				   model = "gamma0",
				   algorithm = "MLE",
					n_classes = 2)

test_that("LC 2-classes", {

	expect_true(abs(result$log.likelihood - (-23054.814483864)) < tol)
	expect_true(abs(result$bic - 46606.771023165) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]][[1]] - 0.7110303) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][2,2] - -8.014577) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["gamma"]][[1,10]] - 3.238039) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["beta_m"]][[1,2]] - 0.04420978) < tol)
})


test_that("Test LC simulations", {
	npols <- 2
	policies <-	CreateBlankPolicies(npols, result$stan_data[["J"]], result$stan_data[["dat_psi"]])
	df_sim <- PrepareSimulationData(result, policies, nsims = 3)

	# Test welfare
	wtp <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 3, sim_type = "welfare")
	sum_wtp <- map(wtp, SummaryWelfare)
	expect_true(sum(abs(sum_wtp[["class1"]][["mean"]]), abs(sum_wtp[["class2"]][["mean"]])) < 1e-3)

	demand <- SimulateMDCEV(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
						 cond_err = 1, nerrs = 3, sim_type = "demand")
	expect_equal(sum(demand[["class2"]][[1]][[2]][1,-1]), sum(result$stan_data$j_quant[1,]))
})


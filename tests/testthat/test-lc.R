context("Test LC")
#Sys.setenv("R_TESTS" = "")
#data(eggs, package = "rmdcev")
library(pacman)

p_load(tidyverse, rstan, rmdcev)

tol <- 0.00001

data(recreation, package = "rmdcev")

result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
				   lc_formula = ~ university,
				   data = data,
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
	df_wtp <- PrepareSimulationData(result, policies, nsims = 3)

	# Test welfare
	wtp <- SimulateWTP(df_wtp$df_indiv, df_common = df_wtp$df_common, sim_options = df_wtp$sim_options, cond_err = 1, nerrs = 3)
	sum_wtp <- map(wtp, SummaryWelfare)
	expect_true(sum(abs(sum_wtp[["class1"]][["mean"]]), abs(sum_wtp[["class2"]][["mean"]])) < 1e-3)

})


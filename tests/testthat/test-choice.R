context("Test MLE")
#Sys.setenv("R_TESTS" = "")
#data(eggs, package = "rmdcev")
library(pacman)

p_load(tidyverse, rstan, rmdcev)

tol <- 0.00001

data(recreation, package = "rmdcev")

test_that("MLE gamma0 specification", {
	result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
								data = data,
								model = "gamma0",
								algorithm = "MLE")

	expect_true(abs(result$log.likelihood - (-23307.7292015328)) < tol)
	expect_true(abs(result$bic - 46857.124680015) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]] - 0.7765478) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][[1]] - -7.258057) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["gamma"]][[2]] - 19.78209) < tol)
})

test_that("MLE gamma specification", {
	result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
					   data = data,
					   model = "gamma",
					   algorithm = "MLE")

	expect_true(abs(result$log.likelihood - (-23299.218163508)) < tol)
	expect_true(abs(result$bic - 46847.007354736) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["alpha"]] - 0.1615235) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][[1]] - -5.497487) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["gamma"]][[2]] - 17.72717) < tol)
})

test_that("MLE alpha specification", {
	result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
					   data = data,
					   model = "alpha",
					   algorithm = "MLE")

	expect_true(abs(result$log.likelihood - (-24221.066687246)) < tol)
	expect_true(abs(result$bic - 48690.70440221) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]] - 0.6882879) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][[1]] - -2.091086) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["alpha"]][[2]] - 0.5536171) < tol)
})

test_that("MLE les specification", {
	result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
					   data = data,
					   model = "les",
					   algorithm = "MLE")

	expect_true(abs(result$log.likelihood - (-23168.772852057)) < tol)
	expect_true(abs(result$bic - 46586.116731833) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]] - 0.6528661) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["alpha"]] - 0.5654564) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["gamma"]][[2]] - 25.28034) < tol)
})

test_that("LC 2-classes", {
	result <- FitMDCEV(psi_formula = ~ factor(activity) -1,
					   lc_formula = ~ university,
					   data = data,
					   model = "gamma0",
					   algorithm = "MLE",
						n_classes = 2)

	expect_true(abs(result$log.likelihood - (-23054.814483864)) < tol)
	expect_true(abs(result$bic - 46606.771023165) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]][[1]] - 0.7110303) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][2,2] - -8.014577) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["gamma"]][[1,10]] - 3.238039) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["beta_m"]][[1,2]] - 0.04420978) < tol)
})

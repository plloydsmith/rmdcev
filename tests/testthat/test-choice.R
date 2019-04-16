context("Test MLE")
library(rmdcev)
library(rstan)
tol <- 0.00001
data(data_rec, package = "rmdcev")
data_rec
test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})

# sprintf("%.10f",result$log.likelihood)
# sprintf("%.10f",result$bic)
# sprintf("%.10f",result[["stan_fit"]][["par"]][["scale"]] )
# sprintf("%.10f",result[["stan_fit"]][["par"]][["psi"]] )
# sprintf("%.10f",result[["stan_fit"]][["par"]][["alpha"]] )
# sprintf("%.10f",wtp_err )


test_that("MLE gamma0 specification", {
	result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
								data = subset(data_rec, id < 1000),
								model = "gamma0",
								algorithm = "MLE",
					   print_iterations = FALSE)
#	SummaryMDCEV(result)
	expect_equal(result$model, "gamma0")
	expect_true(abs(result$log.likelihood - (-23483.6419731299)) < tol)
	expect_true(abs(result$bic - 47209.0203635126) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]] - 0.7602637661) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][[2]] - -8.3287213324) < tol)
	expect_equal(length(result[["stan_fit"]][["par"]][["alpha"]]), 0)
})

test_that("MLE gamma specification", {
	result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 1000),
					   model = "gamma",
					   algorithm = "MLE",
					   print_iterations = FALSE)

	expect_true(abs(result$log.likelihood - (-23468.7243361251)) < tol)
	expect_true(abs(result$bic - 47186.0918442815) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["alpha"]] - 0.1811457085) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][[1]] - -5.3074547474) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["gamma"]][[2]] - 19.4329273237) < tol)
})

test_that("MLE alpha specification", {
	result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 100),
					   model = "alpha",
					   algorithm = "MLE",
					   print_iterations = FALSE)

	expect_true(abs(result$log.likelihood - (-2706.8916415926)) < tol)
	expect_true(abs(result$bic - 5579.2075977901) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]] - 0.6872371856) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["psi"]][[1]] - -0.9909718259) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["alpha"]][[2]] - 0.5348307899) < tol)
	expect_equal(length(result[["stan_fit"]][["par"]][["gamma"]]), 0)
})

test_that("MLE les specification", {
	result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 100),
					   model = "les",
					   algorithm = "MLE",
					   print_iterations = FALSE)

	expect_true(abs(result$log.likelihood - (-2571.1107609921)) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["scale"]] - 0.6104077925) < tol)
	expect_true(abs(result[["stan_fit"]][["par"]][["alpha"]] - 0.6611009275) < tol)
})

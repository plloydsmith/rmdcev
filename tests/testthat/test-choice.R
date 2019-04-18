

library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec

context("Test Data load")

test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})

 sprintf("%.10f",output$log.likelihood)
# sprintf("%.10f",output$bic)
# sprintf("%.10f",output[["stan_fit"]][["par"]][["scale"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["psi"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["alpha"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["gamma"]] )

#skip_on_cran(
#

test_that("MLE names", {
	expect_error(FitMDCEV(psi_formula = ~ factor(good_name) -1,
									 data = data_rec,
									 model = "gamma77",
									 algorithm = "MLE",
									 print_iterations = FALSE))
})


context("MLE gamma0 specification")

test_that("MLE gamma0", {
	output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 100),
					   model = "gamma0",
					   algorithm = "MLE",
					   print_iterations = FALSE)

	output.sum <- SummaryMDCEV(output)
	expect_equal(length(output.sum$Std.err), 35)
	expect_equal(output$model, "gamma0")
	print(output$log.likelihood, digits =10)

	expect_true(abs(output$log.likelihood - (-2596.84))< tol)
	expect_true(abs(output$bic - 5354.50) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["scale"]] - 0.78) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["psi"]][[2]] - -8.05) < tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)
})


context("MLE gamma specification")


test_that("MLE alpha", {
	output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 100),
					   model = "alpha",
					   algorithm = "MLE",
					   print_iterations = FALSE)
	print(output$log.likelihood, digits =10)

	expect_true(abs(output$log.likelihood - (-2706.89)) < tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["gamma"]]), 0)
})

test_that("MLE les", {
	output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 100),
					   model = "les",
					   algorithm = "MLE",
					   print_iterations = FALSE)
	print(output$log.likelihood, digits =10)

	expect_true(abs(output$log.likelihood - (-2571.11)) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["scale"]] - 0.61) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["alpha"]] - 0.66) < tol)
})

test_that("MLE gamma", {
	output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 100),
					   model = "gamma",
					   algorithm = "MLE",
					   print_iterations = FALSE)

	expect_true(abs(output$log.likelihood - (-2596.837)) < tol)
	#	expect_true(abs(output[["stan_fit"]][["par"]][["alpha"]] - 0.15) < tol)
})


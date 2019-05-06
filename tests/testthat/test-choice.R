context("Test Data load")


library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec


test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})

# sprintf("%.10f",output$log.likelihood)
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
	output <- FitMDCEV(psi_formula = ~ 1,
					   data = subset(data_rec, id < 100),
					   model = "gamma0",
					   algorithm = "MLE",
					   print_iterations = FALSE)

	output.sum <- SummaryMDCEV(output)
	expect_equal(length(output.sum$Std.err), 19)
	expect_equal(output$model, "gamma0")
	print(output$log.likelihood, digits =10)

	expect_true(abs(output$log.likelihood - (-2684.757401))< tol)
	expect_true(abs(output$bic - 5456.822) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["scale"]] - 0.8592168) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["psi"]][[1]] - -7.635813) < tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)
})


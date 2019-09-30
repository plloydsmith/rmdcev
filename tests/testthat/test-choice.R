context("Test Data load")

tol <- 0.1

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

data_rec <- mdcev.data(data_rec, subset = id < 100,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

test_that("MLE names", {
	expect_error(mdcev( ~ alt -1,
						 data = data_rec,
						 model = "gamma77",
						 algorithm = "MLE",
						 print_iterations = FALSE))
})

context("MLE hybrid0 specification")

test_that("MLE hybrid0", {

	output <- mdcev( ~ 1,
				   data = data_rec,
				   model = "hybrid0",
				   algorithm = "MLE",
				   std_error = "deltamethod",
				   print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 19)
	expect_equal(output$model, "hybrid0")
	print(output$log.likelihood, digits =10)

	expect_true(abs(output$log.likelihood - (-2684.767401)) < tol)
	expect_true(abs(output$bic - 5456.841) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["scale"]] - 0.8592168) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["psi"]][[1]] - -7.635813) < tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)
})


context("MLE hybrid specification")

test_that("MLE hybrid", {
	output <- mdcev( ~ 1,
					data = data_rec,
					model = "hybrid",
				    algorithm = "MLE",
				    std_error = "deltamethod",
				    print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 20)
})

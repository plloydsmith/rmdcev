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

data_rec <- mdcev.data(data_rec, subset = id <= 100,
					  id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

test_that("MLE names", {
	expect_error(mdcev( ~ 0,
						 data = data_rec,
						 model = "gamma77",
						 algorithm = "MLE",
						 print_iterations = FALSE))
})

context("MLE hybrid0 specification")

test_that("MLE hybrid0", {

	output <- mdcev( ~ alt - 1,
				   data = data_rec,
				   model = "hybrid0",
					 psi_ascs = 0,
				   algorithm = "MLE",
				   print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_equal(output$model, "hybrid0")
	print(output$log.likelihood, digits =10)
	print(output$bic, digits =10)

	expect_true(abs(output$log.likelihood - (-2653.222965)) < tol)
	expect_true(abs(output$bic - 5467.626887) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["scale"]] - 0.7849449) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["psi"]][[1]] - -7.08514) < tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)
})


context("MLE hybrid specification")

test_that("MLE hybrid", {
	output <- mdcev( ~ 0,
					data = data_rec,
					model = "hybrid",
				    algorithm = "MLE",
				    print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
})


context("MLE gamma specification")

test_that("MLE gamma", {
	output <- mdcev( ~ 0,
					 data = data_rec,
					 model = "gamma",
					 algorithm = "MLE",
					 print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_true(abs(output$log.likelihood - (-2628.01)) < tol)
})


context("MLE alpha specification")

test_that("MLE alpha", {
	output <- mdcev( ~ 0,
					 data = data_rec,
					 model = "alpha",
					 algorithm = "MLE",
					 print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_true(abs(output$log.likelihood - (-2764.29)) < tol)
})

context("MLE kt_ee specification")

test_that("MLE kt_ee", {
	output <- mdcev( ~ 0 | 0 | 0,
					 data = data_rec,
					 model = "kt_ee",
					 algorithm = "MLE",
					 print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 19)
	expect_true(abs(output$log.likelihood - (-2775.93)) < tol)
})

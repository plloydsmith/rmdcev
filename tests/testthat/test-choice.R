context("Test Data load")

tol <- 0.1

data(data_rec, package = "rmdcev")
data_rec

test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})

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

	expect_true(abs(output$log.likelihood - (-2653.237031)) < tol)
	expect_true(abs(output$bic - 5467.655) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["scale"]] - 0.7856681) < tol)
	expect_true(abs(output[["stan_fit"]][["par"]][["psi"]][[1]] - -7.096204) < tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)
})


test_that("MLE hybrid0 mvn draws", {

	output <- mdcev( ~ alt - 1,
					 data = data_rec,
					 model = "hybrid0",
					 psi_ascs = 0,
					 algorithm = "MLE",
					 std_errors = "mvn",
					 print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_equal(output$model, "hybrid0")

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

	data_rec$beach = ifelse(data_rec$alt == "beach", 1, 0)

	output <- mdcev( ~ ageindex | 0 | beach,
					 data = data_rec,
					 gamma_ascs = 0,
					 model = "kt_ee",
					 algorithm = "MLE",
					 print_iterations = FALSE)

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 5)
	expect_true(abs(output$log.likelihood - (-2826.694653)) < tol)
})

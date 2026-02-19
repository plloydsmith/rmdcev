tol <- 0.1

data(data_rec, package = "rmdcev")

test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})

data_rec <- mdcev.data(data_rec, subset = id <= 100,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

test_that("MLE names error", {
	expect_error(mdcev(~ 0,
					   data = data_rec,
					   model = "gamma77",
					   algorithm = "MLE",
					   print_iterations = FALSE,
					   backend = "rstan"),
				 regexp = "model")
})

test_that("MLE hybrid0", {

	output <- mdcev(~ alt - 1,
				   data = data_rec,
				   model = "hybrid0",
				   psi_ascs = 0,
				   algorithm = "MLE",
				   print_iterations = FALSE,
				   backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_equal(output$model, "hybrid0")

	expect_equal(output$log.likelihood, -2653.237031, tolerance = tol)
	expect_equal(output$bic, 5467.655, tolerance = tol)
	expect_equal(as.numeric(output[["stan_fit"]][["par"]][["scale"]]), 0.7856681, tolerance = tol)
	expect_equal(as.numeric(output[["stan_fit"]][["par"]][["psi"]][[1]]), -7.096204, tolerance = tol)
	expect_equal(length(output[["stan_fit"]][["par"]][["alpha"]]), 0)

	# S3 method coverage: coef() returns a list of parameter arrays
	cf <- coef(output)
	expect_true(is.list(cf))
	expect_true(length(cf) > 0)

	ll <- logLik(output)
	expect_true(is.numeric(ll))

	# print does not error
	expect_output(print(output))
})

test_that("MLE hybrid0 mvn draws", {

	output <- mdcev(~ alt - 1,
					data = data_rec,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					std_errors = "mvn",
					print_iterations = FALSE,
					backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_equal(output$model, "hybrid0")
})

test_that("MLE hybrid", {
	output <- mdcev(~ 0,
					data = data_rec,
					model = "hybrid",
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
})

test_that("MLE gamma", {
	output <- mdcev(~ 0,
					data = data_rec,
					model = "gamma",
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_equal(output$log.likelihood, -2628.01, tolerance = tol)
})

test_that("MLE alpha", {
	output <- mdcev(~ 0,
					data = data_rec,
					model = "alpha",
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_equal(output$log.likelihood, -2764.29, tolerance = tol)
})

test_that("MLE kt_ee", {
	data_rec$beach <- ifelse(data_rec$alt == "beach", 1, 0)

	output <- mdcev(~ ageindex | 0 | beach,
					data = data_rec,
					gamma_ascs = 0,
					model = "kt_ee",
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 5)
	expect_equal(output$log.likelihood, -2826.694653, tolerance = tol)
})

test_that("MLE hybrid0 trunc_data", {
	output <- mdcev(~ alt - 1,
					data = data_rec,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					trunc_data = 1,
					print_iterations = FALSE,
					backend = "rstan")

	expect_equal(output$model, "hybrid0")
	expect_true(is.numeric(output$log.likelihood))
	expect_true(output$log.likelihood < 0)
})

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
				   psi_ascs = FALSE,
				   algorithm = "MLE",
				   print_iterations = FALSE,
				   backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 35)
	expect_equal(output$model, "hybrid0")

	expect_true(output$log.likelihood < 0)
	expect_snapshot_value(round(output$log.likelihood, 2), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(output$bic, 2), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(as.numeric(output[["stan_fit"]][["par"]][["scale"]]), 4), style = "deparse", cran = FALSE)
	expect_snapshot_value(round(as.numeric(output[["stan_fit"]][["par"]][["psi"]][[1]]), 4), style = "deparse", cran = FALSE)
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
					psi_ascs = FALSE,
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
	expect_true(output$log.likelihood < 0)
	expect_snapshot_value(round(output$log.likelihood, 2), style = "deparse", cran = FALSE)
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
	expect_true(output$log.likelihood < 0)
	expect_snapshot_value(round(output$log.likelihood, 2), style = "deparse", cran = FALSE)
})

test_that("MLE kt_ee", {
	data_rec$beach <- ifelse(data_rec$alt == "beach", 1, 0)

	output <- mdcev(~ ageindex | 0 | beach,
					data = data_rec,
					gamma_ascs = FALSE,
					model = "kt_ee",
					algorithm = "MLE",
					print_iterations = FALSE,
					backend = "rstan")

	output.sum <- summary(output)
	expect_equal(length(output.sum[["CoefTable"]]$Std.err), 5)
	expect_true(output$log.likelihood < 0)
	expect_snapshot_value(round(output$log.likelihood, 2), style = "deparse", cran = FALSE)
})

test_that("MLE hybrid0 trunc_data", {
	output <- mdcev(~ alt - 1,
					data = data_rec,
					model = "hybrid0",
					psi_ascs = FALSE,
					algorithm = "MLE",
					trunc_data = TRUE,
					print_iterations = FALSE,
					backend = "rstan")

	expect_equal(output$model, "hybrid0")
	expect_true(is.numeric(output$log.likelihood))
	expect_true(output$log.likelihood < 0)
})

test_that("MLE gamma parameter recovery from synthetic data", {
	skip_on_cran()
	set.seed(42)

	true_gamma <- c(2, 4, 6, 8, 10)

	sim_data <- GenerateMDCEVData(
		model       = "gamma",
		nobs        = 500,
		nalts       = 5,
		gamma_parms = true_gamma,
		scale_parms = 1.0
	)

	output <- mdcev(
		formula          = ~ b1 + b2 + b3 + b4 + b5 + b6,
		data             = sim_data$data,
		psi_ascs         = 0,
		model            = "gamma",
		algorithm        = "MLE",
		print_iterations = FALSE,
		backend          = "rstan"
	)

	est_gamma <- as.numeric(output[["stan_fit"]][["par"]][["gamma"]])

	expect_true(all(est_gamma > 0), label = "all gamma estimates positive")
	expect_true(est_gamma[5] > est_gamma[1],
				label = "ordering preserved: gamma[5]=10 > gamma[1]=2")
	expect_true(is.finite(output$log.likelihood))
	expect_true(output$log.likelihood < 0)
})

test_that("MLE hybrid0 parameter recovery from synthetic data", {
	skip_on_cran()
	set.seed(123)

	true_psi   <- c(-1.5, 2, -1)
	true_gamma <- rep(5, 5)

	sim_data <- GenerateMDCEVData(
		model       = "hybrid0",
		nobs        = 500,
		nalts       = 5,
		psi_i_parms = true_psi,
		gamma_parms = true_gamma,
		scale_parms = 1.0
	)

	output <- mdcev(
		formula          = ~ b1 + b2 + b3 + b4 + b5 + b6,
		data             = sim_data$data,
		psi_ascs         = 0,
		model            = "hybrid0",
		algorithm        = "MLE",
		print_iterations = FALSE,
		backend          = "rstan"
	)

	est_psi   <- as.numeric(output[["stan_fit"]][["par"]][["psi"]])
	est_gamma <- as.numeric(output[["stan_fit"]][["par"]][["gamma"]])

	# Signs of individual-level psi covariates should match true values
	expect_true(est_psi[length(est_psi) - 2] < 0, label = "first ind-level psi negative")
	expect_true(est_psi[length(est_psi) - 1] > 0, label = "second ind-level psi positive")
	expect_true(all(est_gamma > 0),                label = "all gamma positive")
	expect_true(is.finite(output$log.likelihood))
	expect_true(output$log.likelihood < 0)
})

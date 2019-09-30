context("Test Bayes fixed")

#library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")

data_rec <- mdcev.data(data_rec, subset = id < 500,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

output <- mdcev(~ alt -1,
				data = data_rec,
				model = "hybrid0",
				algorithm = "Bayes",
				random_parameters = "fixed",
				fixed_scale1 = 1,
				print_iterations = FALSE,
				n_cores = 1,
				n_chains = 1,
				n_iterations = 10,
				show_stan_warnings = FALSE)

test_that("Bayes hybrid0 specification", {
	expect_equal(output$parms_info$n_vars$n_parms_total, 34)

	output <- summary(output)
	expect_equal(dim(output[["CoefTable"]]), c(34, 5))
	expect_equal(output[["CoefTable"]]$z.stat[1], Inf)
})

context("Test Bayes rp uncorrelated")

output <- mdcev(formula = ~ alt -1,
				data = data_rec,
				model = "gamma",
				   algorithm = "Bayes",
				   random_parameters = "uncorr",
				   print_iterations = FALSE,
				   n_cores = 1,
				   n_chains = 1,
				   n_iterations = 10,
				   show_stan_warnings = FALSE)

test_that("Bayes gamma uncorr specification", {
	expect_equal(output$parms_info$n_vars$n_parms_total, 71)
})


context("Test Bayes rp correlated")

output <- mdcev(formula = ~ alt -1,
				data = data_rec,
				model = "gamma",
				   algorithm = "Bayes",
				   random_parameters = "corr",
				   print_iterations = FALSE,
				   n_cores = 1,
				   n_chains = 1,
				   n_iterations = 10,
				   show_stan_warnings = FALSE)

test_that("Bayes gamma uncorr specification", {
	expect_equal(output$parms_info$n_vars$n_parms_total, 666)
})

context("Test Bayes rp correlated with fixed gamma/alpha")

output <- mdcev(formula = ~ alt -1,
				data = data_rec,
				model = "gamma",
				   algorithm = "Bayes",
				   random_parameters = "corr",
				   print_iterations = FALSE,
				   gamma_fixed = 1,
				   alpha_fixed = 1,
				   fixed_scale1 = 1,
				   n_cores = 1,
				   n_chains = 1,
				   n_iterations = 10,
				   show_stan_warnings = FALSE)

test_that("Bayes gamma uncorr specification", {
	expect_equal(output$parms_info$n_vars$n_parms_total, 188)
})


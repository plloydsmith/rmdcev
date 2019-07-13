context("Test Bayes fixed")

#library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec

output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
							data = subset(data_rec, id < 201),
							model = "hybrid0",
							algorithm = "Bayes",
						random_parameters = "uncorr",
						fixed_scale1 = 1,
				   print_iterations = FALSE,
				   n_cores = 1,
				   n_chains = 1,
				   n_iterations = 200,
				   show_stan_warnings = FALSE)

test_that("Bayes hybrid0 specification", {
	expect_equal(dim(output$est_pars), c(175, 3))
	expect_equal(output$parms_info$n_vars$n_parms_total, 35)

	output <- SummaryMDCEV(output)
	expect_equal(output$z.stat[1], Inf)
})

context("Test Bayes rp uncorrelated")

output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
				   data = subset(data_rec, id < 500),
				   model = "gamma",
				   algorithm = "Bayes",
				   random_parameters = "uncorr",
				   print_iterations = FALSE,
				   n_cores = 1,
				   n_chains = 1,
				   n_iterations = 10,
				   show_stan_warnings = FALSE)

test_that("Bayes gamma uncorr specification", {
	expect_equal(dim(output$est_pars), c(93805, 3))
	expect_equal(output$parms_info$n_vars$n_parms_total, 71)
})


context("Test Bayes rp correlated")

output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
				   data = subset(data_rec, id < 500),
				   model = "gamma",
				   algorithm = "Bayes",
				   random_parameters = "corr",
				   print_iterations = FALSE,
				   n_cores = 1,
				   n_chains = 1,
				   n_iterations = 10,
				   show_stan_warnings = FALSE)

test_that("Bayes gamma uncorr specification", {
	expect_equal(dim(output$est_pars), c(99930, 3))
	expect_equal(output$parms_info$n_vars$n_parms_total, 666)
})

context("Test Bayes rp correlated")

output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
				   data = subset(data_rec, id < 500),
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
	expect_equal(dim(output$est_pars), c(45565, 3))
	expect_equal(output$parms_info$n_vars$n_parms_total, 188)
})


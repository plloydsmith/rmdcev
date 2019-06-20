library(rmdcev)
context("Test Bayes fixed")

#library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec

output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
							data = subset(data_rec, id < 500),
							model = "hybrid0",
							algorithm = "Bayes",
						random_parameters = "fixed",
				   print_iterations = FALSE,
				   n_cores = 1,
				   n_chains = 1,
				   n_iterations = 10,
				   show_stan_warnings = FALSE)

test_that("Bayes hybrid0 specification", {
	expect_equal(dim(output$est_pars), c(175, 3))
	expect_equal(output$parms_info$n_vars$n_parms_total, 35)

	output <- SummaryMDCEV(output)
	expect_equal(output$z.stat[1], Inf)
})


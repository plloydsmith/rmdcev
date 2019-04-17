context("Test Bayes")

tol <- 0.00001
data(data_rec, package = "rmdcev")
data_rec

# sprintf("%.10f",result$log.likelihood)
# sprintf("%.10f",result$bic)
# sprintf("%.10f",result[["stan_fit"]][["par"]][["scale"]] )
# sprintf("%.10f",result[["stan_fit"]][["par"]][["psi"]] )
# sprintf("%.10f",result[["stan_fit"]][["par"]][["alpha"]] )
# sprintf("%.10f",wtp_err )


test_that("Bayes gamma0 specification", {

	result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
								data = subset(data_rec, id < 1000),
						#	data = data_rec,
								model = "gamma0",
								algorithm = "Bayes",
							random_parameters = "fixed",
					   print_iterations = FALSE,
					   n_chains = 1,
					   n_iterations = 10,
					   show_stan_warnings = FALSE)


	expect_equal(dim(result$est_pars), c(175, 3))

	output <- SummaryMDCEV(result)
	expect_equal(output$z.stat[1], Inf)
})

test_that("Bayes gamma0 specification with correlated RP", {

	result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 1000),
					   #	data = data_rec,
					   model = "gamma0",
					   algorithm = "Bayes",
					   random_parameters = "corr",
					   print_iterations = FALSE,
					   n_chains = 1,
					   n_iterations = 10,
					   show_stan_warnings = FALSE)


	expect_equal(dim(result$est_pars), c(175, 3))
})

test_that("Bayes gamma0 specification with uncorrelated RP", {

	result <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
					   data = subset(data_rec, id < 1000),
					   #	data = data_rec,
					   model = "gamma0",
					   algorithm = "Bayes",
					   random_parameters = "uncorr",
					   print_iterations = FALSE,
					   n_chains = 1,
					   n_iterations = 10,
					   show_stan_warnings = FALSE)


	expect_equal(dim(result$est_pars), c(175, 3))
})

context("Test Bayes fixed")

#library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec

# sprintf("%.10f",output$log.likelihood)
# sprintf("%.10f",output$bic)
# sprintf("%.10f",output[["stan_fit"]][["par"]][["scale"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["psi"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["alpha"]] )
# sprintf("%.10f",wtp_err )


	output <- FitMDCEV(psi_formula = ~ factor(good_name) -1,
								data = subset(data_rec, id < 500),
								model = "hybrid0",
								algorithm = "Bayes",
							random_parameters = "fixed",
					   print_iterations = FALSE,
					   n_chains = 1,
					   n_iterations = 10,
					   show_stan_warnings = FALSE)


	expect_equal(dim(output$est_pars), c(175, 3))
	expect_equal(output$parms_info$n_vars$n_parms_total, 35)

	output <- SummaryMDCEV(output)
test_that("Bayes hybrid0 specification", {
	expect_equal(output$z.stat[1], Inf)
})


context("Test Bayes fixed")

#library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")

data_rec <- mdcev.data(data_rec, subset = id <= 200,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

output <- mdcev(~ 0,
				data = data_rec,
				model = "hybrid0",
				algorithm = "Bayes",
				random_parameters = "fixed",
				fixed_scale1 = 0,
				print_iterations = FALSE,
				n_cores = 1,
				n_chains = 1,
				n_iterations = 10,
				show_stan_warnings = FALSE)

test_that("Bayes hybrid0 specification", {
	expect_equal(output$parms_info$n_vars$n_parms_total, 34)

	output_sum <- summary(output)
	expect_equal(dim(output_sum[["CoefTable"]]), c(34, 5))

	npols <- 2
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)

	# Test welfare
	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
					 cond_err = 1, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < tol)
})

test_that("Bayes gamma specification no gamma/psi ascs", {

	output <- mdcev(~ 0,
					data = data_rec,
					model = "gamma",
					algorithm = "Bayes",
					random_parameters = "fixed",
					psi_ascs = 0,
					gamma_ascs = 0,
					fixed_scale1 = 1,
					print_iterations = FALSE,
					n_cores = 1,
					n_chains = 1,
					n_iterations = 10,
					show_stan_warnings = FALSE)

	expect_equal(output$parms_info$n_vars$n_parms_total, 2)

	output_sum <- summary(output)
	expect_equal(dim(output_sum[["CoefTable"]]), c(2, 5))

	npols <- 2
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)

	# Test welfare

	df_indiv <- lapply(df_sim$df_indiv, "[", c(1:10))
	wtp <- mdcev.sim(df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
					 cond_err = 1, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < tol)
})

context("Test Bayes rp uncorrelated")

output <- mdcev(formula = ~ 0,
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
	expect_equal(output$parms_info$n_vars$n_parms_total, 69)

	output_sum <- summary(output)
	expect_equal(dim(output_sum[["CoefTable"]]), c(69, 5))

	npols <- 2
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)

	# Test welfare
	df_indiv <- lapply(df_sim$df_indiv, "[", c(1:10))
	wtp <- mdcev.sim(df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
					 cond_err = 1, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < tol)
})

context("Test Bayes rp correlated")

output <- mdcev(formula = ~ 0,
				data = data_rec,
				model = "gamma",
				   algorithm = "Bayes",
				   random_parameters = "corr",
				   print_iterations = FALSE,
				   n_cores = 1,
				   n_chains = 1,
				 fixed_scale1 = 1,
				   n_iterations = 10,
				   show_stan_warnings = FALSE)

test_that("Bayes gamma corr specification", {

	output_sum <- summary(output)
	expect_equal(output$parms_info$n_vars$n_parms_total, 629)

	npols <- 2
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)

	# Test welfare
	df_indiv <- lapply(df_sim$df_indiv, "[", c(1:10))
	wtp <- mdcev.sim(df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
					 cond_err = 1, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < tol)
})

context("Test Bayes rp correlated with fixed gamma/alpha")

output <- mdcev(formula = ~ 0,
				data = data_rec,
				model = "gamma",
				algorithm = "Bayes",
				random_parameters = "corr",
				print_iterations = FALSE,
				gamma_nonrandom = 1,
				alpha_nonrandom = 1,
				n_cores = 1,
				n_chains = 1,
				n_iterations = 10,
				show_stan_warnings = FALSE)

test_that("Bayes gamma corr specification with fixed gamma/alpha", {
	output_sum <- summary(output)
	expect_equal(output$parms_info$n_vars$n_parms_total, 171)

	npols <- 2
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)

	# Test welfare
	df_indiv <- lapply(df_sim$df_indiv, "[", c(1:10))
	wtp <- mdcev.sim(df_indiv, df_common = df_sim$df_common, sim_options = df_sim$sim_options,
					 cond_err = 1, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_true(sum(abs(sum_wtp$CoefTable$mean)) < tol)
})

context("Test Bayes kt_ee rp")

output <- mdcev(formula = ~ ageindex|0|0,
				data = data_rec,
				model = "kt_ee",
				algorithm = "Bayes",
				random_parameters = "fixed",
				print_iterations = FALSE,
				gamma_nonrandom = 1,
				alpha_nonrandom = 1,
				fixed_scale1 = 1,
				n_cores = 1,
				n_chains = 1,
				n_iterations = 10,
				show_stan_warnings = FALSE)

test_that("Test Bayes kt_ee rp fixed", {
	output_sum <- summary(output)
	expect_equal(output$parms_info$n_vars$n_parms_total, 19)

	npols <- 2
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)
})

test_that("Test Bayes kt_ee rp uncorrelated", {

	output <- mdcev(formula = ~ ageindex |0|0,
					data = data_rec,
					model = "kt_ee",
					algorithm = "Bayes",
					random_parameters = "uncorr",
					print_iterations = FALSE,
					gamma_ascs = 0,
					fixed_scale1 = 1,
					n_cores = 1,
					n_chains = 1,
					n_iterations = 10,
					show_stan_warnings = FALSE)

		output_sum <- summary(output)
	expect_equal(output$parms_info$n_vars$n_parms_total, 6)

	npols <- 2
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 5)
})

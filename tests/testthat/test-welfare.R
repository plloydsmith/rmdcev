tol <- 0.1

data(data_rec, package = "rmdcev")

data_rec <- mdcev.data(data_rec, subset = id <= 100,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

test_that("MLE hybrid0 unconditional welfare", {

	output <- mdcev(~ alt - 1,
				   data = data_rec,
				   model = "hybrid0",
				   psi_ascs = 0,
				   algorithm = "MLE",
				   std_errors = "mvn",
				   n_draws = 2,
				   print_iterations = FALSE,
				   backend = "rstan")

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 2)

	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = 0, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)
})

test_that("MLE hybrid unconditional welfare", {

	output <- mdcev(~ alt - 1,
					data = data_rec,
					model = "hybrid",
					psi_ascs = 0,
					algorithm = "MLE",
					std_errors = "mvn",
					n_draws = 2,
					print_iterations = FALSE,
					backend = "rstan")

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 2)

	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = 0, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)
})

test_that("MLE gamma unconditional welfare", {
	output <- mdcev(~ alt - 1,
				 data = data_rec,
				 model = "gamma",
				 psi_ascs = 0,
				 algorithm = "MLE",
				 std_errors = "mvn",
				 n_draws = 2,
				 print_iterations = FALSE,
				 backend = "rstan")

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 2)

	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
				 sim_options = df_sim$sim_options,
				 cond_err = 0, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)
})

test_that("MLE alpha unconditional welfare", {

	output <- mdcev(~ alt - 1,
				 data = data_rec,
				 model = "alpha",
				 psi_ascs = 0,
				 algorithm = "MLE",
				 std_errors = "mvn",
				 n_draws = 2,
				 print_iterations = FALSE,
				 backend = "rstan")

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 2)

	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
				 sim_options = df_sim$sim_options,
				 cond_err = 0, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)
})

test_that("MLE kt_ee welfare", {

	data_rec$beach <- ifelse(data_rec$alt == "beach", 1, 0)

	output <- mdcev(~ ageindex | 0 | beach,
					 data = data_rec,
					 gamma_ascs = 0,
					 model = "kt_ee",
					 std_errors = "mvn",
					 n_draws = 30,
					 algorithm = "MLE",
					 print_iterations = FALSE,
					 backend = "rstan")

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 1)

	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = 0, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)
})

test_that("MLE hybrid0 deltamethod welfare", {

	output <- mdcev(~ alt - 1,
					data = data_rec,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					std_errors = "deltamethod",
					print_iterations = FALSE,
					backend = "rstan")

	npols <- 1
	policies <- CreateBlankPolicies(npols, output, price_change_only = TRUE)
	df_sim <- PrepareSimulationData(output, policies, nsims = 1)

	wtp <- mdcev.sim(df_sim$df_indiv, df_common = df_sim$df_common,
					 sim_options = df_sim$sim_options,
					 cond_err = 0, nerrs = 1, sim_type = "welfare")
	sum_wtp <- summary(wtp)

	expect_equal(sum(abs(sum_wtp$CoefTable$mean)), 0, tolerance = tol)
})

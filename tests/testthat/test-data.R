

library(rmdcev)

tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec

context("Test Data load")

test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})


test_that("ChechMdcev data works", {
	expect_error(CheckMdcevData(data_rec[,-1]), regexp = "id column")
})


# sprintf("%.10f",output$log.likelihood)
# sprintf("%.10f",output$bic)
# sprintf("%.10f",output[["stan_fit"]][["par"]][["scale"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["psi"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["alpha"]] )
# sprintf("%.10f",output[["stan_fit"]][["par"]][["gamma"]] )

#skip_on_cran(
#

test_that("MLE names", {
	expect_error(FitMDCEV(psi_formula = ~ factor(good_name) -1,
									 data = data_rec,
									 model = "gamma77",
									 algorithm = "MLE",
									 print_iterations = FALSE))
})


psi_formula = ~ factor(good_name) -1
lc_formula = NULL
weights = NULL
num_price = NULL
model = "gamma0"
n_classes = 1
fixed_scale = 0
trunc_data = 0
seed = "123"
initial.parameters = NULL
algorithm = "MLE"
flat_priors = NULL
print_iterations = TRUE
hessian = TRUE
prior_psi_sd = 10
prior_gamma_sd = 10
prior_alpha_sd = 0.5
prior_scale_sd = 1
prior_beta_m_sd = 10
n_draws = 30
keep_loglik = 0
random_parameters = "fixed"
show_stan_warnings = TRUE
n_iterations = 200
n_chains = 4
n_cores = 1
max_tree_depth = 10
adapt_delta = 0.8
lkj_shape_prior = 4

mle_options <- list(fixed_scale = fixed_scale,
					model = model,
					n_classes = n_classes,
					trunc_data = trunc_data,
					seed = seed,
					print_iterations = print_iterations,
					hessian = hessian,
					n_draws = n_draws,
					keep_loglik = keep_loglik,
					flat_priors = flat_priors,
					prior_psi_sd = prior_psi_sd,
					prior_gamma_sd = prior_gamma_sd,
					prior_alpha_sd = prior_alpha_sd,
					prior_scale_sd = prior_scale_sd,
					prior_beta_m_sd = prior_beta_m_sd)

bayes_options <- list(n_iterations = n_iterations,
					  n_chains = n_chains,
					  n_cores = 1,
					  keep_loglik = keep_loglik,
					  random_parameters = random_parameters,
					  seed = seed,
					  max_tree_depth = max_tree_depth,
					  adapt_delta = adapt_delta,
					  show_stan_warnings = show_stan_warnings,
					  lkj_shape_prior = lkj_shape_prior)

stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula, mle_options)
parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)

test_that("Gamma0 parm names", {
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

test_that("Gamma0 LC 2 parm names", {
	mle_options$n_classes = 2
	lc_formula = ~ university + ageindex
	stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula,  mle_options)
	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})
test_that("Gamma parm names", {
	mle_options$model = "gamma"
	stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula, mle_options)
	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})


test_that("alpha parm names", {
	mle_options$model = "alpha"
	stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula, mle_options)
	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

test_that("les parm names", {
	mle_options$model = "les"
	stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula, mle_options)
	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

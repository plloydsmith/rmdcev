context("Test Data load")
tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec


test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})


test_that("ChechMdcev data works", {
	expect_error(CheckMdcevData(data_rec[,-1]), regexp = "id column")
})


context("Test MLE names")

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
model = "hybrid0"
n_classes = 1
fixed_scale1 = 0
trunc_data = 0
seed = "123"
initial.parameters = NULL
algorithm = "MLE"
flat_priors = NULL
max_iterations = 500
print_iterations = TRUE
hessian = TRUE
prior_psi_sd = 10
prior_gamma_sd = 10
prior_alpha_sd = 0.5
prior_scale_sd = 1
prior_delta_sd = 10
gamma_fixed = 1
alpha_fixed = 1
n_draws = 30
keep_loglik = 0
random_parameters = "fixed"

mle_options <- list(fixed_scale1 = fixed_scale1,
					model = model,
					n_classes = n_classes,
					trunc_data = trunc_data,
					seed = seed,
					max_iterations = max_iterations,
					print_iterations = print_iterations,
					hessian = hessian,
					n_draws = n_draws,
					keep_loglik = keep_loglik,
					flat_priors = flat_priors,
					prior_psi_sd = prior_psi_sd,
					prior_gamma_sd = prior_gamma_sd,
					prior_alpha_sd = prior_alpha_sd,
					prior_scale_sd = prior_scale_sd,
					prior_delta_sd = prior_delta_sd,
					gamma_fixed = gamma_fixed,
					alpha_fixed = alpha_fixed)

stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula, mle_options)
parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)

test_that("hybrid0 parm names", {
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

test_that("hybrid0 LC 2 parm names", {
	mle_options$n_classes = 2
	lc_formula = ~ university + ageindex
	stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula,  mle_options)
	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})
test_that("hybrid parm names", {
	mle_options$model = "hybrid"
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

test_that("gamma parm names", {
	mle_options$model = "gamma"
	stan_data <- processMDCEVdata(data_rec, psi_formula, lc_formula, mle_options)
	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

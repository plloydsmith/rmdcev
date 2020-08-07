context("Test Data load")
tol <- 0.01

data(data_rec, package = "rmdcev")
data_rec

test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})

test_that("CheckMdcev data works", {
	expect_error(mdcev.datacheck(data_rec[,-1]), regexp = "id column")
})

context("Test MLE names")

test_that("MLE names", {
	expect_error(mdcev( ~ alt -1,
									 data = data_rec,
									 model = "gamma77",
									 algorithm = "MLE",
									 print_iterations = FALSE))
})

formula = ~ 0 | university + ageindex
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
gamma_nonrandom = 1
psi_ascs = NULL
gamma_ascs = 1
alpha_nonrandom = 1
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
					psi_ascs = psi_ascs,
					gamma_ascs = gamma_ascs,
					gamma_nonrandom = gamma_nonrandom,
					alpha_nonrandom = alpha_nonrandom)

data_rec <- mdcev.data(data_rec,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

alt_names <- unique(data_rec$alt)

stan_data <- processMDCEVdata(formula, data_rec, mle_options)
parms_info <- rmdcev:::CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)

test_that("hybrid0 parm names", {
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

test_that("hybrid0 LC 2 parm names", {
	mle_options$n_classes = 2
	stan_data <- processMDCEVdata(formula, data_rec, mle_options)
	parms_info <- CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})
test_that("hybrid parm names", {
	mle_options$model = "hybrid"
	stan_data <- processMDCEVdata(formula, data_rec, mle_options)
	parms_info <- CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})


test_that("alpha parm names", {
	mle_options$model = "alpha"
	stan_data <- processMDCEVdata(formula, data_rec, mle_options)
	parms_info <- CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

test_that("gamma parm names", {
	mle_options$model = "gamma"
	stan_data <- processMDCEVdata(formula, data_rec, mle_options)
	parms_info <- CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

test_that("kt parm names", {
	mle_options$model = "kt_les"
	formula = ~ alt -1 | 0 | 0
	stan_data <- processMDCEVdata(formula, data_rec, mle_options)
	parms_info <- CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

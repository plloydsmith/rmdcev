tol <- 0.01

data(data_rec, package = "rmdcev")

test_that("Data ok", {
	expect_equal(data_rec$id[18], 2)
})

test_that("MLE bad model name errors", {
	data_small <- mdcev.data(data_rec, subset = id <= 50,
							 id.var = "id", alt.var = "alt", choice = "quant")
	expect_error(
		mdcev(~ alt - 1,
			  data = data_small,
			  model = "gamma77",
			  algorithm = "MLE",
			  print_iterations = FALSE),
		regexp = "model"
	)
})

test_that("non-id names", {
	skip_on_cran()
	data_test <- data_rec %>%
		rename(id2 = id)

	data_test <- mdcev.data(data_test, subset = id2 <= 100,
							id.var = "id2",
							alt.var = "alt",
							choice = "quant")

	output <- mdcev(~ 0,
					data = data_test,
					model = "hybrid0",
					psi_ascs = 0,
					algorithm = "MLE",
					print_iterations = FALSE)

	output.sum <- summary(output)
	expect_s3_class(output, "mdcev")
})

# ── Module-level setup for CreateParmInfo tests ────────────────────────────
formula <- ~ 0 | university + ageindex
weights <- NULL
num_price <- NULL
model <- "hybrid0"
n_classes <- 1
fixed_scale1 <- 0
single_scale <- 0
trunc_data <- 0
seed <- "123"
initial.parameters <- NULL
algorithm <- "MLE"
flat_priors <- NULL
max_iterations <- 500
print_iterations <- TRUE
hessian <- TRUE
prior_psi_sd <- 10
prior_phi_sd <- 10
prior_gamma_sd <- 10
prior_alpha_shape <- 1
prior_scale_sd <- 1
prior_delta_sd <- 10
gamma_nonrandom <- 1
psi_ascs <- NULL
gamma_ascs <- 1
alpha_nonrandom <- 1
n_draws <- 30
keep_loglik <- 0
random_parameters <- "fixed"
jacobian_analytical_grad <- 1

mle_options <- list(seed = seed,
					max_iterations = max_iterations,
					hessian = hessian,
					print_iterations = print_iterations,
					n_draws = n_draws,
					keep_loglik = keep_loglik,
					n_classes = n_classes)

stan_model_options <- list(fixed_scale1 = fixed_scale1,
						   single_scale = single_scale,
						   model = model,
						   n_classes = n_classes,
						   trunc_data = trunc_data,
						   psi_ascs = psi_ascs,
						   gamma_ascs = gamma_ascs,
						   jacobian_analytical_grad = jacobian_analytical_grad,
						   flat_priors = flat_priors,
						   prior_psi_sd = prior_psi_sd,
						   pior_phi_sd = prior_phi_sd,
						   prior_gamma_sd = prior_gamma_sd,
						   prior_alpha_shape = prior_alpha_shape,
						   prior_scale_sd = prior_scale_sd,
						   prior_delta_sd = prior_delta_sd,
						   gamma_nonrandom = gamma_nonrandom,
						   alpha_nonrandom = alpha_nonrandom)

data_rec <- mdcev.data(data_rec,
					   id.var = "id",
					   alt.var = "alt",
					   choice = "quant")

alt_names <- unique(data_rec$alt)

stan_data <- processMDCEVdata(formula, data_rec, stan_model_options)
parms_info <- rmdcev:::CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)

# ── Parameter-count tests ─────────────────────────────────────────────────

test_that("hybrid0 parm names", {
	expect_equal(parms_info$n_vars$n_parms_total, length(parms_info$parm_names$all_names))
})

test_that("hybrid0 LC 2 parm names", {
	opts2 <- stan_model_options
	opts2$n_classes <- 2
	sd2 <- processMDCEVdata(formula, data_rec, opts2)
	pi2 <- CreateParmInfo(sd2, alt_names, algorithm, random_parameters)
	expect_equal(pi2$n_vars$n_parms_total, length(pi2$parm_names$all_names))
})

test_that("hybrid parm names", {
	opts <- stan_model_options
	opts$model <- "hybrid"
	sd <- processMDCEVdata(formula, data_rec, opts)
	pi <- CreateParmInfo(sd, alt_names, algorithm, random_parameters)
	expect_equal(pi$n_vars$n_parms_total, length(pi$parm_names$all_names))
})

test_that("alpha parm names", {
	opts <- stan_model_options
	opts$model <- "alpha"
	sd <- processMDCEVdata(formula, data_rec, opts)
	pi <- CreateParmInfo(sd, alt_names, algorithm, random_parameters)
	expect_equal(pi$n_vars$n_parms_total, length(pi$parm_names$all_names))
})

test_that("gamma parm names", {
	opts <- stan_model_options
	opts$model <- "gamma"
	sd <- processMDCEVdata(formula, data_rec, opts)
	pi <- CreateParmInfo(sd, alt_names, algorithm, random_parameters)
	expect_equal(pi$n_vars$n_parms_total, length(pi$parm_names$all_names))
})

test_that("kt_ee parm names", {
	opts <- stan_model_options
	opts$model <- "kt_ee"
	formula_ktee <- ~ alt - 1 | 0 | 0
	sd <- processMDCEVdata(formula_ktee, data_rec, opts)
	pi <- CreateParmInfo(sd, alt_names, algorithm, random_parameters)
	expect_equal(pi$n_vars$n_parms_total, length(pi$parm_names$all_names))
})

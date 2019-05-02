#' @title FitMDCEV
#' @description Fit a MDCEV model using MLE or Bayes
#' @param psi_formula Formula for psi
#' @param lc_formula Formula for latent class
#' @param data The (IxJ) data to be passed to Stan including 1) id, 2) good, 3) quant,
#' 4) price, 5) income, and columns for psi variables. Arrange data by id then good.
#' Notes I is number of individuals and J is number of non-numeraire goods.
#' @param weights An optional vector of sampling or frequency weights.
#' @param model A string indicating which model specification is estimated.
#' The options are "alpha","les", "gamma", and "gamma0".
#' @param n_classes The number of latent classes.
#' @param fixed_scale Whether to fix scale at 1.
#' @param trunc_data Whether the estimation should be adjusted for truncation
#' @param seed Random seed.
#' @param algorithm Either "Bayes" for Bayes or "MLE"
#'     for maximum liklihood estimation.
#' @param flat_priors indicator if completely uninformative priors should be specified. If using MLE, the
#' optimizing function will then be equal to log-likelihood. Defaults to 1 if MLE used and 0 if Bayes used.
#' @param print_iterations Whether to print iteration information
#' @param std_errors Compute standard errors using the delta method ("deltamethod")
#' or multivariate normal draws ("mvn"). The default is "mvn" as only mvn parameter draws are required
#' for demand and welfare simulation.
#' @param n_draws The number of MVN draws for standard error calculations
#' @param keep_loglik Whether to keep the log_lik calculations
#' @param hessian Wheter to keep the Hessian matrix
#' @param old test for speed
#' @param initial.parameters Specify initial parameters intead of
#'     starting at random.
#' @param prior_psi_sd standard deviation for normal prior with mean 0.
#' @param prior_gamma_sd standard deviation for normal prior with mean 0.
#' @param prior_alpha_sd standard deviation for normal prior with mean 0.5.
#' @param prior_scale_sd standard deviation for normal prior with mean 1.
#' @param prior_delta_sd standard deviation for normal prior with mean 0.
#' @param lkj_shape_prior Prior for Cholesky matrix
#' @param n_iterations The number of iterations in Bayesian estimation.
#' @param n_chains The number of chains in Bayesian estimation.
#' @param n_cores The number of cores to use in Bayesian estimation.
#' Can set using options(mc.cores = parallel::detectCores()).
#' @param random_parameters The form of the covariance matrix for
#'     Bayes. Can be 'fixed', 'uncorr, 'corr'.
#' @param max_tree_depth
#'     http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#' @param adapt_delta
#'    http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#' @param show_stan_warnings Whether to show warnings from Stan.
#' @param ... Additional parameters to pass on to \code{rstan::stan}
#'     and \code{rstan::sampling}.
#' @return A stanfit object
#' @export
#'
FitMDCEV <- function(data,
					 psi_formula = NULL,
					 lc_formula = NULL,
					 old = 0,
					 weights = NULL,
					 num_price = NULL,
					 model = c("alpha", "les", "gamma", "gamma0"),
					 n_classes = 1,
					 fixed_scale = 0,
					 trunc_data = 0,
					 seed = "123",
					 initial.parameters = NULL,
					 algorithm = c("MLE", "Bayes"),
					 #	std_errors = "draws", # still need to implement
					 flat_priors = NULL,
					 print_iterations = TRUE,
					 #mle_tol = 0.0001,
					 hessian = TRUE,
					 prior_psi_sd = 10,
					 prior_gamma_sd = 10,
					 prior_alpha_sd = 0.5,
					 prior_scale_sd = 1,
					 prior_delta_sd = 10,
					 std_errors = "mvn",
					 n_draws = 50,
					 keep_loglik = 0,
					 random_parameters = "fixed",
					 show_stan_warnings = TRUE,
					 n_iterations = 200,
					 n_chains = 4,
					 n_cores = 4,
					 max_tree_depth = 10,
					 adapt_delta = 0.8,
					 lkj_shape_prior = 4)
{
	CheckMdcevData(data)

	if (algorithm == "Bayes" && n_classes > 1)
		stop("Bayesian estimation can only be used with one class. Switch algorithm to MLE or choose n_classes = 1", "\n")

	if (algorithm == "MLE" && random_parameters != "fixed")
		stop("MLE can only be used with fixed parameters. Switch random_parameters to 'fixed' or change algorithm to Bayes", "\n")

	if (algorithm == "MLE" && is.null(flat_priors)){
		flat_priors = 1
	} else if (algorithm == "Bayes" && is.null(flat_priors))
		flat_priors = 0

	mle_options <- list(fixed_scale = fixed_scale,
						model = model,
						old = old,
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
						prior_delta_sd = prior_delta_sd)

	bayes_options <- list(n_iterations = n_iterations,
						n_chains = n_chains,
						n_cores = n_cores,
						keep_loglik = keep_loglik,
						random_parameters = random_parameters,
						seed = seed,
						max_tree_depth = max_tree_depth,
						adapt_delta = adapt_delta,
						show_stan_warnings = show_stan_warnings,
						lkj_shape_prior = lkj_shape_prior)

	start.time <- proc.time()

	stan_data <- processMDCEVdata(data, psi_formula, lc_formula, mle_options)

	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)

	# If no user supplied weights, replace weights with vector of ones
	if (is.null(weights))
		weights <-  rep(1, stan_data$I)

	stan_data$weights <- as.vector(weights)

	if (algorithm == "Bayes") {
		result <- BayesMDCEV(stan_data, bayes_options,
										initial.parameters,
										keep.samples = FALSE,
										include.stanfit = TRUE)

		# Get parameter estimates in matrix form
		result$est_pars <- extract(result$stan_fit, permuted = TRUE, inc_warmup = FALSE) %>%
				as.data.frame() %>%
				select(-starts_with("log_like"), -starts_with("sum_log_lik"),
					   -starts_with("tau_unif"), -.data$lp__)


	} else if (algorithm == "MLE") {
		result <- maxlikeMDCEV(stan_data, initial.parameters, mle_options)

		# Get parameter estimates in matrix form
		result$est_pars <- tbl_df(result[["stan_fit"]][["theta_tilde"]]) %>%
			select(-starts_with("log_like"), -starts_with("sum_log_lik"))

		result[["stan_fit"]][["theta_tilde"]] <- NULL

	}
	end.time <- proc.time()
	result$parms_info <- parms_info

	if(algorithm == "Bayes" || std_errors == "deltamethod")
		result$n_draws <- NULL

	result$stan_data <- stan_data
	result$algorithm <- algorithm
	result$random_parameters <- random_parameters
	result$psi_formula <- psi_formula
	result$lc_formula <- lc_formula
	result$n_classes <- n_classes
	result$model <- model
	result$algorithm <- algorithm
	result$std_errors <- std_errors
	result$n_draws <- n_draws
	result$n_individuals <- stan_data$I
	result$start.time <- start.time
	result$time.taken <- (end.time - start.time)[3]
	result$effective.sample.size <- sum(stan_data$weights)

	# Rename variables
	if (random_parameters == "fixed"){
		names(result$est_pars)[1:parms_info$n_vars$n_parms_total] <- parms_info$parm_names$all_names
	result$aic <- -2 * result$log.likelihood + 2 * parms_info$n_vars$n_parms_total
	result$bic <- -2 * result$log.likelihood + log(result$effective.sample.size) * parms_info$n_vars$n_parms_total
	}

	result$est_pars <- result$est_pars %>%
		tibble::rowid_to_column("sim_id") %>%
		tidyr::gather(parms, value, -sim_id)

#	class(result) <- "mdcev"

return(result)
}

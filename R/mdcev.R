#' @title mdcev
#' @description Fit a MDCEV model using MLE or Bayes
#' @param formula Formula for the model to be estimated.The formula is divided in
#' two parts, separated by the symbol \code{|}. The first part is reserved for
#' variables in the psi parameter. These can include alternative-specific and
#' individual-specific variables. The second part corresponds for individual-specific
#' variables that enter in the probability assignment in models with latent classes.
#' @param data The (IxJ) data to be passed to Stan including 1) id, 2) alt, 3) quant,
#' 4) price, 5) income, and columns for psi variables. Arrange data by id then alt.
#' Notes I is number of individuals and J is number of non-numeraire alternatives.
#' @param subset an optional vector specifying a subset of observations.
#' @param weights an optional vector of weights. Default to 1.
#' @param na.action a function wich indicated what should happen when the data
#' contains \code{NA}'s.
#' @param model A string indicating which model specification is estimated.
#' The options are "alpha", "gamma", "hybrid" and "hybrid0".
#' @param n_classes The number of latent classes.
#' @param fixed_scale1 Whether to fix scale at 1.
#' @param trunc_data Whether the estimation should be adjusted for truncation
#' @param seed Random seed.
#' @param algorithm Either "Bayes" for Bayes or "MLE"
#'     for maximum liklihood estimation.
#' @param flat_priors indicator if completely uninformative priors should be specified. If using MLE, the
#' optimizing function will then be equal to log-likelihood. Defaults to 1 if MLE used and 0 if Bayes used.
#' @param max_iterations Maximum number of iterations in MLE estimation.
#' @param print_iterations Whether to print iteration information
#' @param std_errors Compute standard errors using the delta method ("deltamethod")
#' or multivariate normal draws ("mvn"). The default is "mvn" as only mvn parameter draws are required
#' for demand and welfare simulation.
#' @param n_draws The number of multivariate normal draws for standard error calculations.
#' @param keep_loglik Whether to keep the log_lik calculations
#' @param hessian Wheter to keep the Hessian matrix
#' @param initial.parameters Specify initial parameters intead of starting at random.
#' Initial parameter values should be included in a named list. For the "hybrid" specification,
#' initial parameters can be specified as:
#' init = list(psi = array(0, dim = c(1, num_psi)),
#'             gamma = array(1, dim = c(1, num_alt)),
#'             alpha = array(0.5, dim = c(1, 0)),
#'             scale = array(1, dim = c(1)))
#' where num_psi is number of psi parameters and num_alt is number of non-numeraire alternatives
#' @param prior_psi_sd standard deviation for normal prior with mean 0.
#' @param prior_gamma_sd standard deviation for normal prior with mean 0.
#' @param prior_alpha_sd standard deviation for normal prior with mean 0.5.
#' @param prior_scale_sd standard deviation for normal prior with mean 1.
#' @param prior_delta_sd standard deviation for normal prior with mean 0.
#' @param alpha_fixed indicator if alpha parameters should be fixed (i.e. not random).
#' @param gamma_fixed indicator if gamma parameters should be fixed (i.e. not random).
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
#' @examples
#' \donttest{
#' data(data_rec, package = "rmdcev")
#'
#' data_rec <- mdcev.data(data_rec, subset = id < 500,
#'                 alt.var = "alt", choice = "quant")
#'
#' mdcev_est <- mdcev( ~ 1,
#' data = data_rec,
#' model = "hybrid0",
#' algorithm = "MLE")
#'}
mdcev <- function(formula = NULL, data, subset, na.action,
					 weights = NULL,
					 model = c("alpha", "gamma", "hybrid", "hybrid0"),
					 n_classes = 1,
					 fixed_scale1 = 0,
					 trunc_data = 0,
					 seed = "123",
					 max_iterations = 2000,
					 initial.parameters = NULL,
					 algorithm = c("MLE", "Bayes"),
					 flat_priors = NULL,
					 print_iterations = TRUE,
					 #mle_tol = 0.0001,
					 hessian = TRUE,
					 prior_psi_sd = 10,
					 prior_gamma_sd = 10,
					 prior_alpha_sd = 0.5,
					 prior_scale_sd = 1,
					 prior_delta_sd = 10,
					 gamma_fixed = 0,
					 alpha_fixed = 0,
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
					 lkj_shape_prior = 4,
				     ...)
{

	start.time <- proc.time()

	# Check models
	if (algorithm == "Bayes" && n_classes > 1)
		stop("Bayesian estimation can only be used with one class. Switch algorithm to MLE or choose n_classes = 1", "\n")

	if (algorithm == "MLE" && random_parameters != "fixed")
		stop("MLE can only be used with fixed parameters. Switch random_parameters to 'fixed' or change algorithm to Bayes", "\n")

	if (!inherits(data, "mdcev.data")) stop("Data must be of class mdcev.data")

	if (algorithm == "MLE" && is.null(flat_priors)){
		flat_priors <- 1
	} else if (algorithm == "Bayes" && is.null(flat_priors))
		flat_priors <- 0

	if (random_parameters == "fixed"){
		gamma_fixed <- 1
		alpha_fixed <- 1
	}

	# Put model options in a list
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

	stan_data <- processMDCEVdata(formula, data, mle_options)

	parms_info <- CreateParmInfo(stan_data, algorithm, random_parameters)

	# If no user supplied weights, replace weights with vector of ones
	if (is.null(weights))
		weights <-  rep(1, stan_data$I)

	stan_data$weights <- as.vector(weights)

	if (algorithm == "Bayes") {
		result <- BayesMDCEV(stan_data,
							 bayes_options,
							 initial.parameters,
							 keep.samples = FALSE,
							 include.stanfit = TRUE,
				 			 ...)

	} else if (algorithm == "MLE") {
		result <- maxlikeMDCEV(stan_data,
							   initial.parameters,
							   mle_options,
							   ...)

#		result[["stan_fit"]][["theta_tilde"]] <- NULL

#		if(length(names(result$est_pars)) == 0)
#			stop("Hessian matrix is not positive definite")

	}
	end.time <- proc.time()

	if(algorithm == "Bayes" || std_errors == "deltamethod")
		result$n_draws <- NULL

	result$parms_info <- parms_info
	result$stan_data <- stan_data
	result$algorithm <- algorithm
	result$random_parameters <- random_parameters
	result$formula <- formula
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
		result$aic <- -2 * result$log.likelihood + 2 * parms_info$n_vars$n_parms_total
		result$bic <- -2 * result$log.likelihood + log(result$effective.sample.size) * parms_info$n_vars$n_parms_total
	}



result <- structure(result,
					class = "mdcev")

return(result)
}

#' @title mdcev
#' @description Fit a MDCEV model using MLE or Bayes
#' @param formula Formula for the model to be estimated. The formula is divided in
#' three parts, separated by the symbol \code{|}. The first part is reserved for
#' alternative-specific and individual-specific variables in the psi parameters.
#' Note that alternative-specific constants are handled by the \code{psi_ascs} argument.
#' The second part corresponds for individual-specific variables that enter in the probability
#' assignment in models with latent classes. The third part is reserved for the $q_k$ variables
#' included in the $phi_k$ parameters in the KT model specification used in environmental economics
#' \code{model = "kt_ee"}.
#' @param data The (IxJ) data to be passed to Stan of class \code{\link[rmdcev]{mdcev.data}}
#'  including 1) id, 2) alt, 3) choice, 4) price, 5) income, and columns for alternative-specific and
#'  individual specific variables. Note: I is number of individuals and J is number of non-numeraire alternatives.
#' @param weights an optional vector of weights. Default to 1.
#' @param model A string indicating which model specification is estimated.
#' The options are "alpha", "gamma", "hybrid", "hybrid0", and "gamma1" for the MDCEV model and "kt_ee" for the
#' environmental economics Kuhn-Tucker specification. The "gamma1" model fixes alpha_0 = 1 (no income effects).
#' @param n_classes The number of latent classes. Note that the LC model is automatically estimated as long as the
#' prespecified number of classes is set greater than 1.
#' @param fixed_scale1 Whether to fix scale at 1.
#' @param single_scale For lc models, whether to estimate a single scale parameter
#' @param trunc_data Whether the estimation should be adjusted for truncation of non-numeraire alternatives.
#' This option is useful if the data only includes individuals with positive non-numeraire consumption levels
#' such as recreation data collected on-site. To account for the truncation of consumption, the likelihood is
#' normalized by one minus the likelihood of observing zero consumption for all non-numeraire alternatives
#' (i.e. likelihood of positive consumption) following Englin, Boxall and Watson (1998) and von Haefen (2003).
#' @param gamma_ascs Indicator to include alternative-specific gammas parameters.
#' @param psi_ascs Whether to include alternative-specific psi parameters. The first alternative is used as
#' the reference category. Only specify to 1 for MDCEV models.
#' @param seed Random seed.
#' @param algorithm Either "Bayes" for Bayes or "MLE" for maximum likelihood estimation.
#' @param max_iterations Maximum number of iterations in MLE.
#' @param jacobian_analytical_grad indicator whether to use analytical gradient method for Jacobian (=1) or numerical
#' gradient method (=0). For "kt_ee" model only,
#' @param print_iterations Whether to print iteration information
#' @param std_errors Compute standard errors using the delta method ("deltamethod")
#' or multivariate normal draws ("mvn"). The default is "deltamethod". Note that mvn parameter draws should be
#' used to incorporate parameter uncertainty for demand and welfare simulation. For maximum likelihood estimation only.
#' @param n_draws The number of multivariate normal draws for standard error calculations if "mvn" is specified.
#' @param keep_loglik Whether to keep the log_lik calculations
#' @param hessian Whether to keep the Hessian matrix
#' @param initial.parameters The default for fixed and random parameter specifications is to use random starting values.
#' (except for the scale parameter with a starting value set to 1).
#' For LC models, the default is to use slightly adjusted MLE point estimates from the single class model.
#' Initial parameter values should be included in a named list. For example, the LC "hybrid" specification
#' initial parameters can be specified as:
#' initial.parameters = list(psi = array(0, dim = c(K, num_psi)),
#'                             gamma = array(1, dim = c(K, num_alt)),
#'                             alpha = array(0.5, dim = c(K, 0)),
#'                             scale = array(1, dim = c(K)))
#' where K is the number of classes (i.e. K = 1 is used for single class models),
#' num_psi is number of psi parameters, and num_alt is number of non-numeraire alternatives.
#' @param flat_priors indicator if completely uninformative priors should be specified. Defaults to 1 if MLE used and 0 if Bayes used. If using MLE and set flat_priors = 0,
#' penalized MLE is used and the optimizing objective is augmented with the priors.
#' @param prior_psi_sd standard deviation for normal prior with mean 0.
#' @param prior_phi_sd standard deviation for normal prior with mean 0.
#' @param prior_gamma_sd standard deviation for half-normal prior with mean 1.
#' @param prior_alpha_shape shape parameter for beta distribution.
#' @param prior_scale_sd standard deviation for half-normal prior with mean 0.
#' @param prior_delta_sd standard deviation for normal prior with mean 0.
#' @param alpha_nonrandom Logical. If \code{TRUE}, alpha parameters are not random (no individual-specific standard deviation).
#' @param gamma_nonrandom Logical. If \code{TRUE}, gamma parameters are not random (no individual-specific standard deviation).
#' @param psi_random A one-sided formula specifying which psi formula terms should be treated as random
#'   coefficients (i.e. have individual-specific draws with estimated population mean and standard deviation).
#'   Only applicable when \code{algorithm = "Bayes"} and \code{random_parameters} is \code{"uncorr"} or
#'   \code{"corr"}. When \code{NULL} (default), all formula terms are random (existing behaviour).
#'   Alternative-specific constants (from \code{psi_ascs = 1}) are always random regardless of this argument.
#'   Formula terms not listed in \code{psi_random} become fixed (pooled) point estimates.
#'   Example: if the psi formula is \code{~ quality + income} and you specify
#'   \code{psi_random = ~ quality}, then \code{quality} gets an individual-specific coefficient
#'   while \code{income} is estimated as a single pooled parameter.
#' @param lkj_shape_prior Prior for Cholesky matrix
#' @param n_iterations The number of iterations to use in Bayesian estimation. The default is for the number of
#' iterations to be split evenly between warmup and posterior draws. The number of warmup draws can be directly controlled using the warmup argument (see \code{rstan::sampling}).
#' @param n_chains The number of independent Markov chains in Bayesian estimation.
#' @param n_cores The number of cores used to execute the Markov chains in parallel in Bayesian estimation.
#' Can set using options(mc.cores = parallel::detectCores()).
#' @param random_parameters The form of the covariance matrix for
#'     Bayes. Can be 'fixed' for no random parameters, 'uncorr' for uncorrelated random parameters, or
#'     'corr' for correlated random parameters.
#' @param max_tree_depth
#'     https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#' @param adapt_delta
#'    https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#' @param show_stan_warnings Whether to show warnings from Stan.
#' @param backend Estimation backend. Either \code{"cmdstanr"} (default) or \code{"rstan"}.
#'   The cmdstanr backend reads \code{.stan} files from \code{inst/stan/} at runtime.
#'   The rstan backend uses pre-compiled C++ in \code{src/stanExports_*.h}.
#' @param ... Additional parameters to pass on to \code{rstan::stan}
#'     and \code{rstan::sampling}.
#' @return A object of class mdcev
#' @export
#' @examples
#' \donttest{
#' data(data_rec, package = "rmdcev")
#'
#' data_rec <- mdcev.data(data_rec, subset = id <= 500, id.var = "id",
#'                 alt.var = "alt", choice = "quant")
#'
#' mdcev_est <- mdcev( ~ 0,
#' data = data_rec,
#' model = "hybrid0",
#' algorithm = "MLE",
#' backend = "rstan")
#'}
mdcev <- function(formula = NULL, data,
				 weights = NULL,
				 model = c("alpha", "gamma", "hybrid", "hybrid0", "kt_ee", "gamma1"),
				 n_classes = 1,
				 fixed_scale1 = FALSE,
				 single_scale = FALSE,
				 trunc_data = FALSE,
				 psi_ascs = NULL,
				 gamma_ascs = TRUE,
				 seed = 123L,
				 max_iterations = 2000,
				 jacobian_analytical_grad = TRUE,
				 initial.parameters = "random",
				 hessian = TRUE,
				 algorithm = c("MLE", "Bayes"),
				 flat_priors = NULL,
				 print_iterations = TRUE,
				 prior_psi_sd = 10,
				 prior_gamma_sd = 10,
				 prior_phi_sd = 10,
				 prior_alpha_shape = 1,
				 prior_scale_sd = 1,
				 prior_delta_sd = 10,
				 gamma_nonrandom = FALSE,
				 alpha_nonrandom = FALSE,
				 psi_random = NULL,
				 std_errors = "deltamethod",
				 n_draws = 50,
				 keep_loglik = FALSE,
				 random_parameters = "fixed",
				 show_stan_warnings = TRUE,
				 n_iterations = 200,
				 n_chains = 4,
				 n_cores = 4,
				 max_tree_depth = 10,
				 adapt_delta = 0.8,
				 lkj_shape_prior = 4,
				 backend = "cmdstanr",
			     ...)
{

	start.time <- proc.time()

	# Check models
	if (!is.element(algorithm, c("MLE", "Bayes"))) stop("algorithm must be 'MLE' or 'Bayes'")
	if (!is.element(random_parameters, c("fixed", "uncorr", "corr"))) stop("random_parameters must be 'fixed', 'uncorr' or 'corr'")
	if (!is.element(model, c("alpha", "gamma", "hybrid", "hybrid0", "kt_ee", "gamma1"))) stop("model must be 'alpha', 'gamma', 'hybrid', 'hybrid0', 'kt_ee', or 'gamma1'")
	if (!is.element(std_errors, c("deltamethod", "mvn"))) stop("std_errors must be 'deltamethod' or 'mvn'")

	if (algorithm == "Bayes" && n_classes > 1)
		stop("Bayesian estimation can only be used with one class. Switch algorithm to MLE or choose n_classes = 1", "\n")

	if (algorithm == "MLE" && random_parameters != "fixed")
		stop("MLE can only be used with fixed parameters. Switch random_parameters to 'fixed' or change algorithm to Bayes", "\n")

	if (!inherits(data, "mdcev.data")) stop("Data must be of class mdcev.data")

	if (isTRUE(fixed_scale1) && isTRUE(single_scale)) stop("Cannot set both fixed_scale1 and single_scale to TRUE")

	if (algorithm == "MLE" && is.null(flat_priors)){
		flat_priors <- TRUE
	} else if (algorithm == "Bayes" && is.null(flat_priors))
		flat_priors <- FALSE

	if (random_parameters == "fixed"){
		gamma_nonrandom <- TRUE
		alpha_nonrandom <- TRUE
	}

	if(algorithm == "Bayes" || std_errors == "deltamethod")
		n_draws <- 0

	# Put model options in a list
	mle_options <- list(initial.parameters = initial.parameters,
						seed = seed,
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
						prior_phi_sd = prior_phi_sd,
						prior_gamma_sd = prior_gamma_sd,
						prior_alpha_shape = prior_alpha_shape,
						prior_scale_sd = prior_scale_sd,
						prior_delta_sd = prior_delta_sd,
						gamma_nonrandom = gamma_nonrandom,
						alpha_nonrandom = alpha_nonrandom)

	bayes_options <- list(initial.parameters = initial.parameters,
						  n_iterations = n_iterations,
						n_chains = n_chains,
						n_cores = n_cores,
						keep_loglik = keep_loglik,
						random_parameters = random_parameters,
						seed = seed,
						max_tree_depth = max_tree_depth,
						adapt_delta = adapt_delta,
						show_stan_warnings = show_stan_warnings,
						lkj_shape_prior = lkj_shape_prior)

	# Need for naming gamma/alpha parameters
	alt_names <- as.character(unlist(unique(attr(data, "index")["alt"])))

	stan_data <- processMDCEVdata(formula, data, stan_model_options)

	# Always initialise fixed-psi fields; override via split when psi_random is set.
	stan_data$NPsi_ij_fixed <- 0L
	stan_data$dat_psi_fixed  <- matrix(0, 0, 0)

	if (algorithm == "Bayes" && random_parameters != "fixed" && !is.null(psi_random)) {
		if (!inherits(psi_random, "formula"))
			stop("psi_random must be a one-sided formula, e.g. ~ quality + income", call. = FALSE)
		stan_data <- .split_psi_for_rp(stan_data, psi_random, data)
	}

	parms_info <- CreateParmInfo(stan_data, alt_names, algorithm, random_parameters)

		# If no user supplied weights, replace weights with vector of ones
	if (is.null(weights))
		weights <-  rep(1, stan_data$I)

	stan_data$weights <- as.vector(weights)

	if (algorithm == "Bayes") {
		result <- BayesMDCEV(stan_data,
							 bayes_options,
							 keep.samples = FALSE,
							 include.stanfit = TRUE,
							 backend = backend,
				 			 ...)

	} else if (algorithm == "MLE") {
		result <- maxlikeMDCEV(stan_data,
							   mle_options,
							   parms_info,
							   backend = backend,
							   ...)

	}
	end.time <- proc.time()

	result$parms_info <- parms_info
	result$stan_data <- stan_data
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

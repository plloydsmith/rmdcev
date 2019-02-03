#' @title FitMDCEV
#' @description Fit a MDCEV model using MLE or HB
#' @param data The data to be passed to Stan. Must include quant (IxJ), price (IxJ), income (Ix1), and dat_psi(IJxNPsi).
#' @param data_class The data for class membership to be passed to Stan.
#' @param weights An optional vector of sampling or frequency weights.
#' @param price_num An optional vector containing price of numeraire or outside good (default is 1).
#' @param model_specification A string indicating which model specification is estimated.
#' The options are "alpha","les", "gamma", and "gamma0".
#' @param n_classes The number of latent classes.
#' @param fixed_scale Whether to fix scale at 1.
#' @param trunc_data Whether the estimation should be adjusted for truncation
#' @param seed Random seed.
#' @param algorithm Either "HB-Stan" for Hierarchical Bayes or "MLE"
#'     for maximum liklihood estimation.
#' @param print_ll Whether to print logliklihood at each iteration
#' @param n_draws The number of MVN draws for standard error calculations
#' @param keep_loglik Whether to keep the log_lik calculations
#' @param hessian Wheter to keep the Hessian matrix
#' @param initial.parameters Specify initial parameters intead of
#'     starting at random.
# #' @param n_iterations The number of iterations in Hierarchical Bayes.
# #' @param n_chains The number of chains in Hierarchical Bayes.
# #' @param hb_random_parameters The form of the covariance matrix for
# #'     Hierarchical Bayes. Can be 'uncorr, 'corr'.
# #' @param hb.iterations The number of iterations in Hierarchical
# #'     Bayes.
# #' @param hb.chains The number of chains in Hierarchical Bayes.
# #' @param hb.max.tree.depth
# #'     http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
# #' @param hb.adapt.delta
# #'     http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# #' @param hb.keep.samples Whether to keep the samples of all the
# #'     parameters in the output.
# #' @param hb.stanfit Whether to include the stanfit property.
# #' @param hb.lkj.prior.shape Real number greater than one; the shape
# #'     hyperparameter for the LKJ prior used for the correlation matrix
# #'     of the respondent coefficients distribution. A value of one gives
# #' equal probability weight to all possible correlation matrices. Larger values
# #' favour less correlation (draws closer to the identity matrix).
# #' @param hb.warnings Whether to show warnings from Stan.
# #' @param hb.beta.draws.to.keep Maximum number of beta draws per
# #'     respondent to return in beta.draws.
# #' @param include.choice.parameters Whether to include
# #'     alternative-specific parameters.
#' @param ... Additional parameters to pass on to \code{rstan::stan}
#'     and \code{rstan::sampling}.
#' @return A stanfit object
#' @export
#'
FitMDCEV <- function(data,
					 data_class = NULL,
					 weights = NULL,
					 price_num = NULL,
					 model_specification = NULL,
					 n_classes = 1,
					 fixed_scale = 0,
					 trunc_data = trunc_data,
					 seed = 123,
					 initial.parameters = NULL,
					 algorithm = "MLE", # HB-Stan", # MLE
					 #	std_errors = "draws", # still need to implement
					 print_ll = 0,
					 #mle_tol = 0.0001,
					 hessian = TRUE,
					 n_draws = 50,
					 keep_loglik = 0)#,
					 #subset = NULL,
#					 hb_random_parameters = hb_random_parameters,
#					 n_iterations = 100, n_chains = 4,
#					 hb.max.tree.depth = 10, hb.adapt.delta = 0.8,
#					 hb.keep.samples = FALSE, hb.stanfit = TRUE,
					 #					 hb.prior.mean = 0, hb.prior.sd = 5,
					 #					 hb.sigma.prior.shape = 1.39435729464721,
					 #					 hb.sigma.prior.scale = 0.39435729464721,
#					 hb.lkj.prior.shape = 4,
#					 hb.warnings = TRUE, hb.beta.draws.to.keep = 0)
#simulated.priors = NULL,
#simulated.priors.from.design = FALSE,
#simulated.sample.size = 300,
#synthetic.priors = NULL,
#synthetic.priors.from.design = NULL,
#synthetic.sample.size = NULL,
#tasks.left.out = 0,
#lc.tolerance = 0.0001,

#include.choice.parameters = TRUE,
#respondent.ids = NULL, ...

{
	if (algorithm == "HB-Stan")
		stop("Not set up for HB yet.")

		if (algorithm == "HB-Stan" && !is.null(weights))
		stop("Weights are not able to be applied for Hierarchical Bayes.")

	if (algorithm == "HB-Stan" && n_classes > 1)
		stop("Hierarchical Bayes can only be used with one class. Switch algorithm to MLE")

	if (identical(dim(data$price), dim(data$quant)) == FALSE)
		stop("Price and quant dimension mismatch. Ensure dim(price) = dim(quant)")

	model_options <- list(fixed_scale = fixed_scale,
						  model_specification = model_specification,
						  n_classes = n_classes,
						  fixed_scale = fixed_scale,
						  trunc_data = trunc_data,
						  print_ll = print_ll,
						  hessian = TRUE,
						  n_draws = 50,
						  keep_loglik = 0)

#	hb_options <- list(n_iterations = n_iterations, n_chains = n_chains,
#						  hb.max.tree.depth = 10, hb.adapt.delta = 0.8,
#						  hb.keep.samples = FALSE, hb.stanfit = TRUE,
						  #					 hb.prior.mean = 0, hb.prior.sd = 5,
						  #					 hb.sigma.prior.shape = 1.39435729464721,
						  #					 hb.sigma.prior.scale = 0.39435729464721,
#						  hb.lkj.prior.shape = 4,
#						  hb.warnings = TRUE, hb.beta.draws.to.keep = 0)

	start.time <- proc.time()

#	dat <- stan.dat

	stan_data <- processMDCEVdata(data, data_class, price_num, model_options)

	# Replace weights with vector of one's if missing
	if (is.null(weights))
		weights <-  rep(1, stan_data$I)

	stan_data$weights <- as.vector(weights)
	stan_data$print_ll <- print_ll

#	if (algorithm == "HB")
#	{
#		result <- HierarchicalBayesMDCEV(stan_data, model_options, hb_options,
#										initial.parameters = initial.parameters,
#										seed = seed,
#										keep.samples = FALSE,
#										include.stanfit = TRUE,
#										hb_random_parameters = hb_random_parameters,
#									  show.stan.warnings = TRUE)
#	} else if (algorithm == "MLE")
#	{

	result <- maxlikeMDCEV(stan_data, initial.parameters, seed, model_options)
						    #lc.tolerance
#	}

	end.time <- proc.time()

	result$stan_data <- stan_data
	result$algorithm <- algorithm
	result$n_classes <- n_classes
	result$weights <- weights
	#	result$weights.description <- if (is.null(weights)) NULL else Labels(weights)
	result$n_respondents <- dat$I
	result$n_alternatives <- dat$n_alternatives
	result$time.taken <- (end.time - start.time)[3]
	class(result) <- "EstimateMDCEV"
	result
}

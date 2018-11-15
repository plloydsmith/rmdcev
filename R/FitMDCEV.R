#' @title FitMDCEV
#' @description Fit a MDCEV model using MLE or HB
#' @param stan.dat The data to be passed to Stan.
#' @param n_classes The number of latent classes.
#' @param weights An optional vector of sampling or frequency weights.
#' @param seed Random seed.
#' @param algorithm Either "HB" for Hierarchical Bayes or "MLE"
#'     for latent class analysis.
#' @param initial.parameters Specify initial parameters intead of
#'     starting at random in latent class analysis. The initial
#'     parameters need to be supplied as list consisting of a matrix
#'     called class.parameters whose columns are the parameters of the
#'     classes, and a vector called class.sizes containing the class
#'     size parameters.
#' @param normal.covariance The form of the covariance matrix for
#'     Hierarchical Bayes. Can be 'Full, 'Spherical', 'Diagonal'.
#' @param hb.iterations The number of iterations in Hierarchical
#'     Bayes.
#' @param hb.chains The number of chains in Hierarchical Bayes.
#' @param hb.max.tree.depth
#'     http://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#' @param hb.adapt.delta
#'     http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#' @param hb.keep.samples Whether to keep the samples of all the
#'     parameters in the output.
#' @param hb.stanfit Whether to include the stanfit property.
#' @param hb.lkj.prior.shape Real number greater than one; the shape
#'     hyperparameter for the LKJ prior used for the correlation matrix
#'     of the respondent coefficients distribution. A value of one gives
#' equal probability weight to all possible correlation matrices. Larger values
#' favour less correlation (draws closer to the identity matrix).
#' @param hb.warnings Whether to show warnings from Stan.
#' @param hb.beta.draws.to.keep Maximum number of beta draws per
#'     respondent to return in beta.draws.
#' @param include.choice.parameters Whether to include
#'     alternative-specific parameters.
#' @param ... Additional parameters to pass on to \code{rstan::stan}
#'     and \code{rstan::sampling}.
#' @export
#'
FitMDCEV <- function(dat,
					 dat_class = NULL,
					 weights = NULL,
					 price_num = NULL,
					 model_specification = NULL,
					 n_classes = 1,
					 fixed_scale = 0,
					 trunc_data = trunc_data,
					 seed = 123,
					 initial.parameters = NULL,
					 algorithm = "HB-Stan", # MLE
					 #	std_errors = "draws", # still need to implement
					 print_ll = 0,
					 #mle_tol = 0.0001,
					 hessian = TRUE,
					 n_draws = 50,
					 keep_loglik = 0,
					 #subset = NULL,
					 hb_random_parameters = hb_random_parameters,
					 n_iterations = 100, n_chains = 4,
					 hb.max.tree.depth = 10, hb.adapt.delta = 0.8,
					 hb.keep.samples = FALSE, hb.stanfit = TRUE,
					 #					 hb.prior.mean = 0, hb.prior.sd = 5,
					 #					 hb.sigma.prior.shape = 1.39435729464721,
					 #					 hb.sigma.prior.scale = 0.39435729464721,
					 hb.lkj.prior.shape = 4,
					 hb.warnings = TRUE, hb.beta.draws.to.keep = 0)
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
	if (algorithm == "HB-Stan" && !is.null(weights))
		stop("Weights are not able to be applied for Hierarchical Bayes.")

	if (algorithm == "HB-Stan" && n_classes > 1)
		stop("Hierarchical Bayes can only be used with one class. Switch algorithm to MLE")

	if (identical(dim(price), dim(quant)) == FALSE)
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

	hb_options <- list(n_iterations = n_iterations, n_chains = n_chains,
						  hb.max.tree.depth = 10, hb.adapt.delta = 0.8,
						  hb.keep.samples = FALSE, hb.stanfit = TRUE,
						  #					 hb.prior.mean = 0, hb.prior.sd = 5,
						  #					 hb.sigma.prior.shape = 1.39435729464721,
						  #					 hb.sigma.prior.scale = 0.39435729464721,
						  hb.lkj.prior.shape = 4,
						  hb.warnings = TRUE, hb.beta.draws.to.keep = 0)

	start.time <- proc.time()

#	dat <- stan.dat

	stan_data <- processMDCEVdata(dat, dat_class, price_num, model_options)

	# Replace weights with vector of one's if missing
	if (is.null(weights))
		weights <-  rep(1, stan_data$I)

	stan_data$weights <- as.vector(weights)
	stan_data$print_ll <- print_ll

	if (algorithm == "HB")
	{
		result <- HierarchicalBayesMDCEV(stan_data, model_options, hb_options,
										initial.parameters = initial.parameters,
										seed = seed,
										keep.samples = FALSE,
										include.stanfit = TRUE,
										hb_random_parameters = hb_random_parameters,
									  show.stan.warnings = TRUE)
	} else if (algorithm == "MLE")
	{

	result <- maxlikeMDCEV(stan_data, seed, model_options,
						   initial.parameters) #lc.tolerance
	}

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

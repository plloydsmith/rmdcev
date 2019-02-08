#' @title FitMDCEV
#' @description Fit a MDCEV model using MLE or HB
#' @param data The data to be passed to Stan including
#' 1) quant (IxJ) where each row is individual and each column is good,
#' 2) price (IxJ) where each row is individual and each column is good,
#' 3) income (Ix1) where each row is individual, and
#' 4) dat_psi(IJxNPsi) where each row is a individual-good combination.
#' Notes I is number of individuals and J is number of non-numeraire goods.
#' @param data_class The data for class membership to be passed to Stan.
#' @param weights An optional vector of sampling or frequency weights.
#' @param num_price An optional vector containing price of numeraire or outside good (default is 1).
#' @param model A string indicating which model specification is estimated.
#' The options are "alpha","les", "gamma", and "gamma0".
#' @param n_classes The number of latent classes.
#' @param fixed_scale Whether to fix scale at 1.
#' @param trunc_data Whether the estimation should be adjusted for truncation
#' @param seed Random seed.
#' @param algorithm Either "HB" for Hierarchical Bayes or "MLE"
#'     for maximum liklihood estimation.
#' @param print_ll Whether to print logliklihood at each iteration
#' @param n_draws The number of MVN draws for standard error calculations
#' @param keep_loglik Whether to keep the log_lik calculations
#' @param hessian Wheter to keep the Hessian matrix
#' @param initial.parameters Specify initial parameters intead of
#'     starting at random.
#' @param n_iterations The number of iterations in Hierarchical Bayes.
#' @param n_chains The number of chains in Hierarchical Bayes.
#' @param hb_random_parameters The form of the covariance matrix for
# #'     Hierarchical Bayes. Can be 'uncorr, 'corr'.
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
					 num_price = NULL,
					 model = c("alpha", "les", "gamma", "gamma0"),
					 n_classes = 1,
					 fixed_scale = 0,
					 trunc_data = trunc_data,
					 seed = 123,
					 initial.parameters = NULL,
					 algorithm = c("MLE", "HB"),
					 #	std_errors = "draws", # still need to implement
					 print_ll = 0,
					 #mle_tol = 0.0001,
					 hessian = TRUE,
					 n_draws = 30,
					 keep_loglik = 0,
					 #subset = NULL,
					 hb_random_parameters = "fixed",
					 n_iterations = 200, n_chains = 4)#,
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
#data <-stan.dat
#include.choice.parameters = TRUE,
#respondent.ids = NULL, ...

{
#	if (algorithm == "HB")
#		stop("Not set up for HB yet.")
	CheckMdcevData(data)

	if (algorithm == "HB" && !is.null(weights))
		stop("Weights are not able to be applied for Hierarchical Bayes.")

	if (algorithm == "HB" && n_classes > 1)
		stop("Hierarchical Bayes can only be used with one class. Switch algorithm to MLE")

	mle_options <- list(fixed_scale = fixed_scale,
						  model = model,
						  n_classes = n_classes,
						  trunc_data = trunc_data,
						  print_ll = print_ll,
						  hessian = hessian,
						  n_draws = n_draws,
						  keep_loglik = keep_loglik)

	hb_options <- list(fixed_scale = fixed_scale,
					   model = model,
					   n_classes = n_classes,
					   trunc_data = trunc_data,
					   n_iterations = n_iterations,
					   n_chains = n_chains,
					   keep_loglik = keep_loglik,
					   hb_random_parameters = hb_random_parameters,
					   seed = seed,
						  hb.max.tree.depth = 10, hb.adapt.delta = 0.8,
						  hb.keep.samples = FALSE, hb.stanfit = TRUE,
						  hb.prior.mean = 0, hb.prior.sd = 5,
						  hb.sigma.prior.shape = 1.39435729464721,
						  hb.sigma.prior.scale = 0.39435729464721,
						  hb.lkj.prior.shape = 4,
						  hb.warnings = TRUE, hb.beta.draws.to.keep = 0)

	start.time <- proc.time()

#	data <- stan.dat

	stan_data <- processMDCEVdata(data, data_class, num_price, mle_options)

	# If no user supplied weights, replace weights with vector of ones
	if (is.null(weights))
		weights <-  rep(1, stan_data$I)

	stan_data$weights <- as.vector(weights)
	stan_data$print_ll <- print_ll

	if (algorithm == "HB") {
		result <- HierarchicalBayesMDCEV(stan_data, hb_options,
										initial.parameters,
										keep.samples = FALSE,
										include.stanfit = TRUE,
									  show.stan.warnings = TRUE)

	###########################################################################
	# Get parameter estimates in matrix form
	result$est_pars <- extract(result$stan_fit, permuted = TRUE, inc_warmup = FALSE) %>%
			as.data.frame() %>%
			select(-starts_with("log_like"), -starts_with("sum_log_lik"), -.data$lp__)

	} else if (algorithm == "MLE") {

	result <- maxlikeMDCEV(stan_data, initial.parameters, seed, mle_options)

	###########################################################################
	# Get parameter estimates in matrix form
	result$est_pars <- tbl_df(result[["stan_fit"]][["theta_tilde"]]) %>%
		select(-starts_with("log_like"), -starts_with("sum_log_lik"))

	result[["stan_fit"]][["theta_tilde"]] <- NULL

	}
	end.time <- proc.time()

	result$est_pars <- result$est_pars %>%
		rowid_to_column("sim_id") %>%
		gather(parms, value, -sim_id)
#		gather_(key_col = 'parms',
#				value_col = 'value', -sim_id, factor_key=TRUE)

	result$stan_data <- stan_data
	result$algorithm <- algorithm
	result$n_classes <- n_classes
	result$n_draws <- n_draws
	#	result$weights.description <- if (is.null(weights)) NULL else Labels(weights)
	result$n_respondents <- data$I
	result$n_alternatives <- data$n_alternatives
	result$time.taken <- (end.time - start.time)[3]
	class(result) <- "mdcev"
	result
}

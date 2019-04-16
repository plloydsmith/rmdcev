#' @title FitMDCEV
#' @description Fit a MDCEV model using MLE or Bayes
#' @param psi_formula Formula for psi
#' @param lc_formula Formula for latent class
#' @param data The (IxJ) data to be passed to Stan including 1) id, 2) good, 3) quant,
#' 4) price, 5) income, and columns for psi variables. Arrange data by id then good.
#' Notes I is number of individuals and J is number of non-numeraire goods.
#' @param weights An optional vector of sampling or frequency weights.
#' @param num_price An optional vector containing price of numeraire or outside good (default is 1).
#' @param model A string indicating which model specification is estimated.
#' The options are "alpha","les", "gamma", and "gamma0".
#' @param n_classes The number of latent classes.
#' @param fixed_scale Whether to fix scale at 1.
#' @param trunc_data Whether the estimation should be adjusted for truncation
#' @param seed Random seed.
#' @param algorithm Either "Bayes" for Hierarchical Bayes or "MLE"
#'     for maximum liklihood estimation.
#' @param no_priors indicator if completely uninformative priors should be specified. If using MLE, the
#' optimizing function will then be equal to log-likelihood. Defaults to 1 if MLE used and 0 if Bayes used.
#' @param print_iterations Whether to print iteration information
#' @param n_draws The number of MVN draws for standard error calculations
#' @param keep_loglik Whether to keep the log_lik calculations
#' @param hessian Wheter to keep the Hessian matrix
#' @param initial.parameters Specify initial parameters intead of
#'     starting at random.
#' @param prior_psi_sd standard deviation for normal prior with mean 0.
#' @param prior_gamma_sd standard deviation for normal prior with mean 0.
#' @param prior_alpha_sd standard deviation for normal prior with mean 0.5.
#' @param prior_scale_sd standard deviation for normal prior with mean 1.
#' @param prior_beta_m_sd standard deviation for normal prior with mean 0.
#' @param n_iterations The number of iterations in Hierarchical Bayes.
#' @param n_chains The number of chains in Hierarchical Bayes.
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
					 no_priors = NULL,
					 print_iterations = TRUE,
					 #mle_tol = 0.0001,
					 hessian = TRUE,
					 prior_psi_sd = 10,
					 prior_gamma_sd = 10,
					 prior_alpha_sd = 0.5,
					 prior_scale_sd = 1,
					 prior_beta_m_sd = 10,
					 n_draws = 30,
					 keep_loglik = 0,
					 random_parameters = "fixed",
					 show_stan_warnings = TRUE,
					 n_iterations = 200,
					 n_chains = 4,
					 max_tree_depth = 10,
					 adapt_delta = 0.8)
{
	CheckMdcevData(data)

	if (algorithm == "Bayes" && !is.null(weights))
		stop("Weights are not able to be applied for Hierarchical Bayes.")

	if (algorithm == "Bayes" && n_classes > 1)
		stop("Hierarchical Bayes can only be used with one class. Switch algorithm to MLE or choose n_classes = 1", "\n")

	if (algorithm == "MLE" && is.null(no_priors)){
		no_priors = 1
	} else if (algorithm == "Bayes" && is.null(no_priors))
		no_priors = 0

	mle_options <- list(fixed_scale = fixed_scale,
						model = model,
						n_classes = n_classes,
						trunc_data = trunc_data,
						seed = seed,
						print_iterations = print_iterations,
						hessian = hessian,
						n_draws = n_draws,
						keep_loglik = keep_loglik,
						no_priors = no_priors,
						prior_psi_sd = prior_psi_sd,
						prior_gamma_sd = prior_gamma_sd,
						prior_alpha_sd = prior_alpha_sd,
						prior_scale_sd = prior_scale_sd,
						prior_beta_m_sd = prior_beta_m_sd)

	bayes_options <- list(n_iterations = n_iterations,
						n_chains = n_chains,
						keep_loglik = keep_loglik,
						random_parameters = random_parameters,
						seed = seed,
						max_tree_depth = max_tree_depth,
						adapt_delta = adapt_delta,
						show_stan_warnings = show_stan_warnings)

	start.time <- proc.time()

	stan_data <- processMDCEVdata(data, psi_formula, lc_formula, num_price, mle_options)

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
				select(-starts_with("log_like"), -starts_with("sum_log_lik"), -.data$lp__)

	} else if (algorithm == "MLE") {

		result <- maxlikeMDCEV(stan_data, initial.parameters, mle_options)

		# Get parameter estimates in matrix form
		result$est_pars <- tbl_df(result[["stan_fit"]][["theta_tilde"]]) %>%
			select(-starts_with("log_like"), -starts_with("sum_log_lik"))

		result[["stan_fit"]][["theta_tilde"]] <- NULL

	}
	end.time <- proc.time()
	# Rename psi variables
	if (n_classes == 1){
		psi.names <- paste0(rep('psi', ncol(stan_data[["dat_psi"]])), sep="_",
							colnames(stan_data[["dat_psi"]]))
		original.names <- colnames(result$est_pars)
		new.names <- c(psi.names, original.names[-c(1:length(psi.names))])
		colnames(result$est_pars) <- new.names
	}
	# LC names still to do
#	class.names <- colnames(stan_est[["stan_data"]][["data_class"]])
#	parms <- c(paste(rep('betam', length(class.names)), 1:length(class.names), sep=""))
#	names <- cbind(parms, paste(rep(class, class.names)
	result$est_pars <- result$est_pars %>%
		tibble::rowid_to_column("sim_id") %>%
		tidyr::gather(parms, value, -sim_id)

	stan_data$M_factorial <- NULL

	result$stan_data <- stan_data
	result$algorithm <- algorithm
	result$psi_formula <- psi_formula
	result$lc_formula <- lc_formula
	result$n_classes <- n_classes
	result$model <- model
	result$algorithm <- algorithm
	result$n_draws <- n_draws
	result$n_respondents <- stan_data$I
	result$start.time <- start.time
	result$time.taken <- (end.time - start.time)[3]

	if(algorithm == "Bayes")
		result$n_draws <- NULL
	class(result) <- "mdcev"

result
}

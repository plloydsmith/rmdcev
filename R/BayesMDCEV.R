#' @title BayesMDCEV
#' @description Fit a MDCEV model using Bayesian estimation and Stan
#' @param bayes_options list of Bayes options
#' @param stan_data data for model
#' @inheritParams FitMDCEV
#' @param keep.samples default is FALSE,
#' @param include.stanfit default isTRUE,
#' @import dplyr
#' @import rstan
#' @keywords internal
BayesMDCEV <- function(stan_data, bayes_options,
								 initial.parameters,
								 keep.samples = FALSE,
								 include.stanfit = TRUE) {

	if (bayes_options$n_iterations <= 0)
		stop("The specified number of iterations must be greater than 0.")

	# allows Stan chains to run in parallel on multiprocessor machines
	options(mc.cores = parallel::detectCores())

	# Create indices for individual level psi parameters
	indexes <- tibble(individual = rep(1:stan_data$I, each = stan_data$J),
						  task = rep(1:stan_data$I, each = stan_data$J),
						  row = 1:(stan_data$I*stan_data$J)) %>%
		group_by(task) %>%
		summarise(task_individual = first(individual),
				  start = first(row),
				  end = last(row))

	stan_data$start = indexes$start
	stan_data$end = indexes$end
	stan_data$task_individual = indexes$task_individual
	stan_data$task = indexes$task
	stan_data$IJ = stan_data$I * stan_data$J
	stan_data$lkj_shape = bayes_options$lkj_shape_prior

	stan_data$K <- 1
	stan_data$L <- 0
	stan_data$data_class <- matrix(0, stan_data$I, 0)

	if (bayes_options$random_parameters == "fixed"){
		message("Using Bayes to estimate MDCEV")
		stan.model <- stanmodels$mdcev
	}else if (bayes_options$random_parameters == "uncorr"){
		message("Using Bayes to estimate uncorrelated RP-MDCEV")
		stan.model <- stanmodels$mdcev_lc
		stan_data$corr <- 0
	}else if (bayes_options$random_parameters == "corr"){
		message("Using Bayes to estimate correlated RP-MDCEV")
		stan.model <- stanmodels$mdcev_lc
		stan_data$corr <- 1
	}


	if (bayes_options$show_stan_warnings == FALSE)
		suppressWarnings(stan_fit <- RunStanSampling(stan_data, stan.model, bayes_options))
	else if(bayes_options$show_stan_warnings == TRUE)
		stan_fit <- RunStanSampling(stan_data, stan.model, bayes_options)

	if(bayes_options$n_chains == 1)
		chain_index <- 1
	else if(bayes_options$n_chains > 1)
		chain_index <- bayes_options$n_chains+1

	result <- list()
	result$stan_fit <- stan_fit
	result$log.likelihood <- rstan::get_posterior_mean(result$stan_fit, pars = "sum_log_lik")[,chain_index]

#	result$stan_fit$par[["theta"]] <- NULL
#	result$stan_fit$par[["beta_m"]] <- NULL
	result
}

#' @title RunStanSampling
#' @description Wrapper function for \code{rstan:stan} and
#' \code{rstan:sampling} to run Stan Bayes analysis.
#' @inheritParams BayesMDCEV
#' @param stan.model Complied Stan model
#' @param ... Additional parameters to pass on to \code{rstan::stan} and
#' \code{rstan::sampling}.
#' @importFrom rstan stan sampling
#' @import Rcpp
#' @return A stanfit object.
#' @keywords internal
RunStanSampling <- function(stan_data, stan.model, bayes_options)
{

#	if (bayes_options$random_parameters == "fixed"){
#		pars = c("psi", "gamma", "alpha", "scale", "log_like", "sum_log_lik", "lp__")
#	}else
#		pars = c("mu", "tau", "L_Omega", "scale", "log_like", "sum_log_lik" )
#	if (is.null(pars))
#		pars <- stanParameters(stan.dat, keep.beta, stan.model)
#	init <- initialParameterValues(stan.dat)
	sampling(stan.model, data = stan_data,
					chains = bayes_options$n_chains,
					cores = bayes_options$n_cores,
#			 pars = pars,
			 iter = bayes_options$n_iterations, seed = bayes_options$seed,
			 control = list(max_treedepth = bayes_options$max_tree_depth,
			 			   adapt_delta = bayes_options$adapt_delta))
#			 init = init,
}

#stanParameters <- function(stan.dat, keep.beta, stan.model)
#{
#	full.covariance <- is.null(stan.dat$U)
#	multiple.classes <- !is.null(stan.dat$P)
#	has.covariates <- !is.null(stan.dat$covariates)
#
#	pars <- c("theta", "sigma")
#
#	if (multiple.classes)
#	{
#		if (has.covariates)
#			pars <- c(pars, "covariates_beta")
#		else
#			pars <- c(pars, "class_weights")
#	}else if (stan.model@model_name == "choicemodelRCdiag")
#		pars <- c("resp_fixed_coef", "sigma", "sig_rc",
#				  "log_likelihood")
#	if (keep.beta)
#		pars <- c(pars, "beta")
#
#	pars
#}


#' @title HierarchicalBayesMDCEV
#' @description Fit a MDCEV model HB
#' @param hb_options list of HB options
#' @param stan_data data for model
#' @inheritParams FitMDCEV
#' @param keep.samples default is FALSE,
#' @param include.stanfit default isTRUE,
#' @param show.stan.warnings = TRUE
#' @import rstan
#' @import dplyr
#' @export
#'

HierarchicalBayesMDCEV <- function(stan_data, hb_options,
								 initial.parameters,
								 keep.samples = FALSE,
								 include.stanfit = TRUE,
								 show.stan.warnings = TRUE)
{
	if (hb_options$n_iterations <= 0)
		stop("The specified number of iterations must be greater than 0.")

	# allows Stan chains to run in parallel on multiprocessor machines
	options(mc.cores = parallel::detectCores())

	# Create indices for individual level psi parameters
	indexes <- data_frame(individual = rep(1:stan_data$I, each = stan_data$J),
						  task = rep(1:stan_data$I, each = stan_data$J),
						  row = 1:(stan_data$I*stan_data$J)) %>%
		group_by_(task) %>%
		summarise(task_individual = first(individual),
				  start = first(row),
				  end = last(row))

	stan_data$start = indexes$start
	stan_data$end = indexes$end
	stan_data$task_individual = indexes$task_individual
	stan_data$task = indexes$task
	stan_data$IJ = stan_data$I * stan_data$J
	stan_data$lkj_shape = hb_options$hb.lkj.prior.shape

#	initial.parameters2 <- list(initial.parameters)#, initial.parameters,initial.parameters,initial.parameters)
#	initial.parameters2 <- list(list(scale = as.array(1, dim = 1)))#, initial.parameters,initial.parameters,initial.parameters)

#	has.covariates <- !is.null(stan_data$covariates)
	stan.model <- stanModel(hb_options$hb_random_parameters)

#	on.warnings <- GetStanWarningHandler(show.stan.warnings)
#	on.error <- GetStanErrorHandler()


#	InterceptExceptions(
#		{
	stan_fit <- RunStanSampling(stan_data, stan.model, hb_options)
#		}, warning.handler = on.warnings, error.handler = on.error)

#		result <- c(result, LogLikelihoodAndBIC(stan.fit, #n.hb.parameters,
#												stan.dat$R,
#												dat$n.questions.left.out,
#												dat$subset))
#		result

	result <- list()
	result$stan_fit <- stan_fit
#	n_parameters <- stan_fit[["par_dims"]][["mu"]] + stan_fit[["par_dims"]][["tau"]]
#	result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
#	result$effective.sample.size <- ess <- sum(weights)
#	n_parameters <- n_classes * n_parameters + n_classes - 1
#	result$bic <- -2 * result$log.likelihood + log(ess) * n_parameters
#	result$log.likelihood
	result
}

#' @title RunStanSampling
#' @description Wrapper function for \code{rstan:stan} and
#' \code{rstan:sampling} to run Stan HB analysis.
#' @inheritParams HierarchicalBayesMDCEV
#' @param stan.model Complied Stan model
#' @param ... Additional parameters to pass on to \code{rstan::stan} and
#' \code{rstan::sampling}.
#' @return A stanfit object.
#' @importFrom rstan stan sampling
#' @import Rcpp
#' @export
RunStanSampling <- function(stan_data, stan.model, hb_options)
{
#	if (is.null(pars))
#		pars <- stanParameters(stan.dat, keep.beta, stan.model)
#	init <- initialParameterValues(stan.dat)
	sampling(stan.model, data = stan_data, chains = hb_options$n_chains,
#			 pars = pars,
			 iter = hb_options$n_iterations, seed = hb_options$seed,
			 control = list(max_treedepth = hb_options$hb.max.tree.depth,
			 			   adapt_delta = hb_options$hb.adapt.delta))
#			 init = init,
}
#stan.model <- stanc("C:/Dropbox/Research/code/rmdcev/src/stan_files/mdcev.stan")

#stan("C:/Dropbox/Research/code/rmdcev/src/stan_files/mdcev.stan",
#	 data = stan_data, chains = n_chains, iter = n_iterations)

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


stanModel <- function(hb_random_parameters)
{
	covariates.error.msg <- paste0("Covariates are not currently implemented ",
								   "for the specified settings.")
#	if (n_classes == 1)
#	{
		if (hb_random_parameters == "fixed")
			stanmodels$mdcev
		else if (hb_random_parameters == "uncorr")
			stanmodels$mdcev_hb_uncorr
		else if (hb_random_parameters == "corr")
			stanmodels$mdcev_hb_corr
#	}
}

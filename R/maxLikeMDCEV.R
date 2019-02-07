#' @title maxlikeMDCEV
#' @description Fit a MDCEV model with MLE
maxlikeMDCEV <- function(data, initial.parameters = NULL,
						 seed,
						 mle_options)
{
	stan.model <- stanmodels$mdcev

	if (is.null(initial.parameters)){
		stan_fit <- optimizing(stan.model, data = data, as_vector = FALSE,
							   draws = mle_options$n_draws, hessian = mle_options$hessian)
	} else {
		stan_fit <- optimizing(stan.model, data = data, as_vector = FALSE, init = initial.parameters,
						   draws = mle_options$n_draws, hessian = mle_options$hessian)
	}

	result <- list()

	if (mle_options$keep_loglik == 0)
		stan_fit <- ReduceStanFitSize(stan_fit)

	result$stan_fit <- stan_fit
	n_parameters <- data$n_parameters
	result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
	result$effective.sample.size <- ess <- sum(data$weights)
#	n_parameters <- n_classes * n_parameters + n_classes - 1
	result$bic <- -2 * result$log.likelihood + log(ess) * n_parameters
	stan_fit <- result$stan_fit

	if (mle_options$n_classes > 1){
		result$mdcev_fit <- stan_fit
		result$mdcev_log.likelihood <- result$log.likelihood
		result$mdcev_bic <- result$bic

		init.par <- stan_fit$par

		# Extract the parameters to use as initial values for LC model
		# Need to ensure to replicate intial values for each class
		init.psi <- init.par$psi

		# add shift to psi values values
		init.shift <- seq(-0.02, 0.02, length.out = data$NPsi)
		for (i in 1:data$NPsi) {
			init.psi[i] <- init.psi[i] + init.shift[i]
		}

		init.psi <- matrix(init.psi, nrow=data$K,  ncol=length(init.psi), byrow=TRUE)

		init = list(psi = init.psi)

		if (data$fixed_scale == 0)
			init$scale <- rep(stan_fit$par[["scale"]], data$K)

		if (data$model_num == 1 || data$model_num == 3){
			init$alpha <- matrix(rep(init.par$alpha, data$K), nrow=data$K, ncol=1)
			init$gamma <- matrix(rep(init.par$gamma, data$K), nrow=data$K, ncol=data$J)
		} else if (data$model_num == 2){
			init$alpha <- matrix(rep(init.par$alpha, data$K), nrow=data$K, ncol=data$J)
		} else if (data$model_num == 4){
#			init$alpha <- matrix(rep(0, data$K), nrow=data$K, ncol=0)
			init$gamma <- init.par$gamma
		}

		stan.model <- stanmodels$mdcev_lc

		stan_fit <- optimizing(stan.model, data = data, as_vector = FALSE, init = init,
							   draws = mle_options$n_draws, hessian = mle_options$hessian)

		if (mle_options$keep_loglik == 0)
			stan_fit <- ReduceStanFitSize(stan_fit)

		result$stan_fit <- stan_fit
		n_parameters <- ncol(stan_fit[["hessian"]])
		result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
		result$effective.sample.size <- ess <- sum(data$weights)
		result$bic <- -2 * result$log.likelihood + log(ess) * n_parameters
		class_probabilities <- exp(t(stan_fit[["par"]][["theta"]]))
		colnames(class_probabilities) <- paste0("class", c(1:mle_options$n_classes))
		result$class_probabilities <- class_probabilities
	}

#	result$class.parameters <- pars$class.parameters
#	result$coef <- createCoefOutput(pars, data$par.names, data$all.names)
#	result$lca.data <- lca.data
	result
}


#' @title ReduceStanFitSize
#' @description This function reduces the size of the stan.fit object to reduce the time
#' it takes to return it from the R server.
#' @param stan_fit A stanfit object.
#' @return A stanfit object with a reduced size.
#' @export
ReduceStanFitSize <- function(stan_fit)
{
	# Replace stanmodel with a dummy as stanmodel makes the output many times larger,
	# and is not required for diagnostic plots.
	stan_fit[["par"]][["log_like"]] <- NULL
	stan_fit[["par"]][["log_like_all"]] <- NULL
	stan_fit[["theta_tilde"]] <- stan_fit[["theta_tilde"]][,1:ncol(stan_fit[["hessian"]])]
	stan_fit
}


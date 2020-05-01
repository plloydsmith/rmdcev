#' @title maxlikeMDCEV
#' @description Fit a MDCEV model with MLE
#' @param stan_data data for model fromatted from processMDCEVdata
#' @inheritParams mdcev
#' @param mle_options modeling options for MLE
#' @param parms_info information on parameters
maxlikeMDCEV <- function(stan_data, initial.parameters,
						 mle_options,
						 parms_info, ...) {

	if (is.null(initial.parameters) || mle_options$n_classes == 1){

		message("Using MLE to estimate MDCEV")

		# ensure single class used for base model
		stan_data_temp <- stan_data
		stan_data_temp$K <- 1
		stan_data_temp$L <- 0
		stan_data_temp$data_class <- matrix(0, stan_data$I, 0)

		if (is.null(initial.parameters)){
			stan_fit <- rstan::optimizing(stanmodels$mdcev, data = stan_data_temp, as_vector = FALSE, seed = mle_options$seed,
										  verbose = mle_options$print_iterations, iter = mle_options$max_iterations,
								   			draws = mle_options$n_draws, hessian = mle_options$hessian, ...)
		} else {
			stan_fit <- rstan::optimizing(stanmodels$mdcev, data = stan_data_temp, as_vector = FALSE, seed = mle_options$seed,
										  verbose = mle_options$print_iterations, iter = mle_options$max_iterations,
										  init = initial.parameters,
							   				draws = mle_options$n_draws, hessian = mle_options$hessian, ...)
		}

		if (mle_options$keep_loglik == 0)
			stan_fit <- ReduceStanFitSize(stan_fit, parms_info)

		result <- list()
		result$stan_fit <- stan_fit
		result$stan_fit$par[["theta"]] <- NULL
		result$stan_fit$par[["delta"]] <- NULL
		result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
		result$effective.sample.size <- sum(stan_data$weights)

	}

	if (mle_options$n_classes > 1){
		if (is.null(initial.parameters)) {
			result$mdcev_fit <- result$stan_fit
			result$mdcev_log.likelihood <- result$log.likelihood

			# Extract the parameters to use as initial values for LC model
			# Need to ensure to replicate intial values for each class
			init.par <- stan_fit$par
			K <- stan_data$K
			# add shift to psi values values
			init.psi <- init.par$psi
			init.shift <- seq(-0.02, 0.02, length.out = stan_data$NPsi)
			for (i in 1:stan_data$NPsi)
				init.psi[i] <- init.psi[i] + init.shift[i]

			init.psi <- matrix(init.psi, nrow=stan_data$K,  ncol=length(init.psi), byrow=TRUE)

			init = list(psi = init.psi)

			if (stan_data$fixed_scale1 == 0)
				init$scale <- rep(stan_fit$par[["scale"]], stan_data$K)

			if (stan_data$model_num == 1 || stan_data$model_num == 3){
				init$alpha <- matrix(rep(init.par$alpha, K), nrow=K, ncol=1)
				init$gamma <- matrix(rep(init.par$gamma, K), nrow=K, ncol=stan_data$J)
			} else if (stan_data$model_num == 2){
				init$alpha <- matrix(rep(init.par$alpha, K), nrow=K, ncol=stan_data$J+1)
			} else if (stan_data$model_num == 4){
				init$gamma <- matrix(rep(init.par$gamma, K), nrow=K, ncol=stan_data$J)
			} else if (stan_data$model_num == 5){
				init$alpha <- matrix(rep(init.par$alpha, K), nrow=K, ncol=1)
				init$phi <- matrix(rep(init.par$phi, K), nrow=K, ncol=stan_data$NPhi)
				if (stan_data$gamma_ascs == 1){
					init$gamma <- matrix(rep(init.par$gamma, K), nrow=K, ncol=stan_data$J)
				} else if (stan_data$gamma_ascs == 0){
					init$gamma <- matrix(rep(init.par$gamma, K), nrow=K, ncol=1)
				}
			}
		} else if (!is.null(initial.parameters)){
			init <- initial.parameters
		}
		message("Using MLE to estimate LC-MDCEV")

		stan_fit <- rstan::optimizing(stanmodels$mdcev, data = stan_data, as_vector = FALSE,
							   seed = mle_options$seed, init = init,
							   verbose = mle_options$print_iterations, iter = mle_options$max_iterations,
							   draws = mle_options$n_draws, hessian = mle_options$hessian)

		if (mle_options$keep_loglik == 0)
			stan_fit <- ReduceStanFitSize(stan_fit, parms_info)

		result$stan_fit <- stan_fit
		result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
		class_probabilities <- exp(t(stan_fit[["par"]][["theta"]]))
		colnames(class_probabilities) <- paste0("class", c(1:mle_options$n_classes))
		result$class_probabilities <- class_probabilities
	}
return(result)
}

#' @title ReduceStanFitSize
#' @description This function reduces the size of the stan.fit object
#' @param stan_fit A stanfit object.
#' @param parms_info information on parameters
#' @return A stanfit object with a reduced size.
ReduceStanFitSize <- function(stan_fit, parms_info) {
	stan_fit[["par"]][["log_like"]] <- NULL
	stan_fit[["theta_tilde"]] <- stan_fit[["theta_tilde"]][,1:parms_info$n_vars$n_parms_total]
	return(stan_fit)
}

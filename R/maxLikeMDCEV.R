#' @title maxlikeMDCEV
#' @description Fit a MDCEV model MLE
#' @import rstan
#' @export

maxlikeMDCEV <- function(dat,
						 n_classes = 1,
						 model_type,
						 n_draws = 50,
						 seed = 123,
						 initial.parameters = NULL,
#						 mle_tol = 0.0001,
						 hessian = TRUE)
{

	stanModel <- stanmodels$mdcev
	stan_fit <- optimizing(stanModel, data = dat,
						   draws = n_draws, hessian = hessian, ...)

	result <- list()
	result$stan_fit <- stan_fit
	n_parameters <- length(stan_fit[["par"]][["psi"]]) + length(stan_fit[["par"]][["gamma"]]) +
		length(stan_fit[["par"]][["alpha"]]) + length(stan_fit[["par"]][["scale"]])
	result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
	result$effective.sample.size <- ess <- sum(weights)
	n_parameters <- n_classes * n_parameters + n_classes - 1
	result$bic <- -2 * result$log.likelihood + log(ess) * n_parameters
	result$log.likelihood

	if (n_classes > 1){
		result$stan_1class <- stan_fit

		# Extract the parameters to use as initial values for LC model
		# Need to ensure to replicate intial values for each class
		init.psi <- unlist(stan_fit$par["psi"])

		# add shift to P_j values
		init.shift <- seq(-0.02, 0.02, length.out = dat$NPsi)
		for (i in 1:dat$NPsi) {
			init.psi[i] <- init.psi[i] + init.shift[i]
		}

		init.gamma <- unlist(stan_fit$par["gamma"])
		init.scale <- unlist(stan_fit$par["scale"])
		init.alpha <- unlist(stan_fit$par["alpha"])

		init.psi <- matrix(init.psi, nrow=dat$K,  ncol=length(init.psi),byrow=TRUE)

		init = list(psi = init.psi,
					# beta_m = matrix(rep(1,ncol(dat_membership))*(K-1),K-1,ncol(dat_membership)),
					alpha = matrix(rep(init.alpha, dat$K), nrow=dat$K, ncol=1),
					scale = rep(init.scale, dat$K),
					gamma = init.gamma)

		if (model_type == 4){
			init$alpha <- matrix(rep(init.alpha, dat$K), nrow=dat$K, ncol=0)
		}


		stanModel <- stanmodels$mdcev_lc
		stan_fit <- optimizing(stanModel, data = dat, init = init,
							   draws = n_draws, hessian = hessian, ...)

		result$stan_fit <- stan_fit
		n_parameters <- length(stan_fit[["par"]][["psi"]]) + length(stan_fit[["par"]][["gamma"]]) +
			length(stan_fit[["par"]][["alpha"]]) + length(stan_fit[["par"]][["scale"]]) + length(stan_fit[["par"]][["beta_m"]])
		result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
		result$effective.sample.size <- ess <- sum(weights)
		result$bic <- -2 * result$log.likelihood + log(ess) * n_parameters
		class_probabilities <- exp(tbl_df(t(stan_fit$par["theta"]$theta)))
		colnames(class_probabilities) = gsub("V", "class", colnames(class_probabilities))
		result$class_probabilities <- class_probabilities
	}

#	result$class.parameters <- pars$class.parameters
#	result$coef <- createCoefOutput(pars, dat$par.names, dat$all.names)
#	result$lca.data <- lca.data
	result
}

#' @title maxlikeMDCEV
#' @description Fit a MDCEV model MLE
#' @import rstan
#' @export

maxlikeMDCEV <- function(dat,
						 n_classes = 1,
						 n_draws = 50,
						 seed = 123,
						 initial.parameters = NULL,
#						 mle_tol = 0.0001,
						 hessian = TRUE,
						keep_loglik)
{

	stan.model <- stanmodels$mdcev
	stan_fit <- optimizing(stan.model, data = dat, as_vector = FALSE,
						   draws = n_draws, hessian = hessian)

	result <- list()

	if (keep_loglik == 0)
		stan_fit <- ReduceStanFitSize(stan_fit)

	result$stan_fit <- stan_fit
	n_parameters <- dat$n_parameters
	result$log.likelihood <- stan_fit[["par"]][["sum_log_lik"]]
	result$effective.sample.size <- ess <- sum(weights)
#	n_parameters <- n_classes * n_parameters + n_classes - 1
	result$bic <- -2 * result$log.likelihood + log(ess) * n_parameters
stan_fit <- result$stan_fit
#dat <- result$processed.data
#dat$K <- n_classes
	if (n_classes > 1){
		result$mdcev_fit <- stan_fit
		result$mdcev_log.likelihood <- result$log.likelihood
		result$mdcev_bic <- result$bic

		init.par <- stan_fit$par

		# Extract the parameters to use as initial values for LC model
		# Need to ensure to replicate intial values for each class
		init.psi <- init.par$psi

		# add shift to psi values values
		init.shift <- seq(-0.02, 0.02, length.out = dat$NPsi)
		for (i in 1:dat$NPsi) {
			init.psi[i] <- init.psi[i] + init.shift[i]
		}

		init.psi <- matrix(init.psi, nrow=dat$K,  ncol=length(init.psi), byrow=TRUE)

		init = list(psi = init.psi)

		if (dat$fixed_scale == 0)
			init$scale <- rep(stan_fit$par[["scale"]], dat$K)

		if (dat$model_type == 1 || dat$model_type == 3){
			init$alpha <- matrix(rep(init.par$alpha, dat$K), nrow=dat$K, ncol=1)
			init$gamma <- init.par$gamma
		} else if (dat$model_type == 2){
			init$alpha <- matrix(rep(init.par$alpha, dat$K), nrow=dat$K, ncol=dat$J)
		} else if (dat$model_type == 4){
#			init$alpha <- matrix(rep(0, dat$K), nrow=dat$K, ncol=0)
			init$gamma <- init.par$gamma
		}

		stan.model <- stanmodels$mdcev_lc

		stan_fit <- optimizing(stan.model, data = dat, as_vector = FALSE, init = init,
							   draws = n_draws, hessian = hessian)

		if (keep_loglik == 0)
			stan_fit <- ReduceStanFitSize(stan_fit)

		result$stan_fit <- stan_fit
		n_parameters <- ncol(stan_fit[["hessian"]])
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


#' @title ReduceStanFitSize
#' @description This function reduces the size of the stan.fit object to reduce the time
#' it takes to return it from the R server.
#' @param stan.fit A stanfit object.
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


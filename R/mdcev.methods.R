############################
# S3 method for rmdcev package
#############################

#' @rdname mdcev
#' @method print mdcev
#' @param x,object an object of class `mdcev`
#' @param digits the number of digits,
#' @param width the width of the printing,
#' @export
print.mdcev <- function(x, digits = max(3, getOption("digits") - 3),
						width = getOption("width"), ...){

	cat("\nCoefficients:\n")
	print.default(format(stats::coef(x), digits = digits), print.gap = 2,
				  quote = FALSE)
	cat("\n")
	invisible(x)
}


#' @rdname mdcev
#' @method summary mdcev
#' @param printCI set to TRUE to print 95\% confidence intervals
#' @export
summary.mdcev <- function(object, printCI=FALSE, ...){
#object <- output
	if(object$algorithm == "MLE"){

		if(object$std_errors == "deltamethod"){

			coefs <- object$stan_fit$par
			coefs$sum_log_lik <- coefs$theta <- NULL
			coefs <- unlist(coefs)
			cov_mat <- solve(-object$stan_fit$hessian)
			std_err <- sqrt(diag(cov_mat))
			parms <- object[["parms_info"]][["parm_names"]][["all_names"]]
			output <- as_tibble(cbind(parms, coefs, std_err)) %>%
				mutate(coefs = as.numeric(coefs),
					   std_err = as.numeric(std_err),
					   Estimate = round(coefs, 3),
					   Std.err = round(ifelse(stringr::str_detect(parms, "alpha"), std_err*exp(-coefs)/((1+exp(-coefs)))^2,
					   					   ifelse(stringr::str_detect(parms, "gamma|scale|phi"), std_err*coefs, std_err)),3),
					   z.stat = round(coefs / Std.err,2),
					   ci_lo95 = round(coefs - 1.96 * Std.err,3),
					   ci_hi95 = round(coefs + 1.96 * Std.err,3)) %>%
				select(parms, Estimate, Std.err, z.stat, ci_lo95, ci_hi95)

		} else if(object$std_errors == "mvn"){

			# Get parameter estimates in matrix form
			output <- as_tibble(object[["stan_fit"]][["theta_tilde"]]) %>%
				dplyr::select(-tidyselect::starts_with("log_like"), -tidyselect::starts_with("sum_log_lik"))

			names(output)[1:object$parms_info$n_vars$n_parms_total] <- object$parms_info$parm_names$all_names

			output <- output %>%
				tibble::rowid_to_column("sim_id") %>%
				tidyr::gather(parms, value, -sim_id)

			output <- output %>%
				mutate(parms = factor(parms,levels=unique(parms))) %>%
				group_by(parms) %>%
				summarise(Estimate = round(mean(value),3),
						  Std.err = round(stats::sd(value),3),
						  z.stat = round(mean(value) / stats::sd(value),2),
						  ci_lo95 = round(stats::quantile(value, 0.025),3),
						  ci_hi95 = round(stats::quantile(value, 0.975),3))
		}
	} else if (object$algorithm == "Bayes"){

		# Get parameter estimates in matrix form
		output <- rstan::extract(object$stan_fit, permuted = TRUE, inc_warmup = FALSE) %>%
			as.data.frame() %>%
			dplyr::select(-tidyselect::starts_with("log_like"), -tidyselect::starts_with("sum_log_lik"),
						  -tidyselect::starts_with("tau_unif"), -.data$lp__)

		if (object$random_parameters == "fixed")
			names(output)[1:object$parms_info$n_vars$n_parms_total] <- object$parms_info$parm_names$all_names

		output <- output %>%
			tibble::rowid_to_column("sim_id") %>%
			tidyr::gather(parms, value, -sim_id)

		bayes_extra <- as_tibble(rstan::summary(object$stan_fit)$summary) %>%
			mutate(parms = row.names(rstan::summary(object$stan_fit)$summary))

		if (object$random_parameters == "fixed"){

			parm.names <- unique(output$parms)

			bayes_extra <- bayes_extra %>%
				filter(!grepl(c("log_lik"), parms)) %>%
				filter(!grepl(c("lp_"), parms)) %>%
				select(parms, n_eff, Rhat) %>%
				mutate(parms = factor(parm.names,levels=unique(parm.names)),
					   n_eff = round(as.numeric(n_eff), 0),
					   Rhat = round(as.numeric(Rhat), 2))

		} else if (object$random_parameters != "fixed"){

			output_non_mu <- output %>%
				filter(stringr::str_detect(parms, "gamma|alpha|tau|scale"))

			output <- output %>%
				filter(stringr::str_detect(parms, "mu")) %>%
				bind_rows(output_non_mu) %>%
				arrange(sim_id)

			output$parms <- rep(c(object[["parms_info"]][["parm_names"]][["all_names"]],
								  object[["parms_info"]][["parm_names"]][["sd_names"]]), max(output$sim_id))

			# Transform estimates
			if(object[["stan_data"]][["gamma_nonrandom"]]==0){
				output <- output %>%
					mutate(value = ifelse(grepl(c("gamma"), parms), exp(value), value))
			}
			if(object[["stan_data"]][["alpha_nonrandom"]]==0){
				output <- output %>%
					mutate(value = ifelse(grepl(c("alpha"), parms), 1 / (1 + exp(-value)), value))
			}

			bayes_extra_non_mu <- bayes_extra %>%
				filter(grepl(c("gamma|alpha|scale|tau"), parms)) %>%
				filter(!grepl(c("tau_unif"), parms))

			bayes_extra <- bayes_extra %>%
				filter(grepl(c("mu"), parms)) %>%
				bind_rows(bayes_extra_non_mu)

			bayes_extra$parms <- c(object[["parms_info"]][["parm_names"]][["all_names"]],
								   object[["parms_info"]][["parm_names"]][["sd_names"]])

			bayes_extra <- bayes_extra %>%
				select(parms, n_eff, Rhat) %>%
				mutate(parms = factor(parms,levels=unique(parms)),
					   n_eff = round(as.numeric(n_eff), 0),
					   Rhat = round(as.numeric(Rhat), 2))
		}

		output <- output %>%
			mutate(parms = factor(parms,levels=unique(parms))) %>%
			group_by(parms) %>%
			summarise(Estimate = round(mean(value),3),
					  Std.dev = round(stats::sd(value),3),
					  z.stat = round(mean(value) / stats::sd(value),2),
					  ci_lo95 = round(stats::quantile(value, 0.025),3),
					  ci_hi95 = round(stats::quantile(value, 0.975),3)) %>%
			left_join(bayes_extra, by = "parms")
	}

	output <- as.data.frame(output)

	rownames(output) <- c(as.character(output$parms))

	dropcolumns=NULL
	if(printCI==FALSE) dropcolumns = c(dropcolumns,5,6)
	dropcolumns = unique(dropcolumns)
	if(length(dropcolumns)>0) output = output[,-dropcolumns, drop=FALSE]

	output$parms <- NULL

	object$CoefTable    <- output
	class(object)       <- c("summary.mdcev", "mdcev")
	return(object)
}

#' @rdname mdcev
#' @method print summary.mdcev
#' @export
print.summary.mdcev <- function(x,...){
	#	x <- output

	rmdcevVersion <- tryCatch(utils::packageDescription("rmdcev", fields = "Version"),
							  warning=function(w) return("rmdcev"),
							  error=function(e) return("rmdcev"))

	cat("Model run using rmdcev for R, version", rmdcevVersion,"\n")
	cat("Estimation method                : ", x$algorithm, "\n", sep="")
	cat("Model type                       : ", x$model," specification", "\n", sep="")
	cat("Number of classes                : ", x$n_classes, "\n", sep="")
	cat("Number of individuals            : ", x$n_individuals,"\n", sep="")
	cat("Number of non-numeraire alts     : ", x$stan_data$J,"\n", sep="")
	cat("Estimated parameters             : ", x$parms_info$n_vars$n_parms_total,"\n", sep="")
	cat("LL                               : ", round(x$log.likelihood,2),"\n", sep="")

	if(x$algorithm == "MLE"){
		cat("AIC                              : ", round(x$aic,2),"\n", sep="")
		cat("BIC                              : ", round(x$bic,2),"\n", sep="")
		if(x$std_errors == "deltamethod"){
			cat("Standard errors calculated using : ", "Delta method","\n", sep="")

		} else if(x$std_errors == "mvn"){
			cat("Standard errors calculated using : ", x$n_draws," MVN draws", "\n", sep="")
		}

		if(x$stan_fit$return_code==0){
			converge <- "successful convergence"
		} else if(x$stan_fit$return_code==1){
			converge <- "unsuccessful convergence"
		}
		cat("Exit of MLE                      : ", converge,"\n", sep="")

	} else if(x$algorithm == "Bayes"){
		if(x$random_parameters != "fixed"){
			cat("Random parameters                : ", x$random_parameters,"elated random parameters","\n", sep="")
		}
		cat("Number of chains                 : ", x[["stan_fit"]]@sim[["chains"]],"\n", sep="")
		cat("Number of warmup draws per chain : ", x[["stan_fit"]]@sim[["warmup"]],"\n", sep="")
		cat("Total post-warmup sample         : ", x[["stan_fit"]]@sim[["chains"]]*(x[["stan_fit"]]@sim[["iter"]]-x[["stan_fit"]]@sim[["warmup"]]),"\n", sep="")
	}
	tmpH <- floor(x$time.taken/60^2)
	tmpM <- floor((x$time.taken-tmpH*60^2)/60)
	tmpS <- round(x$time.taken-tmpH*60^2-tmpM*60,2)
	timeTaken <- paste(formatC(tmpH,width=2,format='d',flag=0),
					   formatC(tmpM,width=2,format='d',flag=0),
					   tmpS,sep=':')
	cat("Time taken (hh:mm:ss)            : ",timeTaken,"\n", sep="")

	cat("\nAverage consumption of non-numeraire alternatives:\n")
	mean_consumption <-  round(colMeans(x$stan_data$quant_j),2)
	names(mean_consumption) <-x$parms_info$alt_names
	print(mean_consumption )
	cat("\n")

	#	cat("\nPsi specification:\n")
	#	cat(paste(x$psi_formula, sep = "\n", collapse = "\n"), "\n", sep = "")

	if (x$n_classes > 1){
		cat("\nClass average probabilities:\n")
		print(round(colMeans(x[["class_probabilities"]]),2))
	}

	cat("Parameter estimates -------------------------------- ","\n")
	if(nrow(x$CoefTable)>options("max.print")) options(max.print=nrow(x$CoefTable)+100)
	print(x$CoefTable)
	if(x$stan_data$fixed_scale1 == 1)
		cat("Note: Scale parameter fixed to 1.",'\n')

	if(x$model == "gamma"){
		cat("Note: All non-numeraire alpha's fixed to 0.",'\n')
	}else if (x$model == "alpha"){
		cat("Note: All non-numeraire gamma's fixed to 1.",'\n')
	}else if (x$model == "hybrid"){
		cat("Note: Alpha parameter is equal for all alternatives.",'\n')
	}else if (x$model == "hybrid0")
		cat("Note: All alpha parameters fixed to 1e-3.",'\n')

	if(x$stan_data$trunc_data == 1)
		cat("Note: Estimation accounts for truncated form of data.",'\n')

	if(x$random_parameters == "corr"){
		cat("Note: The full covariance matrix can be accessed using the print(output$stan_fit, pars = 'Sigma') command where output is name of model output", '\n')
	}

	if(x$algorithm == "Bayes"){
		cat("Note from Rstan: 'For each parameter, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence, Rhat=1)'", '\n')
	}
	if(x$n_classes > 1)
		cat("Note: The membership equation parameters for class 1 are normalized to 0.",'\n')

	cat("\n")

	invisible(x)
}

#' @export
logLik.mdcev <- function(object,...){
	object$log.likelihood
}

#' @export
coef.mdcev <- function(object, ...){
	result <- object$stan_fit[["par"]]
	# first remove the fixed coefficients if required
	result$theta <- NULL
	result$sum_log_lik <- NULL
	ncoefs <- names(result)
	selcoef <- 1:length(result)
	result[selcoef]
}

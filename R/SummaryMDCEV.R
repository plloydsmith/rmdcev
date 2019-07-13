#' @title SummaryMDCEV
#' @description Prints mdcev estimation results to console. The output format is borrowed from the Apollo package
#' @param model Model object returned by function \link{FitMDCEV}.
#' @param printCI set to TRUE to print 95\% confidence intervals
#' @return A matrix of coefficients, s.d. and z-tests (invisible)
#' @export
#' @examples
#' \donttest{
#' data(data_rec, package = "rmdcev")
#' mdcev_est <- FitMDCEV(psi_formula = ~ 1,
#' data = subset(data_rec, id < 500),
#' model = "hybrid0",
#' algorithm = "MLE")
#'
#' SummaryMDCEV(mdcev_est)
#' }
SummaryMDCEV <- function(model, printCI = FALSE){
#model <- mdcev_rp
#	if (model$random_parameters == "corr"){
#		stop("SummaryMDCEV not set up for correlated random parameter models. Use print/traceplot on model.fit$stan_fit to examine output", "\n")
#	}

	rmdcevVersion <- tryCatch(utils::packageDescription("rmdcev", fields = "Version"),
							  warning=function(w) return("rmdcev"),
							  error=function(e) return("rmdcev"))

	cat("Model run using rmdcev for R, version", rmdcevVersion,"\n")
	cat("Estimation method                : ", model$algorithm, "\n", sep="")
	cat("Model type                       : ", model$model," specification", "\n", sep="")
	cat("Number of classes                : ", model$n_classes, "\n", sep="")
	cat("Number of individuals            : ", model$n_individuals,"\n", sep="")
	cat("Number of non-numeraire alts     : ", model$stan_data$J,"\n", sep="")
	cat("Estimated parameters             : ", model$parms_info$n_vars$n_parms_total,"\n", sep="")
	cat("LL                               : ", round(model$log.likelihood,2),"\n", sep="")

	if(model$algorithm == "MLE"){
	cat("AIC                              : ", round(model$aic,2),"\n", sep="")
	cat("BIC                              : ", round(model$bic,2),"\n", sep="")
		if(model$std_errors == "deltamethod"){
    cat("Standard errors calculated using : ", "Delta method","\n", sep="")

		} else if(model$std_errors == "mvn"){
    cat("Standard errors calculated using : ", model$n_draws," MVN draws", "\n", sep="")
		}

		if(model$stan_fit$return_code==0){
			converge <- "successful convergence"
		} else if(model$stan_fit$return_code==0){
			converge <- "unsuccessful convergence"
		}
	cat("Exit of MLE                      : ", converge,"\n", sep="")

	} else if(model$algorithm == "Bayes"){
		cat("Random parameters                : ", model$random_parameters,"elated random parameters","\n", sep="")
		cat("Number of chains                 : ", model[["stan_fit"]]@sim[["chains"]],"\n", sep="")
		cat("Number of warmup draws per chain : ", model[["stan_fit"]]@sim[["warmup"]],"\n", sep="")
		cat("Total post-warmup sample         : ", model[["stan_fit"]]@sim[["chains"]]*(model[["stan_fit"]]@sim[["iter"]]-model[["stan_fit"]]@sim[["warmup"]]),"\n", sep="")
	}
tmpH <- floor(model$time.taken/60^2)
tmpM <- floor((model$time.taken-tmpH*60^2)/60)
tmpS <- round(model$time.taken-tmpH*60^2-tmpM*60,2)
timeTaken <- paste(formatC(tmpH,width=2,format='d',flag=0),
				   formatC(tmpM,width=2,format='d',flag=0),
				   tmpS,sep=':')
	cat("Time taken (hh:mm:ss)            : ",timeTaken,"\n", sep="")

	if(model$algorithm == "MLE"){

		if(model$std_errors == "deltamethod"){

			coefs <- model$stan_fit$par
			coefs$sum_log_lik <- coefs$theta <- NULL
			coefs <- unlist(coefs)
			cov_mat <- solve(-model$stan_fit$hessian)
			std_err <- sqrt(diag(cov_mat))
			parms <- model[["parms_info"]][["parm_names"]][["all_names"]]
			output <- tbl_df(cbind(parms, coefs, std_err)) %>%
				mutate(coefs = as.numeric(coefs),
					   std_err = as.numeric(std_err),
					   Estimate = round(coefs, 3),
					   Std.err = round(ifelse(stringr::str_detect(parms, "alpha"), std_err*exp(-coefs)/((1+exp(-coefs)))^2,
										 ifelse(stringr::str_detect(parms, "gamma|scale"), std_err*coefs, std_err)),3),
					   z.stat = round(coefs / Std.err,2),
					   ci_lo95 = round(coefs - 1.96 * Std.err,3),
					   ci_hi95 = round(coefs + 1.96 * Std.err,3)) %>%
				select(parms, Estimate, Std.err, z.stat, ci_lo95, ci_hi95)

		} else if(model$std_errors == "mvn"){

		output <- model$est_pars %>%
			mutate(parms = factor(parms,levels=unique(parms))) %>%
			group_by(parms) %>%
			summarise(Estimate = round(mean(value),3),
					  Std.err = round(stats::sd(value),3),
					  z.stat = round(mean(value) / stats::sd(value),2),
					  ci_lo95 = round(stats::quantile(value, 0.025),3),
					  ci_hi95 = round(stats::quantile(value, 0.975),3))
		}
	} else if (model$algorithm == "Bayes"){
		output <- model$est_pars

		bayes_extra <- tbl_df(rstan::summary(model$stan_fit)$summary) %>%
			mutate(parms = row.names(rstan::summary(model$stan_fit)$summary))

		if (model$random_parameters == "fixed"){

			parm.names <- unique(output$parms)

			bayes_extra <- bayes_extra %>%
				filter(!grepl(c("log_lik"), parms)) %>%
				filter(!grepl(c("lp_"), parms)) %>%
				select(parms, n_eff, Rhat) %>%
				mutate(parms = factor(parm.names,levels=unique(parm.names)),
					   n_eff = round(as.numeric(n_eff), 0),
					   Rhat = round(as.numeric(Rhat), 2))

		} else if (model$random_parameters != "fixed"){

			output_non_mu <- output %>%
				filter(stringr::str_detect(parms, "gamma|alpha|tau|scale"))

			output <- output %>%
				filter(stringr::str_detect(parms, "mu")) %>%
				bind_rows(output_non_mu) %>%
				arrange(sim_id)

			output$parms <- rep(c(model[["parms_info"]][["parm_names"]][["all_names"]],
						 model[["parms_info"]][["parm_names"]][["sd_names"]]), max(output$sim_id))

			# Transform estimates
			if(model[["stan_data"]][["gamma_fixed"]]==0){
			output <- output %>%
				mutate(value = ifelse(grepl(c("gamma"), parms), exp(value), value))
			}
			if(model[["stan_data"]][["alpha_fixed"]]==0){
				output <- output %>%
					mutate(value = ifelse(grepl(c("alpha"), parms), 1 / (1 + exp(-value)), value))
			}

			bayes_extra_non_mu <- bayes_extra %>%
				filter(grepl(c("gamma|alpha|scale|tau"), parms)) %>%
				filter(!grepl(c("tau_unif"), parms))

			bayes_extra <- bayes_extra %>%
				filter(grepl(c("mu"), parms)) %>%
				bind_rows(bayes_extra_non_mu)

			bayes_extra$parms <- c(model[["parms_info"]][["parm_names"]][["all_names"]],
								  model[["parms_info"]][["parm_names"]][["sd_names"]])

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
					  Std.err = round(stats::sd(value),3),
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

	cat("\nAverage consumption of non-numeraire alternatives:\n")
	mean_consumption <-  round(colMeans(model$stan_data$quant_j),2)
	names(mean_consumption) <-c(1:model$stan_data$J)
	print(mean_consumption )
	cat("\n")

#	cat("\nPsi specification:\n")
#	cat(paste(model$psi_formula, sep = "\n", collapse = "\n"), "\n", sep = "")

	if (model$n_classes > 1){
		cat("\nMembership equation:\n")
		cat(paste(model$lc_formula, sep = "\n", collapse = "\n"), "\n", sep = "")
		cat("\nClass average probabilities:\n")
		print(round(colMeans(model[["class_probabilities"]]),2))
	}

	output$parms <- NULL
	cat("Parameter estimates -------------------------------- ","\n")
	if(nrow(output)>options("max.print")) options(max.print=nrow(output)+100)
	print(output)
	if(model$stan_data$fixed_scale1 == 1)
		cat("Note: Scale parameter fixed to 1.",'\n')

	if(model$model == "gamma"){
		cat("Note: All non-numeraire alpha's fixed to 0.",'\n')
	}else if (model$model == "alpha"){
		cat("Note: All non-numeraire gamma's fixed to 1.",'\n')
	}else if (model$model == "hybrid"){
		cat("Note: Alpha parameter is equal for all goods.",'\n')
	}else if (model$model == "hybrid0")
		cat("Note: All alpha parameters fixed to 1e-3.",'\n')

	if(model$stan_data$trunc_data == 1)
		cat("Note: Estimation accounts for truncated form of data.",'\n')

	if(model$random_parameters == "corr"){
		cat("Note: Full covariance matrix can be accessed using the print(model_est, pars = 'Sigma') command", '\n')
	}

	if(model$algorithm == "Bayes"){
		cat("Note from Rstan: 'For each parameter, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence, Rhat=1)'", '\n')
	}
	if(model$n_classes > 1)
		cat("Note: The membership equation parameters for class 1 are normalized to 0.",'\n')

		cat("\n")

	invisible(output)
}

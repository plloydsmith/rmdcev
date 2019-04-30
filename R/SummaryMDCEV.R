#' @title SummaryMDCEV
#' @description Prints mdcev estimation results to console. The output format is borrowed from the Apollo package
#' @param model Model object returned by function \link{FitMDCEV}.
#' @param printCI set to TRUE to print 95\% confidence intervals
#' @return A matrix of coefficients, s.d. and z-tests (invisible)
#' @export
SummaryMDCEV <- function(model, printCI = FALSE){
#model <- result
	if (model$random_parameters == "corr"){
		stop("SummaryMDCEV not set up for correlated random parameter models. Use print/traceplot on model.fit$stan_fit to examine output", "\n")
	}

	rmdcevVersion <- tryCatch(utils::packageDescription("rmdcev", fields = "Version"),
							  warning=function(w) return("rmdcev"),
							  error=function(e) return("rmdcev"))

	cat("Model run using rmdcev for R, version", rmdcevVersion,"\n")
	cat("Estimation method                : ", model$algorithm, "\n", sep="")
	cat("Model type                       : ", model$model," specification", "\n", sep="")
	cat("Number of classes                : ", model$n_classes, "\n", sep="")
	cat("Number of individuals            : ", model$n_respondents,"\n", sep="")
	cat("Number of non-numeraire alts     : ", model$stan_data$J,"\n", sep="")
	cat("Estimated parameters             : ", model$parms_info$n_vars$n_parms_total,"\n", sep="")
	cat("LL                               : ", round(model$log.likelihood,2),"\n", sep="")

	if(model$algorithm == "MLE"){
	cat("AIC                              : ", round(model$aic,2),"\n", sep="")
	cat("BIC                              : ", round(model$bic,2),"\n", sep="")
		if(model$std_errors == "deltamethod"){
    cat("Standard errors calculated using : ", "Delta method","\n", sep="")

		} else if(model$std_errors == "mvn"){
    cat("Standard errors calculated using : ", model$n_draws,"MVN draws", "\n", sep="")
		}

		if(model$stan_fit$return_code==0){
			converge <- "successful convergence"
		} else if(model$stan_fit$return_code==0){
			converge <- "unsuccessful convergence"
		}
	cat("Exit of MLE                      : ", converge,"\n", sep="")

	} else if(model$algorithm == "Bayes"){
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
		bayes_extra <- tbl_df(summary(model$stan_fit)$summary) %>%
			mutate(parms = row.names(summary(model$stan_fit)$summary))

		if (model$random_parameters == "fixed"){

			bayes_extra <- bayes_extra %>%
				filter(!grepl(c("log_lik"), parms)) %>%
				filter(!grepl(c("lp_"), parms)) %>%
				select(n_eff, Rhat) %>%
				mutate(n_eff = round(n_eff, 0),
					   Rhat = round(Rhat, 2))

		} else if (model$random_parameters == "uncorr"){

			output <- output %>%
				filter(stringr::str_detect(parms, "mu|tau|scale")) %>%
				arrange(sim_id)

			output$parms <- rep(c(model[["parms_info"]][["parm_names"]][["all_names"]],
						 model[["parms_info"]][["parm_names"]][["sd_names"]]), max(output$sim_id))
			bayes_extra <- bayes_extra %>%
					filter(grepl(c("mu|scale|tau"), parms)) %>%
					filter(!grepl(c("tau_unif"), parms)) %>%
					select(n_eff, Rhat) %>%
					mutate(n_eff = round(n_eff, 0),
						   Rhat = round(Rhat, 2))
			}

		output <- output %>%
			mutate(parms = factor(parms,levels=unique(parms))) %>%
			group_by(parms) %>%
			summarise(Estimate = round(mean(value),3),
					  Std.err = round(stats::sd(value),3),
					  z.stat = round(mean(value) / stats::sd(value),2),
					  ci_lo95 = round(stats::quantile(value, 0.025),3),
					  ci_hi95 = round(stats::quantile(value, 0.975),3))

		output <- cbind(output, bayes_extra)
	}

	output <- as.data.frame(output)

	rownames(output) <- c(as.character(output$parms))

	dropcolumns=NULL
	if(printCI==FALSE) dropcolumns = c(dropcolumns,5,6)
	dropcolumns = unique(dropcolumns)
	if(length(dropcolumns)>0) output = output[,-dropcolumns, drop=FALSE]

	cat("\nAverage consumption of non-numeraire alternatives:\n")
	mean_consumption <-  round(colMeans(model$stan_data$j_quant),2)
	names(mean_consumption) <-c(1:model$stan_data$J)
	print(mean_consumption )
	cat("\n")

	cat("\nPsi specification:\n")
	cat(paste(model$psi_formula, sep = "\n", collapse = "\n"), "\n", sep = "")

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
	if(model$stan_data$fixed_scale == 1)
		cat("Note: Scale parameter fixed at 1.",'\n')

	if(model$model == "les"){
		cat("Note: All non-numeraire alpha's set to 0.",'\n')
	}else if (model$model == "alpha"){
		cat("Note: All non-numeraire gamma's set to 1.",'\n')
	}else if (model$model == "gamma"){
		cat("Note: The alpha parameter is equal for all goods.",'\n')
	}else if (model$model == "gamma0")
		cat("Note: The alpha parameter is equal for all goods and set to 1e-6.",'\n')

	if(model$stan_data$trunc_data == 1)
		cat("Note: Estimation accounts for truncated form of data.",'\n')

	if(model$algorithm == "Bayes"){
		cat("Note: For each parameter, n_eff is a crude measure of effective sample size, and Rhat is the potential scale reduction factor on split chains (at convergence, Rhat=1)", '\n')
	}
	if(model$n_classes > 1)
		cat("Note: The membership equation parameters for class 1 are normalized to 0.",'\n')

		cat("\n")

	invisible(output)
}

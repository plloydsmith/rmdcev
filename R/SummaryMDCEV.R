#' @title SummaryMDCEV
#' @description Prints mdcev estimation results to console. The output format is borrowed from the Apollo package
#' @param model Model object returned by function \link{FitMDCEV}.
#' @param printCI set to TRUE to print 95\% confidence intervals
#' @return A matrix of coefficients, s.d. and z-tests (invisible)
#' @export
SummaryMDCEV <- function(model, printCI = FALSE){
#model <- result
	rmdcevVersion <- tryCatch(utils::packageDescription("rmdcev", fields = "Version"),
							  warning=function(w) return("rmdcev"),
							  error=function(e) return("rmdcev"))
	cat("Model run using rmdcev for R, version", rmdcevVersion,"\n")
#	cat("Model run at                     : ", paste(model$start.time[3] ),"\n", sep="")
	cat("Estimation method                : ", model$algorithm, "\n", sep="")
#	if(!apollo_control$HB) cat("Model diagnosis                  : ",model$message,"\n", sep="")
	cat("Model type                       : ", model$model," specification", "\n", sep="")
	cat("Number of classes                : ", model$n_classes, "\n", sep="")
	cat("Number of individuals            : ", model$n_respondents,"\n", sep="")
	cat("Number of non-numeraire alts     : ", model$stan_data$J,"\n", sep="")
	#	cat("Number of observations           : ", model$n_respondents,"\n", sep="")
	cat("Estimated parameters             : ", model$stan_data$n_parameters,"\n", sep="")
	cat("LL                               : ",round(model$log.likelihood,2),"\n")
	cat("AIC                              : ",round(model$aic,2),"\n")
	cat("BIC                              : ",round(model$bic,2),"\n")
	cat("Number of MVN simulation draws   : ",model$n_draws,"\n")
	if(model$algorithm == "MLE"){
		if(model$stan_fit$return_code==0){
			converge <- "successful convergence"
		} else if(model$stan_fit$return_code==0){
			converge <- "unsuccessful convergence"
		}
		cat("Exist of MLE                     : ",converge,"\n")
	}

tmpH <- floor(model$time.taken/60^2)
tmpM <- floor((model$time.taken-tmpH*60^2)/60)
tmpS <- round(model$time.taken-tmpH*60^2-tmpM*60,2)
timeTaken <- paste(formatC(tmpH,width=2,format='d',flag=0),
				   formatC(tmpM,width=2,format='d',flag=0),
				   tmpS,sep=':')
	cat("Time taken (hh:mm:ss)            : ",timeTaken,"\n")
#cat("Iterations                       : ",model$nIter,"\n")
#cat("Number of cores used             : ",model$apollo_control$nCores,"\n")



	output <- model$est_pars %>%
		mutate(parms = gsub("\\[|\\]|\\.", "", parms)) %>%
		mutate(parms = factor(parms,levels=unique(parms))) %>%
		#	mutate(parms = gsub("\\.", "", parms)) %>%
		group_by(parms) %>%
		summarise(Estimate = round(mean(value),3),
				  Std.err = round(stats::sd(value),3),
				  z.stat = round(mean(value) / stats::sd(value),2),
				  cl95_lo = round(stats::quantile(value, 0.025),3),
				  cl95_hi = round(stats::quantile(value, 0.975),3))

	if(model$n_classes > 1){
	psi.names <- cbind("psi",1:ncol(model$stan_data[["dat_psi"]]),
					   paste0(".",colnames(model$stan_data[["dat_psi"]])))
	colnames(psi.names) <- c("parms", "parm_num", "names")
	psi.names <- tbl_df(psi.names)

	beta_m.names <- cbind("beta_m",1:ncol(model$stan_data[["data_class"]]),
					   paste0(".",colnames(model$stan_data[["data_class"]])))
	colnames(beta_m.names) <- c("parms", "parm_num", "names_b")
	beta_m.names <- tbl_df(beta_m.names)

		output <- suppressWarnings(output %>%
		separate(parms, into = c("parms", "parm_num"), sep = ",") %>%
		separate(parms, into = c("parms", "class"), sep = -1) %>%
			mutate(parm_num = ifelse(is.na(parm_num), "", parm_num),
				class = ifelse(parms == "beta_m", as.numeric(class)+1, as.numeric(class))) %>%
			arrange(class) %>%
			left_join(psi.names, by = c("parms", "parm_num")) %>%
			left_join(beta_m.names, by = c("parms", "parm_num")) %>%
			mutate(parm_num = ifelse(parms == "psi", names, parm_num),
				   parm_num = ifelse(parms == "beta_m", names_b, parm_num),
			parms = paste0("class", class, ".", parms, parm_num)) %>%
			select(-class, -parm_num, -names, -names_b))
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

	if(model$n_classes > 1)
		cat("Note: The membership equation parameters for class 1 are normalized to 0.",'\n')

		cat("\n")

	invisible(output)
}

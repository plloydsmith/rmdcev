#' @title SimulateWTP
#' @description Simulate WTP for MDCEV model
#' @param df_indiv Prepared individual level data from PrepareSimulationData
#' @param df_common Prepared common data from PrepareSimulationData
#' @param nerrs Number of error draws for welfare analysis
#' @param cond_error Choose whether to draw errors conditional on actual demand or not.
#' Conditional error draws (=1) or unconditional error draws.
#' @param algo_gen Type of algorhitm for simulation. algo_gen = 0 for Hybrid Approach (i.e. constant alphas,
#' only model 3/4) alog_gen = 1 for General approach (i.e. heterogeneous alpha's, all models)
#' @return wtp a list for each individual holding a nsims x npols matrix of wtp
#' @importFrom stats quantile sd
#' @export
SimulateWTP <- function(df_indiv, df_common,
						nerrs = 30,
						cond_error = 1,
						algo_gen = NULL){

	start.time <- proc.time()

	# Checks on simulation options
	model_num <- df_common$model_num

	if (!is.null(algo_gen)){
		if (model_num < 3 && algo_gen == 0){
			warning("Can't use hybrid algorithm with model_num = 1 or 2. Changing to general approach.")
			algo_gen <- 1
		}
	} else if (is.null(algo_gen)) {
		if (model_num == 3 || model_num == 4)
			algo_gen <- 0
		else if (model_num == 1 || model_num == 2)
			algo_gen <- 1
	}

	if (algo_gen == 1) {
		message("Using general approach to simulation")
	} else if (algo_gen == 0){
		message("Using hybrid approach to simulation")
	}
	# Organize options in list
	sim_options <- list(nerrs = nerrs,
						cond_error = cond_error,
						algo_gen = algo_gen,
						model_num = model_num)

	if(df_common$n_classes == 1){
		wtp <- StanWTP(df_indiv, df_common, sim_options)

	} else if(df_common$n_classes > 1){
		wtp <- map2(df_indiv, df_common, StanWTP, sim_options)
		names(wtp) <- paste0("class", c(1:df_common$n_classes))
	}

	time <- proc.time() - start.time
	n_simulations <- length(unlist(wtp)) * nerrs
	cat("\n", formatC(n_simulations, format = "e", digits = 2), "simulations finished in", round(time[3]/60, 2), "minutes.",
		"(",round(n_simulations/time[3], 0),"per second)")

return(wtp)
}

#' @title SummaryWelfare
#' @description Provide a summary of welfare changes for each policy
#' @param wtp list of welfare changes from SimulateWTP
#' @param ci confidence interval (for 95\% input 0.95)
#' @return wtp_sum summary table of welfare results
#' @export
SummaryWelfare <- function(wtp, ci = 0.95){

	wtp_sum <- apply(simplify2array(wtp),1:2, mean)
	colnames(wtp_sum)<- paste0(rep("policy",ncol(wtp_sum)), 1:ncol(wtp_sum))

	wtp_sum <- tbl_df(wtp_sum) %>%
		gather(policy, wtp) %>%
#		gather_(., "policy","wtp") %>%
#		gather_(key_col = "policy", value_col = "wtp") %>%
		group_by(policy) %>%
		summarise(mean = mean(wtp),
				  sd = sd(wtp),
				  cl_lo = quantile(wtp, (1-ci)/2),
				  cl_hi = quantile(wtp, ci+(1-ci)/2),
				  zstat = mean / sd)
	return(wtp_sum)
}

#' @title StanWTP
#' @description Use Stan functions to simulate WTP
#' @param df_indiv list of inc, quant_j, price_j, psi, and psi_p that vary by individual
#' @param df_common list of parameters that are constant for all individuals
#' @param sim_options list of simualtion options
#' @return wtp list
#' @export
StanWTP <- function(df_indiv, df_common, sim_options){

#	df_indiv <- df_indiv$df_indiv
#	df_common <- df_common$df_common

	message("Compiling simulation code")
#	expose_stan_functions(stanmodels$SimulationFunctions)

	message("Simulating...")

	if (df_common$price_change_only == FALSE){
	wtp <- pmap(df_indiv, CalcWTP_rng,
				price_p=df_common$price_p_list,
				gamma_sim=df_common$gamma_sim_list,
				alpha_sim=df_common$alpha_sim_list,
				scale_sim=df_common$scale_sim,
				nerrs=sim_options$nerrs,
				cond_error=sim_options$cond_error,
				algo_gen=sim_options$algo_gen,
				model_num=sim_options$model_num)
	} else if (df_common$price_change_only == TRUE){
		wtp <- pmap(df_indiv, CalcWTPPriceOnly_rng,
					price_p=df_common$price_p_list,
					gamma_sim=df_common$gamma_sim_list,
					alpha_sim=df_common$alpha_sim_list,
					scale_sim=df_common$scale_sim,
					nerrs=sim_options$nerrs,
					cond_error=sim_options$cond_error,
					algo_gen=sim_options$algo_gen,
					model_num=sim_options$model_num)
	}
	return(wtp)
}

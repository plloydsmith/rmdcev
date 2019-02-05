#' @title SimulateWTP
#' @description Simulate WTP for policies
#' @param stan_est Stan fit model from FitMDCEV
#' @param nsims Number of simulations to use for standard errors
#' @param nerrs Number of error draws
#' @param cond_error Choose whether to draw errors conditional on actual demand or not.
#' Conditional error draws (=1) or unconditional error draws. Default is conditional.
#' @param algo_gen Type of algorhitm for simulation. algo_gen = 0 for Hybrid Approach (i.e. constant alphas,
#' only model 3/4) alog_gen = 1 for General approach (i.e. heterogeneous alpha's, all models)
#' @return wtp list
#' @export
#'
SimulateWTP <- function(stan_est, policies,
						nerrs = 4,
						nsims = 3,
						cond_error = 1,
						algo_gen = NULL){
	nsims <- 3

	model_type <- stan_est$stan_data$model_type

	if (!is.null(algo_gen))
		if (model_type < 3 && algo_gen == 0) stop("Can't use hybrid algorithm with model_type = 1 or 2. Choose algo_gen ==1")

	if (is.null(algo_gen)) {# ensure
		if (model_type == 3 || model_type == 4)
			algo_gen <- 0
		else if (model_type == 1 || model_type == 2)
			algo_gen <- 1
	}

	sim_options <- list(ngoods = stan_est$stan_data[["J"]],
						nerrs = nerrs,
				#		nsims = nsims,
						npols = length(policies[["price_p"]]),
						cond_error = cond_error,
						algo_gen = algo_gen,
						model_type = model_type)

	# Prepare sim data
	sim_welfare <- PrepareSimulationData(stan_est, policies)

	df_common <- sim_welfare
	df_common$df_indiv <- NULL

	df_indiv <- sim_welfare$df_indiv

	wtp <- StanWTP(df_indiv, df_common, sim_options)

	return(wtp)
}


#' @title StanWTP
#' @description Use Stan functions to simulate WTP
#' @param df_indiv list of inc, quant_j, price_j, psi, and psi_p that vary by individual
#' @param df_common list of parameters that are constant for all individuals
#' @param sim_options list of simualtion options
#' @return wtp list
#' @export
#'
StanWTP <- function(df_indiv, df_common, sim_options){

	wtpcppcode <- stanc("src/stan_files/SimulationFunctions.stan",
						model_name = "SimulationFunctions")

	expose_stan_functions(wtpcppcode)

	wtp <- pmap(df_indiv, CalcWTP_rng,
				price_p=df_common$price_p_list,
				gamma_sim=df_common$gamma_sim_list,
				alpha_sim=df_common$alpha_sim_list,
				scale_sim=df_common$scale_sim,
				ngoods=sim_options$ngoods,
				nsims=sim_options$nsims,
				nerrs=sim_options$nerrs,
				npols=sim_options$npols,
				cond_error=sim_options$cond_error,
				algo_gen=sim_options$algo_gen,
				model_type=sim_options$model_type)#,
	# wtpcppcode = wtpcppcode,
	# .progress = TRUE,
	#   .options = future_options(packages=c("rstan"), globals = FALSE))
	return(wtp)
}

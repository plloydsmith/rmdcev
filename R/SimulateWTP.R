#' @title SimulateWTP
#' @description Simulate WTP for MDCEV model
#' @param stan_est Stan fit model from FitMDCEV
#' @param policies list of price_p with additive price increases and dat_psi_p with new psi data
#' @param nsims Number of simulations to use for standard errors
#' @param nerrs Number of error draws
#' @param cond_error Choose whether to draw errors conditional on actual demand or not.
#' Conditional error draws (=1) or unconditional error draws.
#' @param algo_gen Type of algorhitm for simulation. algo_gen = 0 for Hybrid Approach (i.e. constant alphas,
#' only model 3/4) alog_gen = 1 for General approach (i.e. heterogeneous alpha's, all models)
#' @return wtp a list for each individual holding a nsims x npols matrix of wtp
#' @export
#'
SimulateWTP <- function(stan_est, policies,
						nerrs = 30,
						nsims = 30,
						cond_error = 1,
						algo_gen = NULL,
						parralel = FALSE,
						n_workers = 4){

	start.time <- Sys.time()


	# Checks on simulation options
	model_type <- stan_est$stan_data$model_type

	if (!is.null(algo_gen))
		if (model_type < 3 && algo_gen == 0) stop("Can't use hybrid algorithm with model_type = 1 or 2. Choose algo_gen ==1")

	if (is.null(algo_gen)) {# ensure
		if (model_type == 3 || model_type == 4)
			algo_gen <- 0
		else if (model_type == 1 || model_type == 2)
			algo_gen <- 1
	}

	if (nsims > stan_est$n_draws) {# ensure
		nsims <- stan_est$n_draws
		warning("Number of simulations > Number of Draws from stan_est. nsims has been set to: ", nsims)
	}

	if (algo_gen == 1) # ensure
		cat("Using general approach to simulation")
	else if (algo_gen == 1)
		cat("Using hybrid approach to simulation")

	# Organize options in list
	sim_options <- list(nerrs = nerrs,
						cond_error = cond_error,
						algo_gen = algo_gen,
						model_type = model_type)

	# Prepare sim data
	sim_welfare <- PrepareSimulationData(stan_est, policies, nsims)

	df_common <- sim_welfare
	df_common$df_indiv <- NULL

	df_indiv <- sim_welfare$df_indiv


	wtp <- StanWTP(df_indiv, df_common, sim_options, parralel)

	time <- Sys.time() - start.time
	n_simulations <- length(unlist(wtp)) * nerrs
	cat(formatC(n_simulations, format = "e", digits = 2), "simulations finished in", round(time/60, 2), "minutes")

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
StanWTP <- function(df_indiv, df_common, sim_options, parralel){

	wtpcppcode <- stanc("src/stan_files/SimulationFunctions.stan",
						model_name = "SimulationFunctions")

	if (parralel == FALSE){
		expose_stan_functions(wtpcppcode)

	wtp <- pmap(df_indiv, CalcWTP_rng,
				price_p=df_common$price_p_list,
				gamma_sim=df_common$gamma_sim_list,
				alpha_sim=df_common$alpha_sim_list,
				scale_sim=df_common$scale_sim,
				nerrs=sim_options$nerrs,
				cond_error=sim_options$cond_error,
				algo_gen=sim_options$algo_gen,
				model_type=sim_options$model_type)#,

	} else if (parralel == TRUE){
		wtp <- future_pmap(df_indiv, SimulateWTPParallel,
						   price_p=df_common$price_p_list,
						   gamma_sim=df_common$gamma_sim_list,
						   alpha_sim=df_common$alpha_sim_list,
						   scale_sim=df_common$scale_sim,
						   nerrs=sim_options$nerrs,
						   cond_error=sim_options$cond_error,
						   algo_gen=sim_options$algo_gen,
						   model_type=sim_options$model_type,
						   wtpcppcode = wtpcppcode,
						   .progress = TRUE,
						   .options = future_options(packages=c("rstan"), globals = FALSE))
	}
	return(wtp)
}

SimulateWTPParallel <- function(inc, quant_j, price, price_p, psi_p_sim,
								psi_sim, gamma_sim, alpha_sim, scale_sim,
								nerrs, cond_error, algo_gen, model_type,
								wtpcppcode){
	require(rstan)
	expose_stan_functions(wtpcppcode)

	out <- CalcWTP_rng(inc, quant_j, price, price_p, psi_p_sim,
					   psi_sim, gamma_sim, alpha_sim, scale_sim,
					   nerrs, cond_error, algo_gen, model_type)

	return(out)
}


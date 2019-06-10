#' @title SimulateMDCEV
#' @description Simulate welfare or demand for MDCEV model
#' @inheritParams GenerateMDCEVData
#' @param df_indiv Prepared individual level data from PrepareSimulationData
#' @param df_common Prepared common data from PrepareSimulationData
#' @param sim_options Prepared simulation options from PrepareSimulationData
#' @param nerrs Number of error draws for welfare analysis
#' @param cond_error Choose whether to draw errors conditional on actual demand or not.
#' Conditional error draws (=1) or unconditional error draws.
#' @param draw_mlhs Generate draws using Modified Latin Hypercube Sampling algorithm (=1)
#' or uniform (=0)
#' @param algo_gen Type of algorhitm for simulation. algo_gen = 0 for Hybrid Approach (i.e. constant alphas,
#' only model 3/4) alog_gen = 1 for General approach (i.e. heterogeneous alpha's, all models)
#' @param sim_type Either "welfare" or "demand"
#' @return wtp a list for each individual holding a nsims x npols matrix of wtp
#' @export
SimulateMDCEV <- function(df_indiv, df_common, sim_options,
						sim_type = c("welfare", "demand"),
						nerrs = 30,
						cond_error = 1,
						draw_mlhs = 1,
						algo_gen = NULL,
						tol = 1e-20,
						max_loop = 999){

	start.time <- proc.time()

	# Checks on simulation options
	model_num <- sim_options$model_num

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
	sim_options[["nerrs"]] <- nerrs
	sim_options[["cond_error"]] <- cond_error
	sim_options[["draw_mlhs"]] <- draw_mlhs
	sim_options[["algo_gen"]] <- algo_gen
	sim_options[["tol"]] <- tol
	sim_options[["max_loop"]] <- max_loop

	if(sim_options$n_classes == 1){
		if(sim_type == "welfare"){
			out <- StanWelfare(df_indiv, df_common, sim_options)
		} else if(sim_type == "demand"){
			out <- StanDemand(df_indiv, df_common, sim_options)
		}

	} else if(sim_options$n_classes > 1){
		if(sim_type == "welfare"){
			out <- purrr::map2(df_indiv, df_common, StanWelfare, sim_options)
		} else if(sim_type == "demand"){
			out <- purrr::map2(df_indiv, df_common, StanDemand, sim_options)
		}
		names(out) <- paste0("class", c(1:sim_options$n_classes))
	}

	time <- proc.time() - start.time
	n_simulations <- length(unlist(out)) * nerrs
	cat("\n", formatC(n_simulations, format = "e", digits = 2), "simulations finished in", round(time[3]/60, 2), "minutes.",
		"(",round(n_simulations/time[3], 0),"per second)")

return(out)
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
		tidyr::gather(policy, wtp) %>%
		group_by(.data$policy) %>%
		summarise(mean = round(mean(wtp),2),
				  std_dev = round(stats::sd(wtp),2),
				  ci_lo = round(stats::quantile(wtp, (1-ci)/2),2),
				  ci_hi = round(stats::quantile(wtp, ci+(1-ci)/2),2))

	colnames(wtp_sum) <- c("policy", "Mean", "Std.Dev", paste0("ci_lo",(1-ci)/2*100, "%"), paste0("ci_hi",(ci+(1-ci)/2)*100, "%"))

	return(wtp_sum)
}

#' @title SummaryDemand
#' @description Provide a summary of demands each policy
#' @param demand list of welfare changes from SimulateDemand
#' @param ci confidence interval (for 95\% input 0.95)
#' @return demand_sum summary table of demand results
#' @export
SummaryDemand <- function(wtp, ci = 0.95){
	nobs <- length(demand)
	nsims <- length(demand[[1]])
	ngoods <- ncol(demand[[1]][[1]])
	npols <- nrow(demand[[1]][[1]])


	demand <- tibble(demand = unlist(demand),
					  id = rep(1:nobs, each = nsims*ngoods*npols),
					  good = rep(1:ngoods, each = npols, times = nsims*nobs),
					  policy = rep(paste0(rep("policy",npols), 1:npols), times =nobs*nsims*ngoods ),
					  sim_id = rep(1:nsims, each = ngoods*npols, times = nobs)) %>%
		group_by(.data$good, .data$policy, .data$sim_id) %>%
		summarise(demand = mean(demand)) %>%
		group_by(.data$policy, .data$good) %>%
		summarise(mean = round(mean(demand),2),
				  std_dev = round(stats::sd(demand),2),
				  ci_lo = round(stats::quantile(demand, (1-ci)/2),2),
				  ci_hi = round(stats::quantile(demand, ci+(1-ci)/2),2))

	colnames(demand) <- c("policy", "Good", "Mean", "Std.Dev", paste0("ci_lo",(1-ci)/2*100, "%"), paste0("ci_hi",(ci+(1-ci)/2)*100, "%"))

	return(demand)
}

#' @title StanWelfare
#' @description Use Stan functions to simulate Welfare
#' @param df_indiv list of income, quant_j, price, psi, and psi_p that vary by individual
#' @param df_common list of parameters that are constant for all individuals
#' @param sim_options list of simualtion options
#' @return wtp list
#' @export
StanWelfare <- function(df_indiv, df_common, sim_options){

#	df_indiv <- df_wtp$df_indiv
#	df_common <- df_wtp$df_common

	message("Compiling simulation code")
	expose_stan_functions(stanmodels$SimulationFunctions)

	message("Simulating welfare...")

	if (sim_options$price_change_only == FALSE){
		wtp <- purrr::pmap(df_indiv, CalcWTP_rng,
				price_p=df_common$price_p_list,
				gamma_sim=df_common$gamma_sim_list,
				alpha_sim=df_common$alpha_sim_list,
				scale_sim=df_common$scale_sim,
				nerrs=sim_options$nerrs,
				cond_error=sim_options$cond_error,
				draw_mlhs=sim_options$draw_mlhs,
				algo_gen=sim_options$algo_gen,
				model_num=sim_options$model_num,
				tol = sim_options$tol,
				max_loop = sim_options$max_loop)
	} else if (sim_options$price_change_only == TRUE){
		wtp <- purrr::pmap(df_indiv, CalcWTPPriceOnly_rng,
					price_p=df_common$price_p_list,
					gamma_sim=df_common$gamma_sim_list,
					alpha_sim=df_common$alpha_sim_list,
					scale_sim=df_common$scale_sim,
					nerrs=sim_options$nerrs,
					cond_error=sim_options$cond_error,
					draw_mlhs=sim_options$draw_mlhs,
					algo_gen=sim_options$algo_gen,
					model_num=sim_options$model_num,
					tol = sim_options$tol,
					max_loop = sim_options$max_loop)
	}
	return(wtp)
}


#' @title StanDemand
#' @description Use Stan functions to simulate Marshallian demand
#' @param df_indiv list of income, quant_j, price, psi, and psi_p that vary by individual
#' @param df_common list of parameters that are constant for all individuals
#' @param sim_options list of simualtion options
#' @return demand with nsim lists of npolsXngoods+1 matrices
#' @export
StanDemand <- function(df_indiv, df_common, sim_options){

	message("Compiling simulation code")
	expose_stan_functions(stanmodels$SimulationFunctions)

	message("Simulating demand...")

	if (sim_options$price_change_only == FALSE){
		demand <- purrr::pmap(df_indiv, CalcMarshallianDemand_rng,
						   price_p=df_common$price_p_list,
						   gamma_sim=df_common$gamma_sim_list,
						   alpha_sim=df_common$alpha_sim_list,
						   scale_sim=df_common$scale_sim,
						   nerrs=sim_options$nerrs,
						   cond_error=sim_options$cond_error,
						   draw_mlhs=sim_options$draw_mlhs,
						   algo_gen=sim_options$algo_gen,
						   model_num=sim_options$model_num,
						   tol = sim_options$tol,
						   max_loop = sim_options$max_loop)
	} else if (sim_options$price_change_only == TRUE){
		demand <- purrr::pmap(df_indiv, CalcMarshallianDemandPriceOnly_rng,
						   price_p=df_common$price_p_list,
						   gamma_sim=df_common$gamma_sim_list,
						   alpha_sim=df_common$alpha_sim_list,
						   scale_sim=df_common$scale_sim,
						   nerrs=sim_options$nerrs,
						   cond_error=sim_options$cond_error,
						   draw_mlhs=sim_options$draw_mlhs,
						   algo_gen=sim_options$algo_gen,
						   model_num=sim_options$model_num,
						   tol = sim_options$tol,
						   max_loop = sim_options$max_loop)
	}
	return(demand)
}

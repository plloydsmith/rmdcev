#' @title mdcev.sim
#' @description Simulate welfare or demand for MDCEV model
#' @inheritParams GenerateMDCEVData
#' @param df_indiv Prepared individual level data from PrepareSimulationData
#' @param df_common Prepared common data from PrepareSimulationData
#' @param sim_options Prepared simulation options from PrepareSimulationData
#' @param nerrs Number of error draws for welfare analysis
#' @param cond_error Choose whether to draw errors conditional on actual demand or not.
#' Conditional error draws (=1) or unconditional error draws.
#' @param draw_mlhs Generate draws using Modified Latin Hypercube Sampling
#' algorithm (=1) or uniform (=0)
#' @param algo_gen Type of algorhitm for simulation. algo_gen = 0 for Hybrid Approach
#' (i.e. constant alphas, only model 3/4) alog_gen = 1 for General approach
#' (i.e. heterogeneous alpha's, all models)
#' @param sim_type Either "welfare" or "demand"
#' @param suppressTime Supress simulation time calculation
#' @param ... Additional parameters to pass to mdcev.sim
#' @return a object of class mdcev.sim which contains a list for each
#' individual holding either 1) nsims x npols matrix of welfare changes if
#' welfare is being simulated or 2) nsims number of lists of npols x # alternatives
#' matrix of Marshallian demands is demand is being simulated.
#' @export
#' @seealso  [mdcev()] for the estimation of mdcev models.
#' @examples
#' \donttest{
#' data(data_rec, package = "rmdcev")
#'
#' data_rec <- mdcev.data(data_rec, subset = id < 500,
#'                 alt.var = "alt", choice = "quant")
#'
#' mdcev_est <- mdcev( ~ 1, data = data_rec,
#'                model = "hybrid0", algorithm = "MLE")
#'
#' policies <- CreateBlankPolicies(npols = 2,
#' ngoods = mdcev_est[["stan_data"]][["J"]],
#' dat_psi = mdcev_est[["stan_data"]][["dat_psi"]],
#' price_change_only = TRUE)
#'
#' df_sim <- PrepareSimulationData(mdcev_est, policies)
#'
#' wtp <- mdcev.sim(df_sim$df_indiv,
#' df_common = df_sim$df_common,
#' sim_options = df_sim$sim_options,
#' cond_err = 1, nerrs = 5, sim_type = "welfare")
#'}
mdcev.sim <- function(df_indiv, df_common, sim_options,
					  sim_type = c("welfare", "demand"),
					  nerrs = 30,
					  cond_error = 1,
					  draw_mlhs = 1,
					  algo_gen = NULL,
					  tol = 1e-20,
					  max_loop = 999,
					  suppressTime = FALSE,
					  ...){

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

	if(sim_type == "welfare"){
		out <- StanWelfare(df_indiv, df_common, sim_options)
	} else if(sim_type == "demand"){
		out <- StanDemand(df_indiv, df_common, sim_options)
	}

	if(suppressTime == FALSE){
		time <- proc.time() - start.time
		n_simulations <- length(unlist(out)) * nerrs
		message("\n", formatC(n_simulations, format = "e", digits = 2),
				"simulations finished in", round(time[3]/60, 2), "minutes.",
				"(",round(n_simulations/time[3], 0),"per second)")
	}

	out <- structure(out,
					 class = "mdcev.sim")

	return(out)
}

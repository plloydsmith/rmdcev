PrepareSimulationData <- function(result, policies){
	###########################################################################
	# Create policy scenarios (can affect price only at this point)
	###########################################################################

	price_p_list <- policies$price_p

	dat_psi_list <- policies$dat_psi

	###########################################################################
	# Get parameter estimates in matrix form
	###########################################################################
	est_pars <- tbl_df(result[["stan_fit"]][["theta_tilde"]]) %>%
		select(-starts_with("log_like"), -starts_with("sum_log_lik")) %>%
		rowid_to_column("sim_id") %>%
		gather(parms, value, -sim_id, factor_key=TRUE)

	# Sample from draws
	est_sim <- result$est_pars %>%
		distinct(sim_id) %>%
		sample_n(., nsims ) %>%
		left_join(est_pars, by = "sim_id")

	#est_sim <- est_sim %>%
	#	filter(sim_id != 30)

	#nsims = nsims - 1
	parms_sim <- CreateSimulationData(est_sim, nsims, npols, result$stan_data, dat_psi_list)

	# Put in a list for each simulation
	gamma_sim_list <- CreateListsRow(parms_sim[["gamma_sim"]]) # ngoods length

	alpha_sim_list <- CreateListsRow(parms_sim[["alpha_sim"]]) # ngoods+1 length
	#alpha_sim <- matrix(0.000001, nsims, ngoods+1)
	scale_sim <- parms_sim[["scale_sim"]]

	psi_sim <- list(parms_sim[["psi_sim"]])
	names(psi_sim) <- "psi_sim"

	psi_p_sim <- list(parms_sim[["psi_p_sim"]])
	names(psi_p_sim) <- "psi_p_sim"

	###########################################################################
	# Set baseline data into lists
	###########################################################################

	inc <- list(as.list(result$stan_data$inc))
	names(inc) <- "inc"

	quant_j <- list(CreateListsRow(result$stan_data$j_quant))
	names(quant_j) <- "quant_j"

	price <- cbind(1, result$stan_data$j_price) #add numeraire price to price matrix (<-1)
	price <- list(CreateListsRow(price))
	names(price) <- "price"

	###########################################################################
	# Pull individual level data into one list
	###########################################################################

	df_indiv <- c(inc, quant_j, price, psi_sim, psi_p_sim)

	out <- list(df_indiv = df_indiv,
				price_p_list = price_p_list,
				gamma_sim_list = gamma_sim_list,
				alpha_sim_list = alpha_sim_list,
				scale_sim = scale_sim)

	return(out)
}



SimulateWTP <- function(df_indiv, df_common, sim_options,
						wtpcppcode){

	## Check models
	model <- sim_options$model_type
	algo_gen <- sim_options$algo_gen
	if (model < 3 && algo_gen == 0) stop("Can't use hybrid algorithm with model_type = 1 or 2. Choose algo_gen ==1")

wtpcppcode <- stanc("src/stan_files/SimulationFunctions.stan",
				model_name = "SimulationFunctions")

expose_stan_functions(wtpcppcode)

	wtp <- pmap(df_indiv, CalcWTP_rng,
					   price_p=df_common$price_p_list,
					   gamma_sim=df_common$gamma_sim_list,
					   alpha_sim=df_common$alpha_sim_list,
					   scale_sim=df_common$scale_sim,
					   ngoods=sim_options$ngoods, nsims=sim_options$nsims,
					   nerrs=sim_options$nerrs, npols=sim_options$npols,
					   cond_error=sim_options$cond_error, algo_gen=sim_options$algo_gen,
					   model_type=sim_options$model_type)#,
					  # wtpcppcode = wtpcppcode,
					  # .progress = TRUE,
					#   .options = future_options(packages=c("rstan"), globals = FALSE))
	return(wtp)
}


SimulateWTPParallel <- function(inc, quant_j, price, price_p, psi_p_sim,
						psi_j_sim, gamma_j_sim, alpha_sim, scale_sim,
						ngoods, nerrs, nsims, npols,
						cond_error, algo_gen, model_type,
						wtpcppcode){
	require(rstan)
	expose_stan_functions(wtpcppcode)

	out <- CalcWTP_rng(inc, quant_j, price, price_p, psi_p_sim,
					   psi_j_sim, gamma_j_sim, alpha_sim, scale_sim,
					   ngoods, nerrs, npols, nsims,
					   cond_error, algo_gen, model_type)

	return(out)
}

StartWTPParallel <- function(df_input){

	wtp <- future_pmap(df_input, SimulateWTP,
					   price_p=price_p_list,
					   gamma_j_sim=gamma_j_sim_list,
					   alpha_sim=alpha_sim_list,
					   scale_sim=scale_sim,
					   ngoods=ngoods, nerrs=nerrs, nsims=nsims, npols=npols,
					   cond_error=cond_error, algo_gen=algo_gen, model_type=model_type,
					   wtpcppcode = wtpcppcode,
					   .progress = TRUE,
					   .options = future_options(packages=c("rstan"), globals = FALSE))
	return(wtp)
}



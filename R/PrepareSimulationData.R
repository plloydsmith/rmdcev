PrepareSimulationData <- function(stan_est, policies, nsims){

	price_p_list <- policies$price_p

	###########################################################################
	# Get parameter estimates in matrix form
	###########################################################################
	est_pars <- tbl_df(stan_est[["stan_fit"]][["theta_tilde"]]) %>%
		select(-starts_with("log_like"), -starts_with("sum_log_lik")) %>%
		rowid_to_column("sim_id") %>%
		gather(parms, value, -sim_id, factor_key=TRUE)

	# Sample from draws
	est_sim <- stan_est$est_pars %>%
		distinct(sim_id) %>%
		sample_n(., nsims ) %>%
		left_join(est_pars, by = "sim_id")

	#est_sim <- est_sim %>%
	#	filter(sim_id != 30)

	#nsims = nsims - 1
	parms_sim <- CreateSimulationData(est_sim, nsims, npols, stan_est$stan_data, policies$dat_psi)

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

	inc <- list(as.list(stan_est$stan_data$inc))
	names(inc) <- "inc"

	quant_j <- list(CreateListsRow(stan_est$stan_data$j_quant))
	names(quant_j) <- "quant_j"

	price <- cbind(1, stan_est$stan_data$j_price) #add numeraire price to price matrix (<-1)
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

CreateSimulationData <- function(est_sim, nsims, data, data_psi_p){
#data <- result$stan_data
#data_psi_p <- dat_psi_list
gamma_sim <- est_sim %>%
	filter(str_detect(parms, "gamma")) %>%
	spread(sim_id, value) %>%
	select(-parms) %>%
	as.matrix(.) %>%
	t(.)

psi_temp <- est_sim %>%
	filter(str_detect(parms, "psi")) %>%
	spread(sim_id, value) %>%
	select(-parms) %>%
	as.matrix(.)

psi_temp <- lapply(seq_len(ncol(psi_temp)), function(i) psi_temp[,i])

psi_sim <- map(psi_temp, function(x){matrix(data$dat_psi %*% x , nrow = data$I, byrow = TRUE)})

psi_sim <- do.call(cbind, psi_sim)
psi_sim <- lapply(seq_len(nrow(psi_sim)), function(i) psi_sim[i,])
psi_sim <- map(psi_sim, function(x){matrix(x , nrow = nsims, byrow = TRUE)})

#psi_p_sim <- map(data_psi_p, function(data_p){ map(psi_temp, function(x){matrix(data_p %*% x , nrow = data$I, byrow = TRUE)})})
psi_p_sim <- map(psi_temp, function(psi){ map(data_psi_p, function(x){matrix(x %*% psi , nrow = data$I, byrow = TRUE)})})

psi_p_sim <- map(psi_p_sim, function(x){do.call(cbind, x)})
psi_p_sim <- do.call(cbind, psi_p_sim)
psi_p_sim <- lapply(seq_len(nrow(psi_p_sim)), function(i) psi_p_sim[i,])
psi_p_sim <- map(psi_p_sim, function(x){aperm(array(x , dim = c(data$J, npols, nsims)), perm=c(2,1,3))})
psi_p_sim <- map(psi_p_sim, function(x){lapply(seq_len(nsims), function(i) x[,,i])})

#psi_p_sim <- map(psi_temp, function(x){matrix(data$dat_psi %*% x , nrow = data$I, byrow = TRUE)})

#psi_p_sim <- do.call(cbind, psi_p_sim)

#psi_p_sim <- map(psi_p_sim, function(x){lapply(seq_len(data$I), function(i) x[i,])})
#psi_p_sim <- map(psi_p_sim, function(data_p){map(data_p, function(x){matrix(x , nrow = nsims, byrow = TRUE)})})

if (data$model_type != 4){
	alpha_sim <- est_sim %>%
		filter(str_detect(parms, "alpha")) %>%
		spread(sim_id, value) %>%
		select(-parms) %>%
		as.matrix(.) %>%
		t(.)
	if (data$model_type == 1){
		alpha_sim <- cbind(alpha_sim, matrix(0, nsims, data$J) )
	} else if (data$model_type == 3)
		alpha_sim <- matrix(rep(alpha_sim,each=data$J+1), ncol=data$J+1, byrow=TRUE)

} else if (data$model_type ==4)
	alpha_sim = matrix(0, nsims, data$J+1)

if (data$fixed_scale == 0){
	scale_sim <- est_sim %>%
		filter(str_detect(parms, "scale")) %>%
		spread(sim_id, value) %>%
		select(-parms) %>%
		as.matrix(.) %>%
		t(.)
} else if (data$fixed_scale == 1){
	scale_sim = matrix(1, nsims, 1)
}

return(list(psi_sim=psi_sim, psi_p_sim = psi_p_sim, gamma_sim=gamma_sim, alpha_sim=alpha_sim, scale_sim=scale_sim))
}

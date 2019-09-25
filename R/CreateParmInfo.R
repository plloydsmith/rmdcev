#' @title CreateParmInfo
#' @description Create parameter names and number of parameters
#' @param stan_data data for model
#' @inheritParams mdcev
#' @return A list of parameter names and numbers
#' @keywords internal
CreateParmInfo <- function(stan_data, algorithm, random_parameters){

J <- stan_data$J

n_psi <- stan_data$NPsi
psi_names <- paste0(rep('psi', n_psi), sep="_",
					colnames(stan_data[["dat_psi"]]))


if (stan_data$model_num == 1){
	n_alpha <- 1
	n_gamma <- J
} else if (stan_data$model_num == 2){
	n_alpha <- J + 1
	n_gamma <- 0
} else if (stan_data$model_num == 3){
	n_alpha <- 1
	n_gamma <- J
} else if (stan_data$model_num == 4){
	n_alpha <- 0
	n_gamma <- J
}

if (stan_data$fixed_scale1 == 1){
	n_scale <- 0
	scale_names <- NULL
} else {
	n_scale <- 1
	scale_names <- "scale"
}


if (stan_data$model == 2){
	gamma_names <- NULL
} else if (stan_data$model != 2){
	gamma_names <- paste0(rep('gamma', n_gamma), sep = "",
						  c(1:n_gamma))
}
if (stan_data$model == 4){
	alpha_names <- NULL
} else if (stan_data$model != 4){
	alpha_names <- paste0(rep('alpha', n_alpha), sep = "",
						  c(1:n_alpha))
}

parm_names <- list(psi_names=psi_names,
				   gamma_names=gamma_names,
				   alpha_names=alpha_names,
				   scale_names=scale_names,
				   all_names = c(psi_names, gamma_names, alpha_names, scale_names))

n_parameters <- n_psi + n_alpha + n_gamma + n_scale
n_vars = list(n_psi = n_psi, n_alpha = n_alpha, n_gamma = n_gamma, n_scale = n_scale)

if (stan_data$K > 1){
	n_vars <- purrr::map(n_vars, function(x){x* stan_data$K})
	n_vars$n_beta <- stan_data$L * (stan_data$K - 1)

	all_names <- GenClassNames(parm_names$all_names, stan_data$K)

	delta.names <- GenClassNames(colnames(stan_data[["data_class"]]), stan_data$K)
	parm_names$delta.names <-	grep("class1", delta.names, invert=TRUE, value = TRUE)

	parm_names$all_names <- c(all_names, parm_names$delta.names)
}
#standard deviations
if(random_parameters != "fixed"){
	if(stan_data$gamma_fixed == 1){
		n_gamma_rp <- 0
		gamma_sd_names <- NULL
	} else if(stan_data$gamma_fixed == 0){
		n_gamma_rp <- n_gamma
		gamma_sd_names <- gamma_names
	}

	if(stan_data$alpha_fixed == 1){
		n_alpha_rp <- 0
		alpha_sd_names <- NULL
	} else if(stan_data$alpha_fixed == 0){
		n_alpha_rp <- n_alpha
		alpha_sd_names <- alpha_names
	}
	n_vars$n_std_dev <- n_psi + n_gamma_rp + n_alpha_rp
	parm_names$sd_names <- paste0("sd.", c(psi_names, gamma_sd_names, alpha_sd_names))

	if (random_parameters == "corr"){
		n_vars$n_std_corr <- n_vars$n_std_dev * (n_vars$n_std_dev-1) / 2
	}
}

n_vars$n_parms_total <- Reduce("+",n_vars)

parms_info <- list(n_vars = n_vars,
				   parm_names = parm_names)

return(parms_info)
}


#' @title GenClassNames
#' @description Create names for LC
#' @param parms_names list of parameter names
#' @inheritParams FitMDCEV
#' @return A vector of LC names
#' @export
GenClassNames <- function(parms_names, n_classes){
	if(length(parms_names) > 0){
		classes <- paste0(rep('class', n_classes), sep = "",c(1:n_classes))
		parms_names <- paste0(rep(classes, length(parms_names)), sep=".", rep(parms_names,  each = n_classes))
	} else if (length(parms_names) == 0){
		parms_names <- NULL
	}
	return(parms_names)
}

#' @title CreateListsRow
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @return A list
#' @export
#' @examples
#'
#' tmp <- matrix(0, nrow = 10, ncol = 5)
#' tmp_list <- CreateListsRow(tmp)
#'
CreateListsRow <- function(x){
	out <- lapply(seq_len(nrow(x)), function(i) x[i,])
	return(out)
}

#' @title CreateListsCol
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @export
#' @examples
#'
#' tmp <- matrix(0, nrow = 10, ncol = 5)
#' tmp_list <- CreateListsCol(tmp)
#'
CreateListsCol <- function(x){
	if (is.vector(x))
		out <- as.list(x)
	else
		out <- lapply(seq_len(ncol(x)), function(i) x[,i])
	return(out)
}

#' @title MultiplyMatrix
#' @description Convert list to matrix by row
#' @param x matrix B to be multiplied
#' @param mat_temp matrix A to be multipled
#' @param n_rows number of rows for final matrix
#' @keywords internal
MultiplyMatrix <- function(mat_temp, x, n_rows){
	out <- matrix(mat_temp %*% x , nrow = n_rows, byrow = TRUE)
	return(out)
}

#' @title DoCbind
#' @description Convert list to matrix
#' @param x list to be converted
#' @keywords internal
DoCbind <- function(x){
	out <- do.call(cbind, x)
	return(out)
}

#' @title CreateBlankPolicies
#' @description Create 'zero effect' policies that can be modified
#' @param npols Number of policies to simulate
#' @param model Estimated model from mdcev
#' @param price_change_only Logical value for whether to include policy changes to dat_psi. Defaults to TRUE.
#' TRUE implies that only price changes are used in simulation.
#' @export
#' @examples
#' \donttest{
#' data_rec <- mdcev.data(data_rec, subset = id <= 500, id.var = "id",
#'                 alt.var = "alt", choice = "quant")
#'
#' mdcev_est <- mdcev( ~ 0, data = data_rec,
#'                model = "hybrid0", algorithm = "MLE",
#'                std_errors = "mvn")
#' CreateBlankPolicies(npols = 2, mdcev_est)
#'}
CreateBlankPolicies <- function(npols, model, price_change_only = TRUE){

	price_p <- CreateListsRow(matrix(0, nrow = npols, ncol = model[["stan_data"]][["J"]] + 1))
	model_num <- model[["stan_data"]][["model_num"]]

	dat_psi_p <- NULL
	dat_phi_p <- NULL

	if (price_change_only == FALSE){
		if (model_num < 5 && model[["parms_info"]][["n_vars"]][["n_psi"]] == 0)
			stop("No psi variables to vary! Use price_change_only == TRUE option.")

		if (model_num == 5 && model[["parms_info"]][["n_vars"]][["n_phi"]] == 0)
			stop("No phi variables to vary! Use price_change_only == TRUE option.")

		if (model_num < 5)
			dat_psi_p <- lapply(seq_len(npols), function(X) model[["stan_data"]][["dat_psi"]])
		else if (model_num == 5)
			dat_phi_p <- lapply(seq_len(npols), function(X) model[["stan_data"]][["dat_phi"]])
	}

	out <- list(price_p = price_p, dat_psi_p = dat_psi_p, dat_phi_p = dat_phi_p,
				price_change_only = price_change_only)
	return(out)
}

#' @title CreatePsiMatrix
#' @param psi_j matrix (JXn_psi_j) of alt-specific attributes
#' @param psi_i matrix (IXn_psi_i) of people-specific attributes
#' @description Creates the Psi data matrix for use in mdcev model
#' @keywords internal
CreatePsiMatrix <- function(psi_j = NULL, psi_i = NULL){
	if(!is.na(psi_i))
		psi_i <- lapply(psi_i, function(x) {rep(x, each= nrow(psi_j))})
	if(!is.na(psi_j))
		psi_j <- lapply(psi_j, function(x) {rep(x, times=nrow(psi_i))})

	dat_psi <- c(psi_j, psi_i)
return(dat_psi)
}

#' @title GrabParms
#' @param data est_par object from results
#' @param parm_name name of parameter to get simulations
#' @description Pulls out specific mdcev parameter simulations
#' @keywords internal
GrabParms <- function(data, parm_name){
	out <- data %>%
		dplyr::filter(stringr::str_detect(.data$parms, parm_name)) %>%
		tidyr::pivot_wider(sim_id, names_from = "parms", values_from = "value") %>%
		dplyr::select(-sim_id) %>%
		as.matrix(.)
	return(out)
}


#' @title GrabIndividualParms
#' @param est_sim est_sim from results
#' @param parm_name name of parameter to get simulations
#' @description Pulls out specific mdcev parameter simulations
#' @keywords internal
GrabIndividualParms <- function(est_sim, parm_name){
	out <- est_sim %>%
		dplyr::filter(grepl(c(parm_name), parms)) %>%
		dplyr::select(id, sim_id, .data$parm_id, beta) %>%
		tidyr::spread(.data$parm_id, beta) %>%
		dplyr::select(-sim_id) %>%
		dplyr::group_split(id, .keep = F)
	return(out)
}

#' @title CreatePsi
#' @param dat_vars_i psi data for each person
#' @param est_pars_i psi parameter estimates for each person
#' @description Works at individual level and creates the psi
#' variables for each simulation, policy, alternative
#' @keywords internal
CreatePsi = function(dat_vars_i, est_pars_i, J, NPsi_ij, psi_ascs, npols){
#	dat_vars_i = dat_vars[[3]]
#	est_pars_i = psi_sim_temp[[3]]
	lpsi = matrix(0, nrow(est_pars_i), J)
	if (psi_ascs == 1){
		psi_non_ascs_start = J
		psi_non_ascs_end = J+NPsi_ij-1
		if (nrow(est_pars_i) == 1)
			lpsi = lpsi + c(0, est_pars_i[,1:(J-1)])
		else
			lpsi = lpsi + cbind(0, est_pars_i[,1:(J-1)]) ##  alternative specific constants
	} else if (psi_ascs == 0){
		psi_non_ascs_start = 1
		psi_non_ascs_end = NPsi_ij
	}

	if ((NPsi_ij > 0) && (nrow(dat_vars_i) == J)){
		psi_non_ascs = est_pars_i[,psi_non_ascs_start:psi_non_ascs_end, drop=FALSE]
		lpsi = lpsi + as.matrix(psi_non_ascs) %*% t(as.matrix(dat_vars_i))
	}

	if (nrow(dat_vars_i) > J){
		if (NPsi_ij == 0){
			lpsi = CreateListsRow(lpsi)
			lpsi = lapply(lpsi, function(x){
				lpsi= matrix(x,nrow=npols,ncol=length(x),byrow=TRUE)
			})
		} else if (NPsi_ij > 0){
			dat_vars_i = dat_vars_i %>%
				group_split(policy, .keep = F)
			psi_non_ascs <- est_pars_i[,psi_non_ascs_start:psi_non_ascs_end, drop=FALSE]
			lpsi <- lapply(dat_vars_i, function(xx){
				lpsi = lpsi + as.matrix(psi_non_ascs) %*% t(as.matrix(xx))
				return(lpsi)
			})
			lpsi = aperm(array(unlist(lpsi),
							   dim = c(nrow(lpsi[[1]]), ncol(lpsi[[1]]), length(lpsi))),
						 perm=c(1,3,2))
			lpsi = lapply(seq(dim(lpsi)[1]), function(x) lpsi[ x, , ])
		}
	}
	return(lpsi)
}

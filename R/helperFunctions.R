#' @title CreateListsRow
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @return A list
#' @export
#'
CreateListsRow <- function(x){
	out <- lapply(seq_len(nrow(x)), function(i) x[i,])
	return(out)
}

#' @title CreateListsCol
#' @description Convert matrix x to a list with each row as an element
#' @param x matrix to be converted to list
#' @export
CreateListsCol <- function(x){
	out <- lapply(seq_len(ncol(x)), function(i) x[,i])
	return(out)
}

#' @title MultiplyMatrix
#' @description Convert list to matrix by row
#' @param x matrix B to be multiplied
#' @param mat_temp matrix A to be multipled
#' @param z number of rows for final matrix
#' @export
MultiplyMatrix <- function(mat_temp, x, n_rows){
	out <- matrix(mat_temp %*% x , nrow = n_rows, byrow = TRUE)
	return(out)
}

#' @title DoCbind
#' @description Convert list to matrix
#' @param x list to be converted
#' @export
DoCbind <- function(x){
	out <- do.call(cbind, x)
	return(out)
}

#' @title CreateBlankPolicies
#' @description Create 'zero effect' policies that can be modified
#' @param npols Number of policies to simulate
#' @param ngoods Number of non-numeraire goods
#' @param dat_psi Psi data matrix used in estimation
#' @export
CreateBlankPolicies <- function(npols, ngoods, dat_psi){
	price_p <- CreateListsRow(matrix(0, nrow = npols, ncol = ngoods + 1))
	dat_psi_p <- lapply(seq_len(npols), function(X) dat_psi)
	out <- list(price_p = price_p, dat_psi_p = dat_psi_p)
	return(out)
}

#' @title CreatePsiMatrix
#' @param psi_j matrix (JXn_psi_j) of good-specific attributes
#' @param psi_i matrix (IXn_psi_i) of people-specific attributes
#' @description Creates the Psi data matrix for use in mdcev model
#' @export
CreatePsiMatrix <- function(psi_j = NULL, psi_i = NULL){
	if(!is.na(psi_i))
		psi_i <- map(psi_i, function(x) {rep(x, each= nrow(psi_j))})
	if(!is.na(psi_j))
		psi_j <- map(psi_j, function(x) {rep(x, times=nrow(psi_i))})

	dat_psi <- c(psi_j, psi_i)
return(dat_psi)
}

#' @title GrabParms
#' @param data est_par object from results
#' @param parm_name name of parameter to get simulations
#' @description Pulls out specific mdcev parameter simulations
#' @export
GrabParms <- function(data, parm_name){
	out <- data %>%
		filter(str_detect(.data$parms, parm_name)) %>%
		spread_(key_col = 'sim_id',
				value_col = 'value') %>%
		select_(.dots = c('-parms')) %>%
		as.matrix(.)
	return(out)
}



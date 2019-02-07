#' @title CreateListsRow
#' @description Convert matrix x to a list with each row as an element
#' @return A list
#' @export
#'
CreateListsRow <- function(x){
	out <- lapply(seq_len(nrow(x)), function(i) x[i,])
	return(out)
}

#' @title CreateListsCol
#' @description Convert matrix x to a list with each row as an element
#' @export
CreateListsCol <- function(x){
	out <- lapply(seq_len(ncol(x)), function(i) x[,i])
	return(out)
}

#' @title MultiplyMatrix
#' @description Convert list to matrix by row
#' @export
MultiplyMatrix <- function(mat_temp, x, nrows){
	out <- matrix(mat_temp %*% x , nrow = nrows, byrow = TRUE)
	return(out)
}

#' @title DoCbind
#' @description Convert list to matrix
#' @export
DoCbind <- function(x){
	out <- do.call(cbind, x)
	return(out)
}

#' @title CreateBlankPolicies
#' @description Create 'zero effect' policies that can be modified
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
out(dat_psi)
}



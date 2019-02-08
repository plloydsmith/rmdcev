#' @title CheckMdcevData
#' @description Check MDCEV data
#' @inheritParams FitMDCEV
#' @export
CheckMdcevData <- function(data){

	message("Checking data...")
	I <- nrow(data$price)
	J <- ncol(data$price)

	if (identical(dim(data$price), dim(data$quant)) == FALSE)
		stop("price and quant dimension mismatch. Ensure dim(price) = dim(quant)")

	if (identical(nrow(data$dat_psi), I*J )  == FALSE)
		stop("dat_psi not I x J")

	if (identical(length(data$inc), I )  == FALSE)
		stop("inc data not length I")

	message("Data is good")
}

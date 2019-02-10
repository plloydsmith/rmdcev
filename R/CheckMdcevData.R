#' @title CheckMdcevData
#' @description Check MDCEV data
#' @inheritParams FitMDCEV
#' @export
CheckMdcevData <- function(data){

	message("Checking data...")

	if(!"id" %in% colnames(data))
		stop("Data must have id column for individual")

	if(!"good" %in% colnames(data))
		stop("Data must have good column for non-numeraire alternatives")

	if(!"quant" %in% colnames(data))
		stop("Data must have quant column for consumption")

	if(!"price" %in% colnames(data))
		stop("Data must have price column for non-numeraire alternatives")

	if(!"inc" %in% colnames(data))
		stop("Data must have inc column for individual's income")

	check <- tbl_df(data) %>%
		mutate(expend = price * quant) %>%
		group_by(id) %>%
		summarise(numeraire = mean(inc) - sum(expend)) %>%
		select(numeraire)

	if (sum(check$numeraire < 0) > 0)
		stop("Numeraire is less than 0 for at least one individual")

	if (sum(data$price <= 0) > 0)
		stop("Price is less than or equal to 0 for at least one individual good")

	message("Data is good")
}

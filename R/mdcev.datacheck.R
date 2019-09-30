#' @title mdcev.datacheck
#' @description Check mdcev data
#' @inheritParams mdcev
#' @keywords internal
mdcev.datacheck <- function(data_input){

	message("Checking data...")

	if(!"id" %in% colnames(data_input))
		stop("Data must have id column for individual")

	if(!"quant" %in% colnames(data_input))
		stop("Data must have quant column for consumption")

	if(!"price" %in% colnames(data_input))
		stop("Data must have price column for non-numeraire alternatives")

	if(!"income" %in% colnames(data_input))
		stop("Data must have income column for individual's income")

	check <- tbl_df(data_input) %>%
		mutate(expend = .data$price * .data$quant) %>%
		group_by(.data$id) %>%
		summarise(numeraire = mean(.data$income) - sum(.data$expend)) %>%
		dplyr::select(.data$numeraire)

	if (sum(check$numeraire < 0) > 0){
		print(paste0("Numeraire is less than 0 for individuals in rows: ",
					 toString(as.character(which(check$numeraire < 0)))))
		stop()
	}

	if (sum(data_input$price <= 0) > 0){
		print(paste0("Price is less than or equal to 0 for at least one individual alt in rows: ",
					 toString(as.character(which(data_input$price <= 0)))))
		stop()
	}
	message("Data is good")
}

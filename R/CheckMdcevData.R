#' @title CheckMdcevData
#' @description Check MDCEV data
#' @inheritParams FitMDCEV
#' @keywords internal
CheckMdcevData <- function(data_input){

	message("Checking data...")

	if(!"id" %in% colnames(data_input))
		stop("Data must have id column for individual")

	if(!"good" %in% colnames(data_input))
		stop("Data must have good column for non-numeraire alternatives")

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

	check <- tbl_df(data_input) %>%
		arrange(.data$id, .data$good)

	if (identical(identical(check$id, data_input$id), FALSE) || identical(identical(check$good, data_input$good), FALSE) ){
		stop("Data must be arranged by id then good. Can use dplyr command arrange(id, good).")
	}

	if (sum(data_input$price <= 0) > 0){
		print(paste0("Price is less than or equal to 0 for at least one individual good in rows: ",
					 toString(as.character(which(data_input$price <= 0)))))
		stop()
	}
	message("Data is good")
}

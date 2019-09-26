#' data.frame for mdcev model
#'
#' shape a `data.frame` in a suitable form for the use of the
#' `mdcev` function and complete some data checks
#'
#' @name mdcev.data
#' @param data a `data.frame`,
#' @param choice the variable indicating the choice made: it has to be
#'     a numerical vector
#' @param alt.var the name of the variable that contains the
#'     alternative index or the name under which the alternative index
#'     will be stored (the default  name is `alt`),
#' @param alt.levels the name of the alternatives: if null,
#'     they are guessed from the `alt.var` argument,
#' @param id.var the name of the variable that contains the individual
#'     index if any,
#' @param chid.var the name of the variable that contains the choice
#'     index or the name under which the choice index will be stored
#' @param drop.index should the index variables be dropped from the
#'     `data.frame`,
#' @param group.var the name of the variable that contains the group
#'     index if any,
#' @param subset a logical expression which defines the subset of
#'     observations to be selected,
#' @param ... further arguments.
#' @return A `mdcev.data` object, which is a `data.frame` in long
#'     format, *i.e.* one line for each alternative.  It has a `index`
#'     attribute, which is a `data.frame` that contains the index of
#'     the choice made (`chid`), the index of the alternative (`alt`)
#'     and, if any, the index of the individual (`id`) and of the
#'     alternative groups (`group`).
#' @export
mdcev.data <- function(data, choice = NULL, alt.var = NULL, chid.var = NULL,
						alt.levels = NULL, id.var = NULL, group.var = NULL,
						drop.index = FALSE,
						subset = NULL, ...){
	# chid.name, alt.name : the name of the index variables
	# chid, alt : the index variables

	# if a subset argument is prodided, subset the original data frame
	cldata <- match.call(expand.dots = TRUE)
	if (match("subset", names(cldata), 0)){
		m <- match(c("data", "subset"), names(cldata), 0)
		cldata <- cldata[c(1, m)]
		names(cldata)[2] <- "x"
		cldata[[1]] <- as.name("subset")
		data <- eval(cldata, parent.frame())
	}
	if (is.null(chid.var)){
		chid.name <- "chid"
		chid.is.variable <- FALSE
	}
	else{
		chid.name <- chid.var
		chid.is.variable <- ifelse(is.null(data[[chid.var]]), FALSE, TRUE)
	}
	if (is.null(alt.var) && is.null(alt.levels))
		stop("at least one of alt.var and alt.levels should be filled")

	if (! is.null(alt.levels)){
		J <- length(alt.levels)
		n <- nrow(data) / J
		alt <- factor(rep(alt.levels, n), levels = alt.levels)
		if (!is.null(alt.var) && !is.null(data[[alt.var]])){
			warning(paste("variable", alt.var, "exists and will be replaced"))
			alt.is.variable <- TRUE
		}
		else alt.is.variable <- FALSE
		alt.name <- ifelse(is.null(alt.var), "alt", alt.var)
	}
	else{
		alt.name <- alt.var
		alt.is.variable <- TRUE
		if (!is.factor(data[[alt.name]])) data[[alt.name]] <- factor(data[[alt.name]])
		alt.levels <- levels(data[[alt.name]])
		J <- length(alt.levels)
		alt <- data[[alt.name]]

	}
	n <- nrow(data) / J
	if (! chid.is.variable) chid <- rep(1:n, each = J) else chid <- data[[chid.name]]
	if (! is.null(choice)){
		choice.name <- choice
		choice <- data[[choice]]
		data[[choice.name]] <- choice
	}
	chid <- as.factor(chid)
	alt <- as.factor(alt)
	row.names(data) <- paste(chid, alt, sep = ".")

	chidpos <- which(names(data) == chid.name)
	altpos <- which(names(data) == alt.name)
	if (! is.null(id.var)){
		idpos <- which(names(data) == id.var)
		id.var <- as.factor(data[[id.var]])
	}
	if (! is.null(group.var)){
		grouppos <- which(names(data) == group.var)
		group.var <- as.factor(data[[group.var]])
	}

	if (drop.index){
		if (! is.null(id.var)) data <- data[, -c(chidpos, altpos, idpos)]
		else data <- data[, -c(chidpos, altpos)]
	}

	index <- data.frame(chid = chid, alt = alt)
	if (! is.null(id.var)) index <- cbind(index, id = id.var)
	if (! is.null(group.var)) index <- cbind(index, group = group.var)
	rownames(index) <- rownames(data)
	attr(data, "index") <- index
	attr(data, "class") <- c("mdcev.data", "data.frame")
	if (! is.null(choice)) attr(data, "choice") <- choice.name

	mdcev.datacheck(data)
	return(data)
}

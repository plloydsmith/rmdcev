.onAttach <- function(...) {
	ver <- utils::packageVersion("rmdcev")
	packageStartupMessage("This is rmdcev version ", ver)
}

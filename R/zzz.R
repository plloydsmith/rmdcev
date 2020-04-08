.onAttach <- function(...) {
	ver <- utils::packageVersion("rmdcev")
	packageStartupMessage("This is rmdcev version ", ver)
	packageStartupMessage("- This package continues to be under development so may change")
}

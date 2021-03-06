% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperFunctions.R
\name{CreateBlankPolicies}
\alias{CreateBlankPolicies}
\title{CreateBlankPolicies}
\usage{
CreateBlankPolicies(npols, model, price_change_only = TRUE)
}
\arguments{
\item{npols}{Number of policies to simulate}

\item{model}{Estimated model from mdcev}

\item{price_change_only}{Logical value for whether to include policy changes to dat_psi. Defaults to TRUE.
TRUE implies that only price changes are used in simulation.}
}
\description{
Create 'zero effect' policies that can be modified
}
\examples{
\donttest{
data_rec <- mdcev.data(data_rec, subset = id <= 500, id.var = "id",
                alt.var = "alt", choice = "quant")

mdcev_est <- mdcev( ~ 0, data = data_rec,
               model = "hybrid0", algorithm = "MLE",
               std_errors = "mvn")
CreateBlankPolicies(npols = 2, mdcev_est)
}
}

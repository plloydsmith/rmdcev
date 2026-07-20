% Generated manually - do not overwrite with roxygen2
\name{simulation_functions}
\alias{simulation_functions}
\alias{DrawError_rng}
\alias{MarshallianDemand}
\alias{HicksianDemand}
\alias{ComputeUtilJ}
\alias{CalcWTP_rng}
\alias{CalcMarshallianDemand_rng}
\title{Low-level MDCEV simulation functions}
\description{
Compiled C++ routines that underlie \code{\link{mdcev.sim}} and
\code{\link{GenerateMDCEVData}}.  These functions are exported for
advanced users and test code; most users should call the higher-level
wrappers instead.

\describe{
  \item{\code{DrawError_rng}}{Draw preference error terms for one
    individual using the MLHS or independent draw approach.}
  \item{\code{MarshallianDemand}}{Compute Marshallian (uncompensated)
    demand for one individual given marginal utilities at zero.}
  \item{\code{HicksianDemand}}{Compute Hicksian (compensated) demand
    for one individual given a reference utility level.}
  \item{\code{ComputeUtilJ}}{Evaluate the KT utility function for one
    individual.}
  \item{\code{CalcWTP_rng}}{Simulate welfare (WTP) across
    \code{nerrs} error draws for one individual under a price policy.}
  \item{\code{CalcMarshallianDemand_rng}}{Simulate Marshallian demand
    across \code{nerrs} error draws for one individual under a price
    policy.}
}
}
\usage{
DrawError_rng(quant_num, quant_j, price_j, psi_j, phi_j, gamma_j,
              alpha, scale, model_num, nalts, nerrs, cond_error,
              draw_mlhs, base_rng__, pstream__)

MarshallianDemand(income, price, MUzero, phi, gamma, alpha, nalts,
                  algo_gen, model_num, tol_e, max_loop, pstream__)

HicksianDemand(util, price, MUzero, phi, gamma, alpha, nalts,
               algo_gen, model_num, tol_l, max_loop, pstream__)

ComputeUtilJ(income, quant_j, price_j, psi, phi_j, gamma_j, alpha,
             nalts, model_num, pstream__)

CalcWTP_rng(income, quant_j, price, price_p_policy, psi_p_sims,
            phi_p_sims, psi_sims, phi_sims, gamma_sims, alpha_sims,
            scale_sims, nerrs, cond_error, draw_mlhs, algo_gen,
            model_num, price_change_only, tol, max_loop, base_rng__,
            pstream__)

CalcMarshallianDemand_rng(income, quant_j, price, price_p_policy,
                          psi_p_sims, phi_p_sims, psi_sims, phi_sims,
                          gamma_sims, alpha_sims, scale_sims, nerrs,
                          cond_error, draw_mlhs, algo_gen, model_num,
                          price_change_only, tol, max_loop, base_rng__,
                          pstream__)
}
\arguments{
  \item{base_rng__}{External pointer to the Boost RNG created by
    \code{\link{rmdcev_get_rng}}.}
  \item{pstream__}{External pointer to the output stream created by
    \code{\link{rmdcev_get_stream}}.}
  \item{model_num}{Integer model type: 1 = gamma, 2 = alpha,
    3 = hybrid, 4 = hybrid0, 5 = kt_ee, 6 = gamma1.}
  \item{income}{Individual income (scalar).}
  \item{quant_j}{Vector of non-numeraire quantities.}
  \item{quant_num}{Numeraire quantity consumed.}
  \item{price_j}{Vector of non-numeraire prices.}
  \item{price}{Full price vector (numeraire first).}
  \item{price_p_policy}{Matrix of counterfactual prices, one row per
    policy scenario.}
  \item{psi_j, psi}{Psi utility parameters.}
  \item{phi_j, phi}{Phi parameters (kt_ee model).}
  \item{gamma_j, gamma}{Gamma satiation parameters.}
  \item{alpha}{Alpha parameters.}
  \item{scale}{Scale parameter.}
  \item{psi_sims, phi_sims, gamma_sims, alpha_sims, scale_sims}{Lists of
    parameter draws, one element per simulation, for the baseline
    scenario.}
  \item{psi_p_sims, phi_p_sims}{Lists of psi and phi parameter draws for
    the policy scenarios.}
  \item{nalts}{Number of non-numeraire alternatives.}
  \item{nerrs}{Number of error draws.}
  \item{cond_error}{1 = conditional error draws; 0 = unconditional.}
  \item{draw_mlhs}{1 = Modified Latin Hypercube Sampling; 0 = iid.}
  \item{algo_gen}{Algorithm: 0 = hybrid (fast, requires alpha = 0);
    1 = general bisection.}
  \item{price_change_only}{1 = only prices change across policies, so
    the baseline psi and phi draws are reused; 0 = psi and phi also
    change across policies.}
  \item{MUzero}{Vector of marginal utilities at zero consumption.}
  \item{util}{Reference utility level (Hicksian demand).}
  \item{tol, tol_e, tol_l}{Convergence tolerance for bisection.}
  \item{max_loop}{Maximum bisection iterations.}
}
\value{
\code{DrawError_rng} returns a matrix of error draws with \code{nerrs}
rows and \code{nalts + 1} columns.

\code{MarshallianDemand} and \code{HicksianDemand} return a numeric
vector of length \code{nalts + 1} holding the numeraire quantity
followed by the non-numeraire quantities.

\code{ComputeUtilJ} returns a scalar utility value.

\code{CalcWTP_rng} and \code{CalcMarshallianDemand_rng} return a list
with one element per policy scenario, each holding the simulated
welfare measure or demand averaged over the \code{nerrs} error draws.
}
\seealso{
  \code{\link{mdcev.sim}}, \code{\link{PrepareSimulationData}},
  \code{\link{rmdcev_get_rng}}, \code{\link{rmdcev_get_stream}}
}

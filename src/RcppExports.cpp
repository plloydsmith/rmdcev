// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "rmdcev_types.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// Shuffle_rng
Eigen::Matrix<double, Eigen::Dynamic, 1> Shuffle_rng(const Eigen::Matrix<double, Eigen::Dynamic, 1>& inv_temp, const int& nerrs, boost::ecuyer1988& base_rng__, std::ostream* pstream__);
RcppExport SEXP _rmdcev_Shuffle_rng(SEXP inv_tempSEXP, SEXP nerrsSEXP, SEXP base_rng__SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type inv_temp(inv_tempSEXP);
    Rcpp::traits::input_parameter< const int& >::type nerrs(nerrsSEXP);
    Rcpp::traits::input_parameter< boost::ecuyer1988& >::type base_rng__(base_rng__SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(Shuffle_rng(inv_temp, nerrs, base_rng__, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// DrawMlhs_rng
Eigen::Matrix<double, Eigen::Dynamic, 1> DrawMlhs_rng(const int& nerrs, const int& draw_mlhs, boost::ecuyer1988& base_rng__, std::ostream* pstream__);
RcppExport SEXP _rmdcev_DrawMlhs_rng(SEXP nerrsSEXP, SEXP draw_mlhsSEXP, SEXP base_rng__SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type nerrs(nerrsSEXP);
    Rcpp::traits::input_parameter< const int& >::type draw_mlhs(draw_mlhsSEXP);
    Rcpp::traits::input_parameter< boost::ecuyer1988& >::type base_rng__(base_rng__SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(DrawMlhs_rng(nerrs, draw_mlhs, base_rng__, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// DrawError_rng
std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> > DrawError_rng(const double& quant_num, const Eigen::Matrix<double, Eigen::Dynamic, 1>& quant_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& psi_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& gamma_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, const double& scale, const int& model_num, const int& nalts, const int& nerrs, const int& cond_error, const int& draw_mlhs, boost::ecuyer1988& base_rng__, std::ostream* pstream__);
RcppExport SEXP _rmdcev_DrawError_rng(SEXP quant_numSEXP, SEXP quant_jSEXP, SEXP price_jSEXP, SEXP psi_jSEXP, SEXP phi_jSEXP, SEXP gamma_jSEXP, SEXP alphaSEXP, SEXP scaleSEXP, SEXP model_numSEXP, SEXP naltsSEXP, SEXP nerrsSEXP, SEXP cond_errorSEXP, SEXP draw_mlhsSEXP, SEXP base_rng__SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type quant_num(quant_numSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type quant_j(quant_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price_j(price_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type psi_j(psi_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type phi_j(phi_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< const int& >::type nalts(naltsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nerrs(nerrsSEXP);
    Rcpp::traits::input_parameter< const int& >::type cond_error(cond_errorSEXP);
    Rcpp::traits::input_parameter< const int& >::type draw_mlhs(draw_mlhsSEXP);
    Rcpp::traits::input_parameter< boost::ecuyer1988& >::type base_rng__(base_rng__SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(DrawError_rng(quant_num, quant_j, price_j, psi_j, phi_j, gamma_j, alpha, scale, model_num, nalts, nerrs, cond_error, draw_mlhs, base_rng__, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// CalcAltOrder
std::vector<int> CalcAltOrder(const Eigen::Matrix<double, Eigen::Dynamic, 1>& MUzero, const int& nalts, std::ostream* pstream__);
RcppExport SEXP _rmdcev_CalcAltOrder(SEXP MUzeroSEXP, SEXP naltsSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type MUzero(MUzeroSEXP);
    Rcpp::traits::input_parameter< const int& >::type nalts(naltsSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(CalcAltOrder(MUzero, nalts, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// SortParmMatrix
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> SortParmMatrix(const Eigen::Matrix<double, Eigen::Dynamic, 1>& MUzero, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price, const Eigen::Matrix<double, Eigen::Dynamic, 1>& gamma, const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha_phi, const int& nalts, std::ostream* pstream__);
RcppExport SEXP _rmdcev_SortParmMatrix(SEXP MUzeroSEXP, SEXP priceSEXP, SEXP gammaSEXP, SEXP alpha_phiSEXP, SEXP naltsSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type MUzero(MUzeroSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price(priceSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type alpha_phi(alpha_phiSEXP);
    Rcpp::traits::input_parameter< const int& >::type nalts(naltsSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(SortParmMatrix(MUzero, price, gamma, alpha_phi, nalts, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// ComputeE
double ComputeE(const int& M, const double& lambda, const Eigen::Matrix<double, Eigen::Dynamic, 1>& g_price, const Eigen::Matrix<double, Eigen::Dynamic, 1>& b, const Eigen::Matrix<double, Eigen::Dynamic, 1>& c, const Eigen::Matrix<double, Eigen::Dynamic, 1>& d, std::ostream* pstream__);
RcppExport SEXP _rmdcev_ComputeE(SEXP MSEXP, SEXP lambdaSEXP, SEXP g_priceSEXP, SEXP bSEXP, SEXP cSEXP, SEXP dSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type g_price(g_priceSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type b(bSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type d(dSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeE(M, lambda, g_price, b, c, d, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// ComputeKtE
double ComputeKtE(const int& M, const double& lambda, const Eigen::Matrix<double, Eigen::Dynamic, 1>& mu, const Eigen::Matrix<double, Eigen::Dynamic, 1>& g_price__phi, const double& alpha_1, std::ostream* pstream__);
RcppExport SEXP _rmdcev_ComputeKtE(SEXP MSEXP, SEXP lambdaSEXP, SEXP muSEXP, SEXP g_price__phiSEXP, SEXP alpha_1SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type g_price__phi(g_price__phiSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_1(alpha_1SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeKtE(M, lambda, mu, g_price__phi, alpha_1, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// MarshallianDemand
Eigen::Matrix<double, Eigen::Dynamic, 1> MarshallianDemand(const double& income, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price, const Eigen::Matrix<double, Eigen::Dynamic, 1>& MUzero, const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi, const Eigen::Matrix<double, Eigen::Dynamic, 1>& gamma, const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, const int& nalts, const int& algo_gen, const int& model_num, const double& tol_e, const int& max_loop, std::ostream* pstream__);
RcppExport SEXP _rmdcev_MarshallianDemand(SEXP incomeSEXP, SEXP priceSEXP, SEXP MUzeroSEXP, SEXP phiSEXP, SEXP gammaSEXP, SEXP alphaSEXP, SEXP naltsSEXP, SEXP algo_genSEXP, SEXP model_numSEXP, SEXP tol_eSEXP, SEXP max_loopSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price(priceSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type MUzero(MUzeroSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type nalts(naltsSEXP);
    Rcpp::traits::input_parameter< const int& >::type algo_gen(algo_genSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol_e(tol_eSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_loop(max_loopSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(MarshallianDemand(income, price, MUzero, phi, gamma, alpha, nalts, algo_gen, model_num, tol_e, max_loop, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// ComputeUtilJ
double ComputeUtilJ(const double& income, const Eigen::Matrix<double, Eigen::Dynamic, 1>& quant_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& psi_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& gamma_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, const int& nalts, const int& model_num, std::ostream* pstream__);
RcppExport SEXP _rmdcev_ComputeUtilJ(SEXP incomeSEXP, SEXP quant_jSEXP, SEXP price_jSEXP, SEXP psi_jSEXP, SEXP phi_jSEXP, SEXP gamma_jSEXP, SEXP alphaSEXP, SEXP naltsSEXP, SEXP model_numSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type quant_j(quant_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price_j(price_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type psi_j(psi_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type phi_j(phi_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type nalts(naltsSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeUtilJ(income, quant_j, price_j, psi_j, phi_j, gamma_j, alpha, nalts, model_num, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// ComputeUtilM
double ComputeUtilM(const int& M, const double& lambda1, const Eigen::Matrix<double, Eigen::Dynamic, 1>& g_psi_a, const Eigen::Matrix<double, Eigen::Dynamic, 1>& a_a_1, const Eigen::Matrix<double, Eigen::Dynamic, 1>& mu_a_a_1, const Eigen::Matrix<double, Eigen::Dynamic, 1>& psi, const Eigen::Matrix<double, Eigen::Dynamic, 1>& g, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price, const Eigen::Matrix<double, Eigen::Dynamic, 1>& d, const int& model_num, std::ostream* pstream__);
RcppExport SEXP _rmdcev_ComputeUtilM(SEXP MSEXP, SEXP lambda1SEXP, SEXP g_psi_aSEXP, SEXP a_a_1SEXP, SEXP mu_a_a_1SEXP, SEXP psiSEXP, SEXP gSEXP, SEXP priceSEXP, SEXP dSEXP, SEXP model_numSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type g_psi_a(g_psi_aSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type a_a_1(a_a_1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type mu_a_a_1(mu_a_a_1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type g(gSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price(priceSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeUtilM(M, lambda1, g_psi_a, a_a_1, mu_a_a_1, psi, g, price, d, model_num, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// ComputeKtUtilM
double ComputeKtUtilM(const int& M, const double& lambda1, const Eigen::Matrix<double, Eigen::Dynamic, 1>& psi, const Eigen::Matrix<double, Eigen::Dynamic, 1>& mu, const double& alpha_1, std::ostream* pstream__);
RcppExport SEXP _rmdcev_ComputeKtUtilM(SEXP MSEXP, SEXP lambda1SEXP, SEXP psiSEXP, SEXP muSEXP, SEXP alpha_1SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha_1(alpha_1SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeKtUtilM(M, lambda1, psi, mu, alpha_1, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// HicksianDemand
Eigen::Matrix<double, Eigen::Dynamic, 1> HicksianDemand(const double& util, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price, const Eigen::Matrix<double, Eigen::Dynamic, 1>& MUzero, const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi, const Eigen::Matrix<double, Eigen::Dynamic, 1>& gamma, const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, const int& nalts, const int& algo_gen, const int& model_num, const double& tol_l, const int& max_loop, std::ostream* pstream__);
RcppExport SEXP _rmdcev_HicksianDemand(SEXP utilSEXP, SEXP priceSEXP, SEXP MUzeroSEXP, SEXP phiSEXP, SEXP gammaSEXP, SEXP alphaSEXP, SEXP naltsSEXP, SEXP algo_genSEXP, SEXP model_numSEXP, SEXP tol_lSEXP, SEXP max_loopSEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type util(utilSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price(priceSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type MUzero(MUzeroSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type nalts(naltsSEXP);
    Rcpp::traits::input_parameter< const int& >::type algo_gen(algo_genSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol_l(tol_lSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_loop(max_loopSEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(HicksianDemand(util, price, MUzero, phi, gamma, alpha, nalts, algo_gen, model_num, tol_l, max_loop, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// CalcWTP_rng
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> CalcWTP_rng(const double& income, const Eigen::Matrix<double, Eigen::Dynamic, 1>& quant_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& price_p_policy, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& psi_p_sims, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& phi_p_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& psi_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& phi_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& gamma_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& alpha_sims, const Eigen::Matrix<double, Eigen::Dynamic, 1>& scale_sims, const int& nerrs, const int& cond_error, const int& draw_mlhs, const int& algo_gen, const int& model_num, const int& price_change_only, const double& tol, const int& max_loop, boost::ecuyer1988& base_rng__, std::ostream* pstream__);
RcppExport SEXP _rmdcev_CalcWTP_rng(SEXP incomeSEXP, SEXP quant_jSEXP, SEXP priceSEXP, SEXP price_p_policySEXP, SEXP psi_p_simsSEXP, SEXP phi_p_simsSEXP, SEXP psi_simsSEXP, SEXP phi_simsSEXP, SEXP gamma_simsSEXP, SEXP alpha_simsSEXP, SEXP scale_simsSEXP, SEXP nerrsSEXP, SEXP cond_errorSEXP, SEXP draw_mlhsSEXP, SEXP algo_genSEXP, SEXP model_numSEXP, SEXP price_change_onlySEXP, SEXP tolSEXP, SEXP max_loopSEXP, SEXP base_rng__SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type quant_j(quant_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price(priceSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& >::type price_p_policy(price_p_policySEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& >::type psi_p_sims(psi_p_simsSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& >::type phi_p_sims(phi_p_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type psi_sims(psi_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type phi_sims(phi_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type gamma_sims(gamma_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type alpha_sims(alpha_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type scale_sims(scale_simsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nerrs(nerrsSEXP);
    Rcpp::traits::input_parameter< const int& >::type cond_error(cond_errorSEXP);
    Rcpp::traits::input_parameter< const int& >::type draw_mlhs(draw_mlhsSEXP);
    Rcpp::traits::input_parameter< const int& >::type algo_gen(algo_genSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< const int& >::type price_change_only(price_change_onlySEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_loop(max_loopSEXP);
    Rcpp::traits::input_parameter< boost::ecuyer1988& >::type base_rng__(base_rng__SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(CalcWTP_rng(income, quant_j, price, price_p_policy, psi_p_sims, phi_p_sims, psi_sims, phi_sims, gamma_sims, alpha_sims, scale_sims, nerrs, cond_error, draw_mlhs, algo_gen, model_num, price_change_only, tol, max_loop, base_rng__, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// CalcMarshallianDemand_rng
std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> > CalcMarshallianDemand_rng(const double& income, const Eigen::Matrix<double, Eigen::Dynamic, 1>& quant_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& price_p_policy, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& psi_p_sims, const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& phi_p_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& psi_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& phi_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& gamma_sims, const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& alpha_sims, const Eigen::Matrix<double, Eigen::Dynamic, 1>& scale_sims, const int& nerrs, const int& cond_error, const int& draw_mlhs, const int& algo_gen, const int& model_num, const int& price_change_only, const double& tol, const int& max_loop, boost::ecuyer1988& base_rng__, std::ostream* pstream__);
RcppExport SEXP _rmdcev_CalcMarshallianDemand_rng(SEXP incomeSEXP, SEXP quant_jSEXP, SEXP priceSEXP, SEXP price_p_policySEXP, SEXP psi_p_simsSEXP, SEXP phi_p_simsSEXP, SEXP psi_simsSEXP, SEXP phi_simsSEXP, SEXP gamma_simsSEXP, SEXP alpha_simsSEXP, SEXP scale_simsSEXP, SEXP nerrsSEXP, SEXP cond_errorSEXP, SEXP draw_mlhsSEXP, SEXP algo_genSEXP, SEXP model_numSEXP, SEXP price_change_onlySEXP, SEXP tolSEXP, SEXP max_loopSEXP, SEXP base_rng__SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type quant_j(quant_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price(priceSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1> >& >::type price_p_policy(price_p_policySEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& >::type psi_p_sims(psi_p_simsSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> >& >::type phi_p_sims(phi_p_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type psi_sims(psi_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type phi_sims(phi_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type gamma_sims(gamma_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& >::type alpha_sims(alpha_simsSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type scale_sims(scale_simsSEXP);
    Rcpp::traits::input_parameter< const int& >::type nerrs(nerrsSEXP);
    Rcpp::traits::input_parameter< const int& >::type cond_error(cond_errorSEXP);
    Rcpp::traits::input_parameter< const int& >::type draw_mlhs(draw_mlhsSEXP);
    Rcpp::traits::input_parameter< const int& >::type algo_gen(algo_genSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< const int& >::type price_change_only(price_change_onlySEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_loop(max_loopSEXP);
    Rcpp::traits::input_parameter< boost::ecuyer1988& >::type base_rng__(base_rng__SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(CalcMarshallianDemand_rng(income, quant_j, price, price_p_policy, psi_p_sims, phi_p_sims, psi_sims, phi_sims, gamma_sims, alpha_sims, scale_sims, nerrs, cond_error, draw_mlhs, algo_gen, model_num, price_change_only, tol, max_loop, base_rng__, pstream__));
    return rcpp_result_gen;
END_RCPP
}
// CalcmdemandOne_rng
Eigen::Matrix<double, Eigen::Dynamic, 1> CalcmdemandOne_rng(const double& income, const Eigen::Matrix<double, Eigen::Dynamic, 1>& price, const Eigen::Matrix<double, Eigen::Dynamic, 1>& psi_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& phi_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& gamma_j, const Eigen::Matrix<double, Eigen::Dynamic, 1>& alpha, const double& scale, const int& nerrs, const int& model_num, const int& algo_gen, const double& tol, const int& max_loop, boost::ecuyer1988& base_rng__, std::ostream* pstream__);
RcppExport SEXP _rmdcev_CalcmdemandOne_rng(SEXP incomeSEXP, SEXP priceSEXP, SEXP psi_jSEXP, SEXP phi_jSEXP, SEXP gamma_jSEXP, SEXP alphaSEXP, SEXP scaleSEXP, SEXP nerrsSEXP, SEXP model_numSEXP, SEXP algo_genSEXP, SEXP tolSEXP, SEXP max_loopSEXP, SEXP base_rng__SEXP, SEXP pstream__SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type income(incomeSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type price(priceSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type psi_j(psi_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type phi_j(phi_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type gamma_j(gamma_jSEXP);
    Rcpp::traits::input_parameter< const Eigen::Matrix<double, Eigen::Dynamic, 1>& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double& >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< const int& >::type nerrs(nerrsSEXP);
    Rcpp::traits::input_parameter< const int& >::type model_num(model_numSEXP);
    Rcpp::traits::input_parameter< const int& >::type algo_gen(algo_genSEXP);
    Rcpp::traits::input_parameter< const double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_loop(max_loopSEXP);
    Rcpp::traits::input_parameter< boost::ecuyer1988& >::type base_rng__(base_rng__SEXP);
    Rcpp::traits::input_parameter< std::ostream* >::type pstream__(pstream__SEXP);
    rcpp_result_gen = Rcpp::wrap(CalcmdemandOne_rng(income, price, psi_j, phi_j, gamma_j, alpha, scale, nerrs, model_num, algo_gen, tol, max_loop, base_rng__, pstream__));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4mdcev_mod();
RcppExport SEXP _rcpp_module_boot_stan_fit4mdcev_rp_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_rcpp_module_boot_stan_fit4mdcev_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mdcev_mod, 3},
    {"_rcpp_module_boot_stan_fit4mdcev_rp_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4mdcev_rp_mod, 3},
    {"_rmdcev_Shuffle_rng", (DL_FUNC) &_rmdcev_Shuffle_rng, 4},
    {"_rmdcev_DrawMlhs_rng", (DL_FUNC) &_rmdcev_DrawMlhs_rng, 4},
    {"_rmdcev_DrawError_rng", (DL_FUNC) &_rmdcev_DrawError_rng, 15},
    {"_rmdcev_CalcAltOrder", (DL_FUNC) &_rmdcev_CalcAltOrder, 3},
    {"_rmdcev_SortParmMatrix", (DL_FUNC) &_rmdcev_SortParmMatrix, 6},
    {"_rmdcev_ComputeE", (DL_FUNC) &_rmdcev_ComputeE, 7},
    {"_rmdcev_ComputeKtE", (DL_FUNC) &_rmdcev_ComputeKtE, 6},
    {"_rmdcev_MarshallianDemand", (DL_FUNC) &_rmdcev_MarshallianDemand, 12},
    {"_rmdcev_ComputeUtilJ", (DL_FUNC) &_rmdcev_ComputeUtilJ, 10},
    {"_rmdcev_ComputeUtilM", (DL_FUNC) &_rmdcev_ComputeUtilM, 11},
    {"_rmdcev_ComputeKtUtilM", (DL_FUNC) &_rmdcev_ComputeKtUtilM, 6},
    {"_rmdcev_HicksianDemand", (DL_FUNC) &_rmdcev_HicksianDemand, 12},
    {"_rmdcev_CalcWTP_rng", (DL_FUNC) &_rmdcev_CalcWTP_rng, 21},
    {"_rmdcev_CalcMarshallianDemand_rng", (DL_FUNC) &_rmdcev_CalcMarshallianDemand_rng, 21},
    {"_rmdcev_CalcmdemandOne_rng", (DL_FUNC) &_rmdcev_CalcmdemandOne_rng, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_rmdcev(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

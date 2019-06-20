
int I; // number of Individuals
int J; // number of non-numeraire alternatives
int NPsi; // number of alt covariates
matrix[I * J, NPsi] dat_psi; // alt characteristics
matrix[I, J] price_j; // non-numeraire price
matrix[I, J] quant_j; // non-numeraire consumption
vector[I] income;
int<lower=0, upper=1> flat_priors; //indicator to include flat priors for all parameters
real prior_psi_sd;
real prior_gamma_sd;
real prior_alpha_sd;
real prior_scale_sd;
int<lower = 1, upper = 4> model_num; // 1 is les, 2 is alpha, 3 gamma (one alpha for all), 4 alpha's set to 1e-03
int<lower=0, upper=1> fixed_scale1; // indicator to fix scale at 1
int<lower=0, upper=1> trunc_data; //indicator to correct estimation for truncated data
vector[I] weights; // user supplied weights

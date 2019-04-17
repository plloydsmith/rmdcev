
int I; // number of Individuals
int J; // number of non-numeraire alternatives
int NPsi; // number of alt covariates
int K; // number of mixtures
matrix[I * J, NPsi] dat_psi; // alt characteristics
matrix[I, J] j_price; // non-numeraire price
matrix[I, J] j_quant; // non-numeraire consumption
vector[I] income;
vector[I] num_price; // numeraire price
vector[I] M_factorial; // (M-1)!
real prior_psi_sd;
real prior_gamma_sd;
real prior_alpha_sd;
real prior_scale_sd;
int<lower = 1, upper = 4> model_num; // 1 is les, 2 is alpha, 3 gamma (one alpha for all), 4 alpha's set to 1e-06
int<lower=0, upper=1> fixed_scale; // indicator to fix scale
int<lower=0, upper=1> trunc_data; //indicator to correct estimation for truncated data
int<lower=0, upper=1> no_priors; //indicator to include priors or not
vector[I] weights; // user supplied weights

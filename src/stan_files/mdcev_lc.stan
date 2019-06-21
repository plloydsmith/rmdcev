// saved as mdcev_hb_corr.stan
functions {
#include /common/mdcev_ll.stan
}


data {
#include /common/mdcev_data.stan
  int<lower=0, upper=1> corr;
  int task[I]; // index for tasks
  int task_individual[I]; // index for individual
  int start[I]; // the starting observation for each task
  int end[I]; // the ending observation for each task
  real<lower=1> lkj_shape; // shape parameter for LKJ prior
//  vector[NPsi] psi_ndx;
  int<lower=0, upper=1>  gamma_fixed;
  int<lower=0, upper=1>  alpha_fixed;
}

transformed data {
	int RP;
	int RP_g;
	int RP_a;
#include /common/mdcev_tdata.stan
	{
	int n_gamma_rp = gamma_fixed == 0 ? Gamma : 0;
	int n_alpha_rp = alpha_fixed == 0 ? A : 0;
	RP = NPsi + n_gamma_rp + n_alpha_rp; // number of random parameters
	RP_g = NPsi + 1; // location of first gamma
	RP_a = NPsi + n_gamma_rp + 1; // location of first alpha
	}
//	n_psi_rp = sum(psi_ndx);
//	n_psi_fixed = NPsi - sum(psi_ndx);
}

parameters {
//	vector[n_psi_fixed] psi;
	vector<lower=0 >[gamma_fixed == 1 ? Gamma : 0] gamma;
	vector<lower=0, upper=1>[alpha_fixed == 1 ? A : 0] alpha;
	vector[RP] mu;                                // means for beta
  	matrix[I, RP] z;                             // std normal draws
	cholesky_factor_corr[corr == 1 ? RP : 0] L_Omega;                // cholesky factors
  	vector<lower=0,upper=pi()/2>[RP] tau_unif;
	vector<lower=0>[fixed_scale1 == 0 ? 1 : 0] scale;
}

transformed parameters {
	vector[I] log_like;
//	cholesky_factor_cov[RP] L;                       // cholesky factors
  	vector<lower=0>[RP] tau;   	// diagonal of the part-worth covariance matrix

  	{
	matrix[I, RP] beta;             // utility parameters (individual level)
	matrix[I, J] lpsi;
   	matrix[I, NPsi] psi_individual;
  	matrix[I, J+1] alpha_full;
  	vector[I] alpha_individual_1;
  	matrix[I, J] alpha_individual_j;
  	matrix[I, J] gamma_individual;
	real scale_full = fixed_scale1 == 0 ? scale[1] : 1.0;

	for (rp in 1:RP) tau[rp] = 2.5 * tan(tau_unif[rp]);

	// individual level parameters
	if (corr == 1){
		beta = rep_matrix(mu', I) + (z * diag_pre_multiply(tau, L_Omega));
	} else if (corr == 0){
		beta = rep_matrix(mu', I) + diag_post_multiply(z, tau);
	}

	if (alpha_fixed == 0){
		if (model_num < 4){
			alpha_individual_1 = inv(1 + exp(-col(beta, RP_a)));
			if (model_num == 1)
			  	alpha_individual_j = rep_matrix(0, I, J);
			else if (model_num == 2)
			  	alpha_individual_j = inv(1 + exp(-block(beta, 1, RP_a + 1, I, J)));
			else
				alpha_individual_j = rep_matrix(alpha_individual_1, J);
		} else {
			alpha_individual_1 = rep_vector(1e-03, I);
			alpha_individual_j = rep_matrix(1e-03, I, J);
		}
	} else if (alpha_fixed == 1){
		alpha_full = alpha_ll(alpha, I, J, model_num);
		alpha_individual_1 = col(alpha_full, 1);
		alpha_individual_j = block(alpha_full, 1, 2, I, J);
	}

	if (gamma_fixed == 0){
		if (model_num == 2)
	  		gamma_individual = rep_matrix(1, I, J);
		else
	  		gamma_individual = exp(block(beta, 1, RP_g, I, J));
	} else if (gamma_fixed == 1){
		gamma_individual = gamma_ll(gamma, I, J, model_num);
	}

	psi_individual = block(beta, 1, 1, I, NPsi);

	for(t in 1:I){
		row_vector[J] util;
		util = psi_individual[task_individual[t]] * dat_psi[start[t]:end[t]]';
		lpsi[t] = util;
	}

	log_like = mdcev_ll(quant_j, price_j, log_num, log_inc, M, log_M_fact, // data
			lpsi, gamma_individual, alpha_individual_1, alpha_individual_j, scale_full, 						// parameters
			I, J, nonzero, trunc_data);
	}
}

model {
  // priors on the parameters
  	gamma ~ normal(0, prior_gamma_sd);
	alpha ~ normal(.5, prior_alpha_sd);
	to_vector(z) ~ normal(0, 1);
	to_vector(mu) ~ normal(0, 10);
	L_Omega ~ lkj_corr_cholesky(lkj_shape);                 // lkj prior
	scale ~ normal(1, 1);
	// no priors for tau because already constrained to uniform
  target += sum(log_like .* weights);//objective to target
}

generated quantities{
   matrix[RP, RP] Sigma;                            // cov matrix
   real<upper=0> sum_log_lik = 0;// log_lik for each sample

	{
	   matrix[RP, RP] L;
	//	vector[RP] tau_trans = append_col(tau[1:NPsi], append_col(exp(tau[RP_g:(RP_g+Gamma)]),
	//									exp(tau[RP_a:(RP_a+A)])/(1+exp(tau[RP_a:(RP_a+A)]))));
		if (corr == 1){
			L = diag_pre_multiply(tau, L_Omega);
		} else if (corr == 0){
			L = diag_matrix(tau);
		}
		Sigma = tcrossprod(L);
	}

	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i];
	}
}

// saved as mdcev_hb_corr.stan
functions {
#include /common/mdcev_ll_matrix.stan
}


data {
#include /common/mdcev_data.stan
  int<lower=0, upper=1> corr;
  int task[I]; // index for tasks
  int task_individual[I]; // index for individual
  int start[I]; // the starting observation for each task
  int end[I]; // the ending observation for each task
  real<lower=1> lkj_shape; // shape parameter for LKJ prior
}

transformed data {
	int RP;
	int RP_g;
	int RP_a;
	int N_Omega;
#include /common/mdcev_tdata.stan
 	RP = NPsi + Gamma + A; // number of random parameters
	RP_g = NPsi + 1; // location of first gamma
	RP_a = NPsi + Gamma + 1; // location of alpha

	if (corr == 1){
	  	N_Omega = RP;
	} else if (corr == 0){
	  	N_Omega = 0;
	}
}

parameters {
	vector[RP] mu;                                // means for beta
  	matrix[I, RP] z;                             // std normal draws
	cholesky_factor_corr[N_Omega] L_Omega;                // cholesky factors
  	vector<lower=0,upper=pi()/2>[RP] tau_unif;
	vector<lower=0>[fixed_scale == 0 ? 1 : 0] scale;
}

transformed parameters {
	vector[I] log_like;
	cholesky_factor_cov[RP] L;                       // cholesky factors
  	vector<lower=0>[RP] tau;   	// diagonal of the part-worth covariance matrix
	matrix[I, RP] beta;             // utility parameters (individual level)
  	{
	matrix[I, J] lpsi;
  	matrix[I, NPsi] psi_individual;
  	matrix[I, G] alpha_individual;
  	matrix[I, G] gamma_individual;
	real scale_full = fixed_scale == 0 ? scale[1] : 1.0;

	for (rp in 1:RP) tau[rp] = 2.5 * tan(tau_unif[rp]);

	// individual level parameters
	if (corr == 1){
	  	L = diag_pre_multiply(tau, L_Omega);
	} else if (corr == 0){
	 	L = diag_matrix(tau);
	}
	beta = rep_matrix(mu', I) + (z * L);

	if (model_num == 1)
	  	alpha_individual = append_col(1 - exp(col(beta, RP_a)), rep_matrix(0, I, J));
	else if (model_num == 2)
	  	alpha_individual = 1 - exp(block(beta, 1, RP_a, I, G));
	else if (model_num == 3)
		alpha_individual = rep_matrix(1 - exp(col(beta, RP_a)), G);
	else
		alpha_individual = rep_matrix(1e-06, I, G);

	if (model_num == 2)
	  gamma_individual = append_col(rep_vector(0, I), rep_matrix(1, I, J));
	else
	  gamma_individual = append_col(rep_vector(0, I), exp(block(beta, 1, RP_g, I, J)));

	psi_individual = block(beta, 1, 1, I, NPsi);

	for(t in 1:I){
		row_vector[J] util;
		util = psi_individual[task_individual[t]] * dat_psi[start[t]:end[t]]';
		lpsi[t] = util;
	}

//	gamma_full = gamma_ll(gamma[1], I, J, G, model_num);
//	alpha_full = alpha_ll(alpha[1], I, J, G, model_num);

	log_like = mdcev_ll_matrix(j_quant, quant_full, log_price, log_num, price_full,
			log_inc, dat_psi, M, M_factorial, weights, // data
			lpsi, gamma_individual, alpha_individual, scale_full, 						// parameters
			I, J, G, ones_g, nonzero, model_num, fixed_scale, trunc_data);
	}
}

model {
  // priors on the parameters
	to_vector(z) ~ normal(0, 1);
	to_vector(mu) ~ normal(0, 10);
	L_Omega ~ lkj_corr_cholesky(lkj_shape);                 // lkj prior
	scale ~ normal(1, 1);
	// no priors for tau because already constrained to uniform
  target += sum(log_like);//objective to target
}

generated quantities{
  cov_matrix[RP] Sigma;                            // cov matrix
  real<upper=0> sum_log_lik = 0;                          // log_lik for each sample
  Sigma = tcrossprod(L);

	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i];
	}
}

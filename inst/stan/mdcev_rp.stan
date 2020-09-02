// saved as mdcev_lc.stan
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
  int<lower=0, upper=1>  gamma_nonrandom;
  int<lower=0, upper=1>  alpha_nonrandom;
}

transformed data {
	int RP;
	int RP_Phi;
	int RP_g;
	int RP_a;

#include /common/mdcev_tdata.stan
	{
	int n_gamma_rp = gamma_nonrandom == 0 ? Gamma : 0;
	int n_alpha_rp = alpha_nonrandom == 0 ? A : 0;
	RP = NPsi + NPhi + n_gamma_rp + n_alpha_rp; // number of random parameters
	RP_Phi = NPsi + 1; // location of first phi
	RP_g = RP_Phi + NPhi; // location of first gamma
	RP_a = RP_g + n_gamma_rp; // location of first alpha
	}
//	n_psi_rp = sum(psi_ndx);
//	n_psi_fixed = NPsi - sum(psi_ndx);
}

parameters {
//	vector[n_psi_fixed] psi;
	vector<lower=0 >[gamma_nonrandom == 1 ? Gamma : 0] gamma;
	vector<lower=0, upper=1>[alpha_nonrandom == 1 ? A : 0] alpha;
	vector[RP] mu;                                // means for random parameters
  	matrix[I, RP] z;                             // std normal draws
	cholesky_factor_corr[corr == 1 ? RP : 0] L_Omega; // cholesky factors
  	vector<lower=0,upper=pi()/2>[RP] tau_unif;
	vector<lower=0>[fixed_scale1 == 0 ? 1 : 0] scale;
}

transformed parameters {
	vector[I] log_like;
  	vector<lower=0>[RP] tau;   	// diagonal of the covariance matrix

  	{
	matrix[I, RP] beta;             // individual level parameters
	matrix[I, J] lpsi;
  	vector[I] alpha_individual_1;
  	matrix[I, J] alpha_individual_j;
  	matrix[I, J] gamma_individual;
	real scale_full = fixed_scale1 == 0 ? scale[1] : 1.0;

	for (rp in 1:RP) tau[rp] = 2.5 * tan(tau_unif[rp]);

	// individual level parameters
	if (corr == 1){
		beta = rep_matrix(mu', I) + z * diag_pre_multiply(tau, L_Omega);
	} else if (corr == 0){
		beta = rep_matrix(mu', I) + diag_post_multiply(z, tau);
	}

	if (alpha_nonrandom == 0){
		if (model_num != 4){
			alpha_individual_1 = inv_logit(col(beta, RP_a));
			if (model_num == 1)
			  	alpha_individual_j = rep_matrix(0, I, J);
			else if (model_num == 2)
			  	alpha_individual_j = inv_logit(block(beta, 1, RP_a + 1, I, J));
			else if  (model_num == 3)
				alpha_individual_j = rep_matrix(alpha_individual_1, J);
		} else {
			alpha_individual_1 = rep_vector(1e-03, I);
			alpha_individual_j = rep_matrix(1e-03, I, J);
		}
	} else if (alpha_nonrandom == 1){
		matrix[I, J+1] alpha_full = alpha_ll(alpha, I, J, model_num);
		alpha_individual_1 = col(alpha_full, 1);
		alpha_individual_j = block(alpha_full, 1, 2, I, J);
	}

	if (gamma_nonrandom == 0 && model_num != 2){
		if (gamma_ascs == 1)
			gamma_individual = exp(block(beta, 1, RP_g, I, Gamma));
		else if (gamma_ascs == 0)
			gamma_individual = rep_matrix(exp(col(beta, RP_g)), J);
	} else
		gamma_individual = gamma_ll(gamma, I, J, Gamma);

	if (psi_ascs == 1){
		lpsi = append_col(rep_vector(0, I), block(beta, 1, 1, I, J-1)); // pull out acs
		if (NPsi_ij > 0)
			for(i in 1:I)
				lpsi[i] = lpsi[i] + sub_row(beta, task_individual[i], J, NPsi_ij) * dat_psi[start[i]:end[i]]';
	} else if (psi_ascs == 0){
		if (NPsi_ij > 0)
			for(i in 1:I)
				lpsi[i] = sub_row(beta, task_individual[i], 1, NPsi_ij) * dat_psi[start[i]:end[i]]';
		else if (NPsi_ij == 0)
			lpsi = rep_matrix(0, I, J);
	}

	if (model_num < 5){
		log_like = mdcev_ll(quant_j, price_j, log_num, income, M, log_M_fact, // data
			lpsi, gamma_individual, alpha_individual_1, alpha_individual_j, scale_full, 						// parameters
			I, J, nonzero, trunc_data);
	} else if (model_num == 5){
		matrix[I, J] phi_ij;
		if (NPhi == 0)
			phi_ij = rep_matrix(1, I, J);
		else if (NPhi > 0)
			for(i in 1:I)
				phi_ij[i] = exp(sub_row(beta, task_individual[i], RP_Phi, NPhi) * dat_phi[start[i]:end[i]]');

		log_like = kt_ll(quant_j, price_j, log_num, income,
  					lpsi, phi_ij, gamma_individual, alpha_individual_1, scale_full,
  					I, J, nonzero, trunc_data, jacobian_analytical_grad);
	}
  	}
}

model {
  // priors on the parameters
  	gamma ~ normal(1, prior_gamma_sd);
	alpha ~ beta(prior_alpha_shape, prior_alpha_shape);
	to_vector(z) ~ std_normal();
	to_vector(mu) ~ normal(0, 10);
	L_Omega ~ lkj_corr_cholesky(lkj_shape);                 // lkj prior
	scale ~ normal(0, 1);
	// no priors for tau because already constrained to uniform
  target += dot_product(log_like, weights);//objective to target
}

generated quantities{
   matrix[RP, RP] Sigma;                            // cov matrix
   real<upper=0> sum_log_lik = 0;// log_lik for each sample

	{
	   matrix[RP, RP] L;
		if (corr == 1){
			matrix[RP, RP] Omega;
			Omega = multiply_lower_tri_self_transpose(L_Omega);  // correlation matrix
			Sigma = quad_form_diag(Omega, tau);               // var-covar matrix
		} else if (corr == 0){
			Sigma = diag_matrix(tau);
		}
	}

	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i] * weights[i];
	}
}

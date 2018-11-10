// saved as mdcev_fixed.stan
data {
  int I; // number of Individuals
  int J; // number of non-numeraire alternatives
  int IJ; // number of rows
  int NPsi; // number of alt covariates
  matrix[IJ, NPsi] dat_psi; // alt characteristics
  matrix[I, J] j_price; // non-numeraire price
  matrix[I, J] j_quant; // non-numeraire consumption
  vector[I] num_price; // numeraire price
  vector[I] num_quant; // numeraire consumption
  int<lower = 1, upper = 4> model_type; // 1 is les, 2 is alpha, 3 gamma (one alpha for all), 4 alpha's set to 0
  int<lower=0, upper=1> fixed_scale; // indicator to fix scale
  vector[I] weights; // user supplied weights
	vector[I] M_factorial;
	matrix[I, J + 1] nonzero;
	matrix[I, J + 1] zero;
	vector[I] M;	//  Number of consumed goods (including numeraire)
  int task[I]; // index for tasks
  int task_individual[I]; // index for individual
  int start[I]; // the starting observation for each task
  int end[I]; // the ending observation for each task
}

transformed data {
	int G = J + 1;
	int A;
	int Gamma;
	int RP;
	int RP_g;
	int RP_a;
	vector[G] ones_g = rep_vector(1, G);
	vector[NPsi] ones_npsi = rep_vector(1, NPsi);
	matrix[I, G] price_full = append_col(num_price, j_price);
	matrix[I, G] quant_full = append_col(num_quant, j_quant);
	matrix[I, J] log_price = log(j_price);
	vector[I] log_num = log(num_quant ./ num_price);

 	if (model_type == 1 || model_type == 3){
 		A = 1;
 		Gamma = J;
	}else if (model_type == 2){
 		A = G;
 		Gamma = 0;
	}else if (model_type == 4){
		A = 0;
 		Gamma = J;
	}
 	RP = NPsi + Gamma + A; // number of random parameters
	RP_g = NPsi + 1; // location of first gamma
	RP_a = NPsi + Gamma + 1; // location of alpha

}

parameters {
	vector[RP] mu;                                // means for beta
  	matrix[I, RP] z;                             // std normal draws
	cholesky_factor_corr[RP] L_Omega;                // cholesky factors
  	vector<lower=0,upper=pi()/2>[RP] tau_unif;
	vector<lower=0>[fixed_scale == 0 ? 1 : 0] scale;
}

transformed parameters {
	vector[I] log_like;
	cholesky_factor_cov[RP] L;                       // cholesky factors
  	vector<lower=0>[RP] tau;   	// diagonal of the part-worth covariance matrix
  	{
	matrix[I, J] lpsi;
	matrix[I, G] f;
	matrix[I, G] v;
	matrix[I, G] vf;
	vector[I] sumv;
	vector[I] pf;
	vector[I] prodvf;
	vector[G] gamma_full;
	vector[G] alpha_full;
	matrix[I, RP] beta;             // utility parameters (individual level)
  	matrix[I, NPsi] psi_individual;
  	matrix[I, G] alpha_individual;
  	matrix[I, G] gamma_individual;
 //	vector[I] scale_individual;
	real scale_full;
	scale_full = fixed_scale == 0 ? scale[1] : 1.0;

	// individual level parameters
	for (rp in 1:RP) tau[rp] = 2.5 * tan(tau_unif[rp]);

  	L = diag_pre_multiply(tau, L_Omega);

	// individual level parameters
	beta = rep_matrix(mu', I) + (z * L);

	if (model_type == 1)
	  	alpha_individual = append_col(1 - exp(col(beta, RP_a)), rep_matrix(0, I, J));
	else if (model_type == 2)
	  	alpha_individual = 1 - exp(block(beta, 1, RP_a, I, G));
	else if (model_type == 3)
		alpha_individual = rep_matrix(1 - exp(col(beta, RP_a)), G);
	else
		alpha_individual = rep_matrix(0, I, G);

	if (model_type == 2)
	  gamma_individual = append_col(rep_vector(0, I), rep_matrix(1, I, J));
	else
	  gamma_individual = append_col(rep_vector(0, I), exp(block(beta, 1, RP_g, I, J)));

	psi_individual = block(beta, 1, 1, I, NPsi);

	for(t in 1:I){
		row_vector[J] util;
		util = psi_individual[task_individual[t]] * dat_psi[start[t]:end[t]]';
		lpsi[t] = util;
	}

	lpsi = lpsi + (block(alpha_individual, 1, 2, I, J) - 1) .* log(j_quant ./ block(gamma_individual, 1, 2, I, J) + 1) - log_price;
	f = (quant_full + gamma_individual) ./ (1 - alpha_individual);
	v = append_col((col(alpha_individual, 1) - 1) .* log_num, lpsi);
	v = exp(v / scale_full);
	sumv = v * ones_g;

	vf = nonzero .* v ./ f + (1 - nonzero);
	pf = (nonzero .* price_full .* f) * ones_g;

	for(i in 1:I){
		sumv[i] = pow(sumv[i], M[i]) * pow(scale_full, M[i] - 1);
		prodvf[i] = prod(vf[i]);
	}

	log_like = log((prodvf .* pf .* M_factorial) ./ sumv) .* weights;
	}
}

model {
  // priors on the parameters
	to_vector(z) ~ normal(0, 1);
	to_vector(mu) ~ normal(0, 1);
	L_Omega ~ lkj_corr_cholesky(4);                 // lkj prior
	scale ~ normal(1, 1);

  target += sum(log_like);//objective to target
}

generated quantities{
	real sum_log_lik = 0;
	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i];
	}
}

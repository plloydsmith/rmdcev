// saved as mdcev.stan
data {
	int K; // number of mixtures
	int I; // number of Individuals
	int J; // number of non-numeraire alternatives
	int NPsi; // number of psi covariates
	int L; // number of predictors in membership equation
	matrix[I * J, NPsi] dat_psi; // alt characteristics
	matrix[I, J] j_price; // non-numeraire price
	matrix[I, J] j_quant; // non-numeraire consumption
	vector[L] data_class[I];   // predictors for component membership
	vector[I] income;
	vector[I] num_price; // numeraire price
	vector[I] M_factorial;
	real prior_psi_sd;
	real prior_gamma_sd;
	real prior_alpha_sd;
	real prior_scale_sd;
	real prior_beta_m_sd;
	int<lower = 1, upper = 4> model_num; // 1 is les, 2 is alpha, 3 gamma (one alpha for all), 4 alpha's set to 1e-6
	int<lower=0, upper=1> fixed_scale; // indicator to fix scale
	int<lower=0, upper=1> trunc_data; //indicator to correct estimation for truncated data
    int<lower=0, upper=1> no_priors; //indicator to include priors or not
	vector[I] weights; // user supplied weights
}

transformed data {
	int G = J + 1;
	int A;
	vector[G] ones_g = rep_vector(1, G);
	matrix[I, G] price_full = append_col(num_price, j_price);
	matrix[I, G] quant_full;
	matrix[I, J] log_price = log(j_price);
	vector[I] log_inc = log(income);
  	vector[I] num_quant; // numeraire consumption
	vector[I] log_num;
	matrix[I, G] nonzero = rep_matrix(rep_vector(1, G)',I);
	vector[I] M;	//  Number of consumed goods (including numeraire)

	for(i in 1:I){
		num_quant[i] = (income[i] - j_price[i] * j_quant[i]') / num_price[i];
		for(g in 2:G){
			nonzero[i,g] = j_quant[i,g - 1] > 0 ? 1 : 0;
		}
  		M[i] = sum(nonzero[i]);
	}

	log_num = log(num_quant ./ num_price);
	quant_full = append_col(num_quant, j_quant);

 	if (model_num == 1 || model_num == 3)
 		A = 1;
	else if (model_num == 2)
	 	A = G;
	else if (model_num == 4)
		A = 0;
}

parameters {
	vector[NPsi] psi[K];
	vector<lower=0 >[model_num == 2 ? 0 : J] gamma[K];
	vector<lower=0, upper=1>[A] alpha[K];
	vector<lower=0>[fixed_scale == 0 ? K : 0] scale;
//	vector<lower=0>[K] scale;
  	matrix[K - 1, L] beta_m;  // mixture regression coeffs
}

transformed parameters {
	vector[I] log_like[K];
	vector[I] log_like_all;

	for (k in 1:K){
		matrix[I, J] lpsi = to_matrix(dat_psi * psi[k], I, J, 0);
		matrix[I, G] f;
		matrix[I, G] v;
		matrix[I, J] v_j;
		matrix[I, G] vf;
		vector[I] sumv;
		vector[I] pf;
		vector[I] prodvf;
		vector[G] gamma_full;
		vector[G] alpha_full;
		real scale_full;
		scale_full = fixed_scale == 0 ? scale[k] : 1.0;

		if (model_num == 1)
	  		alpha_full = append_row(alpha[k], rep_vector(0, J));
		else if (model_num == 2)
	  		alpha_full = alpha[k];
		else if (model_num == 3)
	  		alpha_full = rep_vector(alpha[k, 1], G);
		else
	  		alpha_full = rep_vector(1e-6, G);

		if (model_num == 2)
	  		gamma_full = append_row(0, rep_vector(1, J));
		else
	  		gamma_full = append_row(0, gamma[k]);

		v_j = lpsi + rep_matrix(alpha_full[2:G]'- 1, I) .* log(j_quant ./ rep_matrix(gamma_full[2:G]', I) + 1) - log_price;
		f = (quant_full + rep_matrix(gamma_full', I)) ./ rep_matrix((1 - alpha_full)', I);
		v = append_col((alpha_full[1] - 1) * log_num, v_j);
		v = exp(v / scale_full);
		sumv = v * ones_g;

		vf = nonzero .* v ./ f + (1 - nonzero);
		pf = (nonzero .* price_full .* f) * ones_g;

		for(i in 1:I){
			sumv[i] = pow(sumv[i], M[i]) * pow(scale_full, M[i] - 1);
			prodvf[i] = prod(vf[i]);
		}

		if (trunc_data == 1){
			matrix[I, G] v_1;
			vector[I] like_cond;
			vector[I] like_trunc;
			like_cond = prodvf .* pf .* M_factorial ./ sumv;

			v_1 = append_col((alpha_full[1] - 1) * log_inc, lpsi - log_price);
			v_1 = exp(v_1 / scale_full);
			sumv = v_1 * ones_g;

			like_trunc = col(v_1, 1) ./ sumv;

			for(i in 1:I)
				like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

			log_like[k] = log(like_cond ./ (1 - like_trunc)) .* weights;

		} else if (trunc_data == 0){
			log_like[k] = log((prodvf .* pf .* M_factorial) ./ sumv) .* weights;
		}
	}

	for(i in 1:I){
		vector[K] ltheta = log_softmax(append_row(0, beta_m * data_class[i])); // class membership equation
		vector[K] lps;
		for (k in 1:K){
			lps[k] = ltheta[k] + log_like[k,i];
		}
		log_like_all[i] = log_sum_exp(lps);
	}
}

model {
	if(no_priors == 0){
		scale ~ normal(1, prior_scale_sd);
	//	theta ~ dirichlet(rep_vector(2.0, K)); // no predictors
		to_vector(beta_m) ~ normal(0, prior_beta_m_sd);
		for (k in 1:K){
			to_vector(gamma[k]) ~ normal(0, prior_gamma_sd);
			to_vector(psi[k]) ~ normal(0, prior_psi_sd);
			to_vector(alpha[k]) ~ normal(.5, prior_alpha_sd);
		}
	}
  target += sum(log_like_all);

}

generated quantities{
	real sum_log_lik = 0;
	vector[I] theta[K];
	for(i in 1:I){
  		vector[K] theta1 = log_softmax(append_row(0, beta_m * data_class[i]));
		sum_log_lik = sum_log_lik + log_like_all[i];
		for(k in 1:K)
  			theta[k,i] =theta1[k];
	}
}

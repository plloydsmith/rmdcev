// mdcev model

functions {
#include /common/mdcev_ll.stan
}

data {
// declares I J NPsi K dat_psi j_price j_quant income num_price M_factorial
// prior_psi_sd prior_gamma_sd prior_alpha_sd prior_scale_sd
// model_num fixed_scale trunc_data no priors weights
#include /common/mdcev_data.stan

// additional LC data that can be set to 0 if not used
	int L; // number of predictors in membership equation
	vector[L] data_class[I];   // predictors for component membership
	real prior_beta_m_sd;

}

transformed data {
#include /common/mdcev_tdata.stan
}

parameters {
	vector[NPsi] psi[K];
	vector<lower=0 >[model_num == 2 ? 0 : J] gamma[K];
	vector<lower=0, upper=1>[A] alpha[K];
	vector<lower=0>[fixed_scale == 0 ? K : 0] scale;
  	matrix[K - 1, L] beta_m;  // mixture regression coeffs

//	vector[NPsi] psi;
//	vector<lower=0>[model_num == 2 ? 0 : J] gamma;
//	vector<lower=0, upper=1>[A] alpha;
//	vector<lower=0>[fixed_scale == 0 ? 1 : 0] scale;
}

transformed parameters {
	vector[I] log_like;

  	if(K == 1){
		vector[G] gamma_full;
		vector[G] alpha_full;
		real scale_full;
		scale_full = fixed_scale == 0 ? scale[1] : 1.0;

		if (model_num == 1)
		  alpha_full = append_row(alpha[1], rep_vector(0, J));
		else if (model_num == 2)
		  alpha_full = alpha[1];
		else if (model_num == 3)
		  alpha_full = rep_vector(alpha[1,1], G);
		else
		  alpha_full = rep_vector(1e-06, G);

		if (model_num == 2)
		  gamma_full = append_row(0, rep_vector(1, J));
		else
		  gamma_full = append_row(0, gamma[1]);

		log_like = mdcev_ll(j_quant, quant_full, log_price, log_num, price_full,
				log_inc, dat_psi, M, M_factorial, weights, // data
				psi[1], gamma_full, alpha_full, scale_full, 						// parameters
				I, J, G, ones_g, nonzero, model_num, fixed_scale, trunc_data);

	} else if (K > 1){
		vector[I] log_like_util[K];
		for (k in 1:K){
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

			log_like_util[k] = mdcev_ll(j_quant, quant_full, log_price, log_num, price_full,
				log_inc, dat_psi, M, M_factorial, weights, // data
				psi[k], gamma_full, alpha_full, scale_full, 						// parameters
				I, J, G, ones_g, nonzero, model_num, fixed_scale, trunc_data);
		}
		for(i in 1:I){
			vector[K] ltheta = log_softmax(append_row(0, beta_m * data_class[i])); // class membership equation
			vector[K] lps;
			for (k in 1:K){
				lps[k] = ltheta[k] + log_like_util[k,i];
			}
			log_like[i] = log_sum_exp(lps);
		}
	}
}

model {
  // priors on the parameters
  if(no_priors == 0){ // no priors ensure optimizing returns MLE
		scale ~ normal(1, prior_scale_sd);

  	if(K == 1){
	  psi[1] ~ normal(0, prior_psi_sd);
	  gamma[1] ~ normal(0, prior_gamma_sd);
	  alpha[1] ~ normal(.5, prior_alpha_sd);
	} else if (K > 1){
		to_vector(beta_m) ~ normal(0, prior_beta_m_sd);
		for (k in 1:K){
			to_vector(psi[k]) ~ normal(0, prior_psi_sd);
			to_vector(gamma[k]) ~ normal(0, prior_gamma_sd);
			to_vector(alpha[k]) ~ normal(.5, prior_alpha_sd);
		}
	}
  }

target += sum(log_like);//objective to target
}

generated quantities{
	real sum_log_lik = 0;
	vector[K > 1 ? I : 0] theta[K];

	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i];
		if (K > 1){
  			vector[K] theta1 = log_softmax(append_row(0, beta_m * data_class[i]));
			for(k in 1:K)
  				theta[k,i] = theta1[k];
		}
	}
}

// mdcev model

functions {
#include /common/mdcev_ll.stan
}

data {
// declares I J NPsi K dat_psi price_j quant_j income num_price M_factorial
// prior_psi_sd prior_gamma_sd prior_alpha_sd prior_scale_sd
// model_num fixed_scale1 trunc_data flat_priors weights
#include /common/mdcev_data.stan
	int K; // number of mixtures

// additional LC data that can be set to 0 if not used
	int L; // number of predictors in membership equation
	vector[L] data_class[I];   // predictors for component membership
	real prior_delta_sd;
}

transformed data {
#include /common/mdcev_tdata.stan
}

parameters {
	vector[NPsi] psi[K];
	vector<lower=0 >[Gamma] gamma[K];
	vector<lower=0, upper=1>[A] alpha[K];
	vector<lower=0>[fixed_scale1 == 0 ? K : 0] scale;
  	matrix[K - 1, L] delta;  // mixture regression coeffs
}

transformed parameters {
	vector[I] log_like;

  	if(K == 1){
  		matrix[I, J] lpsi = to_matrix(dat_psi[] * psi[1], I, J, 0);
		matrix[I, J] gamma_j = gamma_ll(gamma[1], I, J, model_num);
		matrix[I, J+1] alpha_full = alpha_ll(alpha[1], I, J, model_num);
		real scale_full = fixed_scale1 == 0 ? scale[1] : 1.0;

		log_like = mdcev_ll(quant_j, price_j, log_num, log_inc, M, log_M_fact, // data
				lpsi, gamma_j, col(alpha_full, 1), block(alpha_full, 1, 2, I, J), scale_full, // parameters
				I, J, nonzero, trunc_data);

	} else if (K > 1){
		vector[I] log_like_util[K];
		for (k in 1:K){
				matrix[I, J] lpsi = to_matrix(dat_psi[] * psi[k], I, J, 0);
				matrix[I, J] gamma_full = gamma_ll(gamma[k], I, J, model_num);
				matrix[I, J+1] alpha_full = alpha_ll(alpha[k], I, J, model_num);
				real scale_full = fixed_scale1 == 0 ? scale[k] : 1.0;

			log_like_util[k] = mdcev_ll(quant_j, price_j, log_num, log_inc, M, log_M_fact, // data
				lpsi, gamma_full, col(alpha_full, 1), block(alpha_full, 1, 2, I, J), scale_full, 						// parameters
				I, J, nonzero, trunc_data);
		}
		for(i in 1:I){
			vector[K] ltheta = log_softmax(append_row(0, delta * data_class[i])); // class membership equation
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
  if(flat_priors == 0){ // no priors ensure optimizing returns MLE
	scale ~ normal(1, prior_scale_sd);

  	if(K == 1){
	  psi[1] ~ normal(0, prior_psi_sd);
	  gamma[1] ~ normal(0, prior_gamma_sd);
	  alpha[1] ~ normal(.5, prior_alpha_sd);
	} else if (K > 1){
		to_vector(delta) ~ normal(0, prior_delta_sd);
		for (k in 1:K){
			to_vector(psi[k]) ~ normal(0, prior_psi_sd);
			to_vector(gamma[k]) ~ normal(1, prior_gamma_sd);
			to_vector(alpha[k]) ~ normal(.5, prior_alpha_sd);
		}
	}
  }

target += sum(log_like .* weights);//objective to target
}

generated quantities{
	real sum_log_lik = 0;
	vector[K > 1 ? I : 0] theta[K];

	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i];
		if (K > 1){
  			vector[K] theta1 = log_softmax(append_row(0, delta * data_class[i]));
			for(k in 1:K)
  				theta[k,i] = theta1[k];
		}
	}
}

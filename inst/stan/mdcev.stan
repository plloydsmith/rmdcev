// mdcev model

functions {
#include common/mdcev_ll.stan
}

data {
// declares I J NPsi K dat_psi price_j quant_j income num_price M_factorial
// prior_psi_sd prior_gamma_sd prior_alpha_sd prior_scale_sd
// model_num fixed_scale1 trunc_data flat_priors weights
#include common/mdcev_data.stan
	int K; // number of mixtures

// additional LC data that can be set to 0 if not used
	int L; // number of predictors in membership equation
//	array[I] vector[L] data_class;   // predictors for component membership
	matrix[I, L] data_class;   // predictors for component membership
	int<lower=0, upper=1> single_scale; // indicator to estimate one scale for lc
	real prior_delta_sd;
}

transformed data {
#include common/mdcev_tdata.stan

if (single_scale == 0 && fixed_scale1 == 0)
	S = K;
}

parameters {
	array[K] vector[NPsi] psi;
	array[K] vector[NPhi] phi;
	array[K] vector<lower=0 >[Gamma] gamma;
	array[K] vector<lower=0, upper=1>[A] alpha;
	vector<lower=0>[S] scale;
  	matrix[K - 1, L] delta;  // mixture regression coeffs
}

transformed parameters {
	vector[I] log_like;
	{
	array[K] vector[I] log_like_util;
	for (k in 1:K){
		matrix[I, J] gamma_full = gamma_ll(gamma[k], I, J, Gamma);
		real scale_full;
		matrix[I, J] lpsi;
		vector[I] alpha_1 = alpha_1_ll(alpha[k], I, model_num);
		vector[NPsi] psi_k = psi[k];

		if (S == 1)
  			scale_full = scale[1];
  		else if (S == K)
   			scale_full = scale[k];
  		else
  			scale_full = 1.0;

		if (psi_ascs == 1){
			lpsi = rep_matrix(append_row(0, head(psi_k, J-1))', I); //  alternative specific constants
			if (NPsi_ij > 0)
    			lpsi = lpsi + to_matrix(dat_psi[] * segment(psi_k, J, NPsi_ij), I, J, 0);
		} else if (psi_ascs == 0){
			if (NPsi_ij > 0)
    			lpsi = to_matrix(dat_psi[] * psi_k, I, J, 0);
			else if (NPsi_ij == 0)
    			lpsi = rep_matrix(0, I, J);
		}

		if(model_num < 5){
			matrix[I, J] alpha_j = alpha_j_ll(alpha[k], I, J, model_num);

			if(K == 1) {
			log_like = mdcev_ll(quant_j, price_j, log_num, income, M, log_M_fact, // data
				lpsi, gamma_full, alpha_1, alpha_j, scale_full, 						// parameters
				I, J, nonzero, trunc_data);
			} else if (K > 1){
			log_like_util[k] = mdcev_ll(quant_j, price_j, log_num, income, M, log_M_fact, // data
				lpsi, gamma_full, alpha_1, alpha_j, scale_full, 						// parameters
				I, J, nonzero, trunc_data);
			}
		} else if (model_num == 5){
			matrix[I, J] phi_ij;
			if(NPhi > 0)
	    		phi_ij = exp(to_matrix(dat_phi[] * phi[k], I, J, 0));
			else if(NPhi == 0)
	    		phi_ij = rep_matrix(1, I, J);

			if(K == 1) {
			log_like = kt_ll(quant_j, price_j, log_num, income,
  				lpsi, phi_ij, gamma_full, alpha_1, scale_full,
  				I, J, nonzero, trunc_data, jacobian_analytical_grad);
			} else if (K > 1){
			log_like_util[k] = kt_ll(quant_j, price_j, log_num, income,
  				lpsi, phi_ij, gamma_full, alpha_1, scale_full,
  				I, J, nonzero, trunc_data, jacobian_analytical_grad);
			}
		}
	}

	if (K > 1){
//		matrix[I, K] theta = append_col(rep_vector(0, I), data_class * delta'); // class membership equation
		for(i in 1:I){
			//vector[K] theta = softmax(append_row(0, delta * data_class[i]')); // class membership equation
//			vector[K] theta_temp = softmax(theta[i]); // class membership equation
		//	log_like[i] = log_mix(theta, log_like_util[i]);
			vector[K] ltheta = log_softmax(append_row(0, delta * data_class[i]')); // class membership equation
			vector[K] lps;
			for (k in 1:K){
				lps[k] = ltheta[k] + log_like_util[k,i];
			}
			log_like[i] = log_sum_exp(lps);
		}
	}
	}
}

model {
  // priors on the parameters
  if(flat_priors == 0){ // no priors ensure optimizing returns MLE
	scale ~ normal(0, prior_scale_sd);

  	if(K == 1){
	  psi[1] ~ normal(0, prior_psi_sd);
	  phi[1] ~ normal(0, prior_phi_sd);
	  gamma[1] ~ normal(1, prior_gamma_sd);
	  alpha[1] ~ beta(prior_alpha_shape, prior_alpha_shape);

	} else if (K > 1){
		to_vector(delta) ~ normal(0, prior_delta_sd);
		for (k in 1:K){
			to_vector(psi[k]) ~ normal(0, prior_psi_sd);
			to_vector(phi[k]) ~ normal(0, prior_phi_sd);
			to_vector(gamma[k]) ~ normal(1, prior_gamma_sd);
			to_vector(alpha[k]) ~ beta(prior_alpha_shape, prior_alpha_shape);
		}
	}
  }

target += dot_product(log_like, weights);//objective to target
}

generated quantities{
	real sum_log_lik = 0;
	matrix[K > 1 ? I : 0, K] theta;

	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i] * weights[i];
		if (K > 1)
  			theta[i] = softmax(append_row(0, delta * data_class[i]'))';
	}
}

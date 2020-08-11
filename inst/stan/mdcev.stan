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
	vector<lower=0 >[NPhi] phi[K];
	vector<lower=0 >[Gamma] gamma[K];
	vector<lower=0, upper=1>[A] alpha[K];
	vector<lower=0>[fixed_scale1 == 0 ? K : 0] scale;
  	matrix[K - 1, L] delta;  // mixture regression coeffs
}

transformed parameters {
	vector[I] log_like;
	{
	vector[I] log_like_util[K];
	for (k in 1:K){
		matrix[I, J] gamma_full = gamma_ll(gamma[k], I, J, Gamma);
		real scale_full = fixed_scale1 == 0 ? scale[k] : 1.0;
		matrix[I, J] lpsi;
		vector[NPsi] psi_k = psi[k];
		if (psi_ascs == 1){
			lpsi = rep_matrix(append_row(0, segment(psi_k, 1, J-1))', I); //  alternative specific constants
			if (NPsi_ij > 0)
    			lpsi = lpsi + to_matrix(dat_psi[] * segment(psi_k, J, NPsi_ij), I, J, 0);
		} else if (psi_ascs == 0){
			if (NPsi_ij > 0)
    			lpsi = to_matrix(dat_psi[] * psi_k, I, J, 0);
			else if (NPsi_ij == 0)
    			lpsi = rep_matrix(0, I, J);
		}

		if(model_num < 5){
			matrix[I, J+1] alpha_full = alpha_ll(alpha[k], I, J, model_num);

			if(K == 1) {
			log_like = mdcev_ll(quant_j, price_j, log_num, income, M, log_M_fact, // data
				lpsi, gamma_full, col(alpha_full, 1), block(alpha_full, 1, 2, I, J), scale_full, 						// parameters
				I, J, nonzero, trunc_data);
			} else if (K > 1){
			log_like_util[k] = mdcev_ll(quant_j, price_j, log_num, income, M, log_M_fact, // data
				lpsi, gamma_full, col(alpha_full, 1), block(alpha_full, 1, 2, I, J), scale_full, 						// parameters
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
  				lpsi, phi_ij, gamma_full, rep_vector(alpha[k, 1], I), scale_full,
  				I, J, nonzero, trunc_data, jacobian_analytical_grad);
			} else if (K > 1){
			log_like_util[k] = kt_ll(quant_j, price_j, log_num, income,
  				lpsi, phi_ij, gamma_full, rep_vector(alpha[k, 1], I), scale_full,
  				I, J, nonzero, trunc_data, jacobian_analytical_grad);
			}
		}
	}

	if (K > 1){
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
}

model {
  // priors on the parameters
  if(flat_priors == 0){ // no priors ensure optimizing returns MLE
	scale ~ normal(1, prior_scale_sd);

  	if(K == 1){
	  psi[1] ~ normal(0, prior_psi_sd);
	  phi[1] ~ normal(0, prior_phi_sd);
	  gamma[1] ~ normal(1, prior_gamma_sd);
	  alpha[1] ~ normal(.5, prior_alpha_sd);
	} else if (K > 1){
		to_vector(delta) ~ normal(0, prior_delta_sd);
		for (k in 1:K){
			to_vector(psi[k]) ~ normal(0, prior_psi_sd);
			to_vector(phi[k]) ~ normal(0, prior_phi_sd);
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
		sum_log_lik = sum_log_lik + log_like[i] * weights[i];
		if (K > 1){
  			vector[K] theta1 = log_softmax(append_row(0, delta * data_class[i]));
			for(k in 1:K)
  				theta[k,i] = theta1[k];
		}
	}
}

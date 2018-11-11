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
	vector[L] dat_class[I];   // predictors for component membership
	vector[I] num_price; // psi variables
	vector[I] num_quant; // psi variables
	int<lower=0, upper=1> print_ll; //indicator to print log_lik at each iteration
  	int<lower = 1, upper = 4> model_type; // 1 is les, 2 is alpha, 3 gamma (one alpha for all), 4 alpha's set to 0
  	vector[I] weights; // user supplied weights
  	vector[I] M_factorial;
	matrix[I, J + 1] nonzero;
	vector[I] M;	//  Number of consumed goods (including numeraire)
}

transformed data {
	int G = J + 1;
	int A;
	vector[G] ones_g = rep_vector(1, G);
	matrix[I, G] price_full = append_col(num_price, j_price);
	matrix[I, G] quant_full = append_col(num_quant, j_quant);
	matrix[I, J] log_price = log(j_price);
	vector[I] log_num = log(num_quant ./ num_price);

 	if (model_type == 1 || model_type == 3)
 		A = 1;
	else if (model_type == 2)
 		A = G;
	else if (model_type == 4)
		A = 0;
}

parameters {
	vector[NPsi] psi[K];
	vector<lower=0 >[model_type == 2 ? 0 : J] gamma;
	vector<lower=0, upper=1>[A] alpha[K];
	vector<lower=0>[K] scale;
//	simplex[K] theta;  // membership coefficients (no predictors)
  matrix[K - 1, L] beta_m;  // mixture regression coeffs
}

transformed parameters {
//	real scale = 1;
//	vector[K] ltheta = log(theta); // cache log calculation
	vector[I] log_like[K];
	vector[I] log_like_all;

	for (k in 1:K){
		matrix[I, J] lpsi = to_matrix(dat_psi[] * psi[k], I, J, 0);
		matrix[I, G] f;
		matrix[I, G] v;
		matrix[I, G] vf;
		vector[I] sumv;
		vector[I] pf;
		vector[I] prodvf;
		vector[G] gamma_full;
		vector[G] alpha_full;

		if (model_type == 1)
	  		alpha_full = append_row(alpha[k], rep_vector(0, J));
		else if (model_type == 2)
	  		alpha_full = alpha[k];
		else if (model_type == 3)
	  		alpha_full = rep_vector(alpha[k,1], G);
		else
	  		alpha_full = rep_vector(0, G);

		if (model_type == 2)
	  		gamma_full = append_row(0, rep_vector(1, J));
		else
	  		gamma_full = append_row(0, gamma);

		lpsi = lpsi + rep_matrix(alpha_full[2:G]'- 1, I) .* log(j_quant ./ rep_matrix(gamma_full[2:G]', I) + 1) - log_price;
		f = (quant_full + rep_matrix(gamma_full', I)) ./ rep_matrix((1 - alpha_full)', I);
		v = append_col((alpha_full[1] - 1) * log_num, lpsi);
		v = exp(v / scale[k]);
		sumv = v * ones_g;

		vf = nonzero .* v ./ f + (1 - nonzero);
		pf = (nonzero .* price_full .* f) * ones_g;

		for(i in 1:I){
			sumv[i] = pow(sumv[i], M[i]) * pow(scale[k], M[i] - 1);
			prodvf[i] = prod(vf[i]);
		}

		log_like[k] = log((prodvf .* pf .* M_factorial) ./ sumv) .* weights;
	}

	for(i in 1:I){
		vector[K] ltheta = log_softmax(append_row(0, beta_m * dat_class[i])); // class membership equation
		vector[K] lps;
		for (k in 1:K){
			lps[k] = ltheta[k] + log_like[k,i];
		}
		log_like_all[i] = log_sum_exp(lps);
	}
}

model {
	gamma ~ normal(0, 10);
	scale ~ normal(1, 1);
//	theta ~ dirichlet(rep_vector(2.0, K)); // no predictors
	to_vector(beta_m) ~ normal(0,10);
	for (k in 1:K){
		to_vector(psi[k]) ~ normal(0, 10);
		to_vector(alpha[k]) ~ normal(.5 ,.5);
	}
  target += sum(log_like_all);

}

generated quantities{
	real sum_log_lik = 0;
	vector[I] theta[K];
	for(i in 1:I){
  		vector[K] theta1 = log_softmax(append_row(0, beta_m * dat_class[i]));
		sum_log_lik = sum_log_lik + log_like_all[i];
		for(k in 1:K)
  			theta[k,i] =theta1[k];
	}
	if(print_ll==1){
		print("sum_log_lik =", sum_log_lik);
	}
}

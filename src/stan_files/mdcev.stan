// saved as mdcev.stan
data {
  int I; // number of Individuals
  int J; // number of non-numeraire alternatives
  int NPsi; // number of alt covariates
  matrix[I * J, NPsi] dat_psi; // alt characteristics
  matrix[I, J] j_price; // non-numeraire price
  matrix[I, J] j_quant; // non-numeraire consumption
  vector[I] income;
  vector[I] num_price; // numeraire price
  vector[I] M_factorial; // (M-1)!
  real prior_psi_sd;
  real prior_gamma_sd;
  real prior_alpha_sd;
  real prior_scale_sd;
  int<lower = 1, upper = 4> model_num; // 1 is les, 2 is alpha, 3 gamma (one alpha for all), 4 alpha's set to 1e-06
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
	vector[NPsi] psi;
	vector<lower=0>[model_num == 2 ? 0 : J] gamma;
	vector<lower=0, upper=1>[A] alpha;
	vector<lower=0>[fixed_scale == 0 ? 1 : 0] scale;
}

transformed parameters {
	vector[I] log_like;
	{
	matrix[I, J] lpsi = to_matrix(dat_psi[] * psi, I, J, 0);
	matrix[I, J] v_j;
	matrix[I, G] f;
	matrix[I, G] v;
	matrix[I, G] vf;
	vector[I] sumv;
	vector[I] pf;
	vector[I] prodvf;
	vector[G] gamma_full;
	vector[G] alpha_full;
	real scale_full;
	scale_full = fixed_scale == 0 ? scale[1] : 1.0;

	if (model_num == 1)
	  alpha_full = append_row(alpha, rep_vector(0, J));
	else if (model_num == 2)
	  alpha_full = alpha;
	else if (model_num == 3)
	  alpha_full = rep_vector(alpha[1], G);
	else
	  alpha_full = rep_vector(1e-06, G);

	if (model_num == 2)
	  gamma_full = append_row(0, rep_vector(1, J));
	else
	  gamma_full = append_row(0, gamma);

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

		log_like = log(like_cond ./ (1 - like_trunc)) .* weights;

	} else if (trunc_data == 0){
		log_like = log((prodvf .* pf .* M_factorial) ./ sumv) .* weights;
	}
	}
}

model {
  // priors on the parameters
  if(no_priors == 0){
	  psi ~ normal(0, prior_psi_sd);
	  gamma ~ normal(0, prior_gamma_sd);
	  alpha ~ normal(.5, prior_alpha_sd);
	  scale ~ normal(1, prior_scale_sd);
	}

  target += sum(log_like);//objective to target
}

generated quantities{
	real sum_log_lik = 0;
	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i];
	}
}

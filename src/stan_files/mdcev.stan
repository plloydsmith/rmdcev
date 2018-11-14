// saved as mdcev_fixed.stan
data {
  int I; // number of Individuals
  int J; // number of non-numeraire alternatives
  int NPsi; // number of alt covariates
  matrix[I * J, NPsi] dat_psi; // alt characteristics
  matrix[I, J] j_price; // non-numeraire price
  matrix[I, J] j_quant; // non-numeraire consumption
  vector[I] income;
  vector[I] num_price; // numeraire price
  vector[I] num_quant; // numeraire consumption
  int<lower=0, upper=1> print_ll; //indicator to print log_lik at each iteration.
  int<lower = 1, upper = 4> model_type; // 1 is les, 2 is alpha, 3 gamma (one alpha for all), 4 alpha's set to 0
  int<lower=0, upper=1> fixed_scale; // indicator to fix scale
  vector[I] weights; // user supplied weights
	vector[I] M_factorial;
	matrix[I, J + 1] nonzero;
	vector[I] M;	//  Number of consumed goods (including numeraire)
  int<lower=0, upper=1> trunc_data; //indicator to correct estimation for truncated data
}

transformed data {
	int G = J + 1;
	int A;
	vector[G] ones_g = rep_vector(1, G);
	matrix[I, G] price_full = append_col(num_price, j_price);
	matrix[I, G] quant_full = append_col(num_quant, j_quant);
	matrix[I, J] log_price = log(j_price);
	vector[I] log_num = log(num_quant ./ num_price);
	vector[I] log_inc = log(income);

 	if (model_type == 1 || model_type == 3)
 		A = 1;
	else if (model_type == 2)
	 	A = G;
	else if (model_type == 4)
		A = 0;
}

parameters {
	vector[NPsi] psi;
	vector<lower=0>[model_type == 2 ? 0 : J] gamma;
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
	matrix[I, G] v_m;
	vector[I] sumv;
	vector[I] pf;
	vector[I] prodvf;
	vector[I] prodv;
	vector[G] gamma_full;
	vector[G] alpha_full;
	real scale_full;
	scale_full = fixed_scale == 0 ? scale[1] : 1.0;

	if (model_type == 1)
	  alpha_full = append_row(alpha, rep_vector(0, J));
	else if (model_type == 2)
	  alpha_full = alpha;
	else if (model_type == 3)
	  alpha_full = rep_vector(alpha[1], G);
	else
	  alpha_full = rep_vector(0, G);

	if (model_type == 2)
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
		vector[I] like_temp;
		vector[I] like_trunc;
		like_temp = prodvf .* pf .* M_factorial ./ sumv;

		v = append_col((alpha_full[1] - 1) * log_inc, lpsi - log_price);
		v = exp(v / scale_full);
		sumv = v * ones_g;

		v_m = nonzero .* v + (1 - nonzero);

		for(i in 1:I){
			sumv[i] = pow(sumv[i], M[i]) * pow(scale_full, M[i] - 1);
			prodv[i] = prod(v_m[i]);
		}

		like_trunc = prodv ./ sumv .* M_factorial;

		for(i in 1:I)
			like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 0;

		log_like = log(like_temp ./ ( 1- like_trunc)) .* weights;

	} else
		log_like = log((prodvf .* pf .* M_factorial) ./ sumv) .* weights;
	}
}

model {
  // priors on the parameters
  psi ~ normal(0, 10);
  gamma ~ normal(0, 10);
  alpha ~ normal(.5, .5);
  scale ~ normal(1, 1);

  target += sum(log_like);//objective to target
}

generated quantities{
	real sum_log_lik = 0;
	for(i in 1:I){
		sum_log_lik = sum_log_lik + log_like[i];
	}

	if(print_ll==1){
		print("sum_log_lik =", sum_log_lik);
	}
}

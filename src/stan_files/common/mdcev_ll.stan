
/* These functions (without the underscores) are all documented in R
   See also Appendix C of Pinheiro and Bates
https://books.google.com/books?id=3TVDAAAAQBAJ&lpg=PR3&dq=Pinheiro%20and%20Bates&pg=PA511#v=onepage&q&f=false
  These functions may be numerically unstable
*/

matrix gamma_ll(vector gamma, int I, int J, int model_num) {

	matrix[I, J] gamma_j;

	if (model_num == 2)
	  gamma_j = rep_matrix(1, I, J);
	else
	  gamma_j = rep_matrix(gamma', I);

return(gamma_j);
}

matrix alpha_ll(vector alpha, int I, int J, int model_num) {

	matrix[I, J+1] alpha_full;

	if (model_num == 1)
	  alpha_full = append_col(rep_vector(alpha[1], I), rep_matrix(0, I, J));
	else if (model_num == 2)
	  alpha_full = rep_matrix(alpha', I);
	else if (model_num == 3)
	  alpha_full = rep_matrix(alpha[1], I, J+1);
	else
	  alpha_full = rep_matrix(1e-06, I, J+1);

return(alpha_full);
}

vector mdcev_ll(matrix j_quant, matrix j_price, vector log_num, vector log_inc,
				vector M, vector log_M_fact, // data
				matrix lpsi, matrix gamma_j, vector alpha1, matrix alpha_j, real scale_full,// parameters
				int I, int J, matrix nonzero, int trunc_data)  {//options

	vector[I] log_like;
	vector[J] ones_j = rep_vector(1, J);
	matrix[I, J] v_j= lpsi + (alpha_j - 1) .* log(j_quant ./ gamma_j + 1) - log(j_price);
	vector[I] v1 = (alpha1 - 1) .* log_num / scale_full;
	matrix[I, J] logf = log(1 - alpha_j) - log(j_quant + gamma_j);
	vector[I] logf1 = log(1 - alpha1) - log_num;
	v_j = v_j / scale_full;

	if (trunc_data == 1){
		matrix[I, J+1] v_1;
		vector[I] like_cond;
		vector[I] like_trunc;
		vector[I] sumv;

		like_cond = exp((1 - M) * log(scale_full) + logf1 + v1 + (nonzero .* (logf + v_j)) * ones_j +
			log(inv(exp(logf1)) + (nonzero .* j_price ./ exp(logf)) * ones_j) -
			M .* log(exp(v1) + exp(v_j) * ones_j) + log_M_fact);

		v_1 = append_col((alpha1 - 1) .* log_inc, lpsi - log(j_price));
		v_1 = exp(v_1 / scale_full);
		sumv = v_1 * rep_vector(1, J+1);

		like_trunc = col(v_1, 1) ./ sumv;

		for(i in 1:I)
			like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

		log_like = log(like_cond ./ (1 - like_trunc));

	} else if (trunc_data == 0){
		log_like = (1 - M) * log(scale_full) + logf1 + v1 + (nonzero .* (logf + v_j)) * ones_j +
			log(inv(exp(logf1)) + (nonzero .* j_price ./ exp(logf)) * ones_j) -
			M .* log(exp(v1) + exp(v_j) * ones_j) + log_M_fact;
	}

return(log_like);
}

vector mdcev_ll_old(matrix j_quant, matrix j_price, vector log_num, vector log_inc,
				vector M, vector log_M_fact, // data
				matrix lpsi, matrix gamma_j, vector alpha1, matrix alpha_j, real scale_full,// parameters
				int I, int J, matrix nonzero, int trunc_data)  {//options

	vector[I] log_like;
	matrix[I, J] v_j;
	vector[I] f1;
	matrix[I, J] fj;
	matrix[I, J+1] f;
	matrix[I, J+1] v;
	matrix[I, J+1] vf;
	vector[I] sumv;
	vector[I] pf;
	vector[I] prodvf;
	vector[J+1] ones_g = rep_vector(1, J+1);
	matrix[I, J+1] price_full = append_col(rep_vector(1, I), j_price);
	matrix[I, J+1] nonzero_full = append_col(rep_vector(1, I), nonzero);

	v_j= lpsi + (alpha_j - 1) .* log(j_quant ./ gamma_j + 1) - log(j_price);
	f1 = exp(log_num) ./ (1 - alpha1);
	fj = (j_quant + gamma_j) ./ (1 - alpha_j);
	f = append_col(f1, fj);
	v = append_col((alpha1 - 1) .* log_num, v_j);
	v = exp(v / scale_full);
	sumv = v * ones_g;

	vf = nonzero_full .* v ./ f + (1 - nonzero_full);
	pf = (nonzero_full .* price_full .* f) * ones_g;

	for(i in 1:I){
		sumv[i] = pow(sumv[i], M[i]) * pow(scale_full, M[i] - 1);
		prodvf[i] = prod(vf[i]);
	}

	if (trunc_data == 1){
		matrix[I, J+1] v_1;
		vector[I] like_cond;
		vector[I] like_trunc;
		like_cond = prodvf .* pf .* exp(log_M_fact) ./ sumv;

		v_1 = append_col((alpha1 - 1) .* log_inc, lpsi - log(j_price));
		v_1 = exp(v_1 / scale_full);
		sumv = v_1 * ones_g;

		like_trunc = col(v_1, 1) ./ sumv;

		for(i in 1:I)
			like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

		log_like = log(like_cond ./ (1 - like_trunc));

	} else if (trunc_data == 0){
		log_like = log((prodvf .* pf .* exp(log_M_fact)) ./ sumv);
	}
return(log_like);
}

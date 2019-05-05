
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

vector mdcev_ll(matrix quant_j, matrix price_j, vector log_num, vector log_inc,
				vector M, vector log_M_fact, // data
				matrix lpsi, matrix gamma_j, vector alpha1, matrix alpha_j, real scale_full,// parameters
				int I, int J, matrix nonzero, int trunc_data)  {//options

	vector[I] log_like;
	vector[J] ones_j = rep_vector(1, J);
	matrix[I, J] v_j= lpsi + (alpha_j - 1) .* log(quant_j ./ gamma_j + 1) - log(price_j);
	vector[I] v1 = (alpha1 - 1) .* log_num / scale_full;
	matrix[I, J] logf = log(1 - alpha_j) - log(quant_j + gamma_j);
	vector[I] logf1 = log(1 - alpha1) - log_num;
	v_j = v_j / scale_full;

	if (trunc_data == 0){
		log_like = (1 - M) * log(scale_full) + logf1 + v1 + (nonzero .* (logf + v_j)) * ones_j +
			log(inv(exp(logf1)) + (nonzero .* price_j ./ exp(logf)) * ones_j) -
			M .* log(exp(v1) + exp(v_j) * ones_j) + log_M_fact;

	} else if (trunc_data == 1){
		matrix[I, J+1] v_1;
		vector[I] like_cond;
		vector[I] like_trunc;
		vector[I] sumv;

		like_cond = exp((1 - M) * log(scale_full) + logf1 + v1 + (nonzero .* (logf + v_j)) * ones_j +
			log(inv(exp(logf1)) + (nonzero .* price_j ./ exp(logf)) * ones_j) -
			M .* log(exp(v1) + exp(v_j) * ones_j) + log_M_fact);

		v_1 = append_col((alpha1 - 1) .* log_inc, lpsi - log(price_j));
		v_1 = exp(v_1 / scale_full);
		sumv = v_1 * rep_vector(1, J+1);

		like_trunc = col(v_1, 1) ./ sumv;

		for(i in 1:I)
			like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

		log_like = log(like_cond ./ (1 - like_trunc));
	}

return(log_like);
}

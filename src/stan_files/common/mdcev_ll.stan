
/* These functions (without the underscores) are all documented in R
   See also Appendix C of Pinheiro and Bates
https://books.google.com/books?id=3TVDAAAAQBAJ&lpg=PR3&dq=Pinheiro%20and%20Bates&pg=PA511#v=onepage&q&f=false
  These functions may be numerically unstable
*/

matrix gamma_ll(vector gamma, int I, int J, int G, int model_num) {

	matrix[I, G] gamma_full;

	if (model_num == 2)
	  gamma_full = append_col(rep_vector(0, I), rep_matrix(1, I, J));
	else
	  gamma_full = append_col(rep_vector(0, I), rep_matrix(gamma', I));

return(gamma_full);
}

matrix alpha_ll(vector alpha, int I, int J, int G, int model_num) {

	matrix[I, G] alpha_full;

	if (model_num == 1)
	  alpha_full = append_col(rep_vector(alpha[1], I), rep_matrix(0, I, J));
	else if (model_num == 2)
	  alpha_full = rep_matrix(alpha', I);
	else if (model_num == 3)
	  alpha_full = rep_matrix(alpha[1], I, G);
	else
	  alpha_full = rep_matrix(1e-06, I, G);

return(alpha_full);
}

vector mdcev_ll(matrix j_quant, matrix quant_full, matrix log_price,
				vector log_num, matrix price_full, vector log_inc,
				matrix dat_psi, vector M, vector M_factorial, vector weights, // data
				matrix lpsi, matrix gamma_full, matrix alpha_full, real scale_full,// parameters
				int I, int J, int G, vector ones_g, matrix nonzero,
				int model_num, int fixed_scale, int trunc_data)  {//options

	vector[I] log_like;
	matrix[I, J] v_j;
	matrix[I, G] f;
	matrix[I, G] v;
	matrix[I, G] vf;
	vector[I] sumv;
	vector[I] pf;
	vector[I] prodvf;

	v_j = lpsi + (block(alpha_full, 1, 2, I, J) - 1) .* log(j_quant ./ block(gamma_full, 1, 2, I, J) + 1) - log_price;
	f = (quant_full + gamma_full) ./ (1 - alpha_full);
	v = append_col((col(alpha_full, 1) - 1) .* log_num, v_j);
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

		v_1 = append_col((col(alpha_full, 1) - 1) .* log_inc, lpsi - log_price);
		v_1 = exp(v_1 / scale_full);
		sumv = v_1 * ones_g;

		like_trunc = col(v_1, 1) ./ sumv;

		for(i in 1:I)
			like_trunc[i] = like_trunc[i] < 1 ? like_trunc[i] : 1;

		log_like = log(like_cond ./ (1 - like_trunc)) .* weights;

	} else if (trunc_data == 0){
		log_like = log((prodvf .* pf .* M_factorial) ./ sumv) .* weights;
	}
return(log_like);
}

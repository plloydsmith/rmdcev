
/* This computes the LL for the MDCEV model
*/

vector mdcev_ll(matrix j_quant, matrix quant_full, matrix log_price, vector log_num, matrix price_full,
				vector log_inc, matrix dat_psi, vector M, vector M_factorial, vector weights, // data
				vector psi, vector gamma_full, vector alpha_full, real scale_full, 						// parameters
				int I, int J, int G, vector ones_g, matrix nonzero, int model_num, int fixed_scale, int trunc_data)  {//options

	vector[I] log_like;
	matrix[I, J] lpsi = to_matrix(dat_psi[] * psi, I, J, 0);
	matrix[I, J] v_j;
	matrix[I, G] f;
	matrix[I, G] v;
	matrix[I, G] vf;
	vector[I] sumv;
	vector[I] pf;
	vector[I] prodvf;

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
return(log_like);
}

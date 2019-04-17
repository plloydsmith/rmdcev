
int G = J + 1;
int A;
int Gamma;
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

if (model_num == 1 || model_num == 3){
 	A = 1;
 	Gamma = J;
} else if (model_num == 2){
 	A = G;
 	Gamma = 0;
} else if (model_num == 4){
	A = 0;
 	Gamma = J;
}

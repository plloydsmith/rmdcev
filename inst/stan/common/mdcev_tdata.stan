
int A;
int Gamma;
int NPsi;
vector[I] log_num;
matrix[I, J] nonzero;
vector[I] M;	//  Number of consumed alts (including numeraire)
vector[I] log_M_fact;

for(i in 1:I){
	log_num[i] = log(income[i] - price_j[i] * quant_j[i]');
	for(j in 1:J){
		nonzero[i,j] = quant_j[i,j] > 0 ? 1 : 0;
	}
  	M[i] = sum(nonzero[i])+1; // add 1 for numeraire
}

log_M_fact = lgamma(M); // lgamma(M) = log((M-1)!)

// Set number of Psi parameters
if (psi_ascs == 1){
 	NPsi = J - 1 + NPsi_ij;
} else if (psi_ascs == 0){
 	NPsi = NPsi_ij;
}

if (model_num == 1 || model_num == 3 || model_num == 5){
 	A = 1;
 	Gamma = J;
} else if (model_num == 2){
 	A = J + 1;
 	Gamma = 0;
} else if (model_num == 4){
	A = 0;
 	Gamma = J;
}

if (model_num != 2 && gamma_ascs == 0){
 	Gamma = 1;
}

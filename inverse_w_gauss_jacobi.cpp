#include <iostream>
#include "gauss_jacobi.cpp"


std::vector<std::vector<double>> inversa_gauss_seidel(std::vector<std::vector<double>> A){
	int l = A.size();
	int c = A[0].size();
	std::vector<double> B(c);

	std::vector<std::vector<double>> InverseA(l);
	for (int j=0;j<l;j++){
		InverseA[j].reserve(c);
		std::vector<double> aux(c);
		for (int i=0;i<c;i++){
			aux[i] = 0;
		}
		InverseA[j] = aux;
	}
	std::cout << InverseA.size() << '\n' << InverseA[0].size() << '\n' ;
	std::cout << "vector Inverse Alocado" << '\n';
	for (int x = 0; x < l; x++ ){
		std::vector<double> result;
		for (int i = 0; i < c; i++){
			if (i == x){
				B[i] = 1;
			}
			else{
				B[i] = 0;
			}
		}
		result = gauss_jacobi(A, B, 0.00001, 100);
		std::vector<double> temp;
		for (int y = 0; y < l ; y++){
			InverseA[y][x] += result[y];
			
		}
	}
	std::cout << "vector Inverse Criado" << '\n';


	
	


	return (InverseA);
}
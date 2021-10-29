#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../include/Auxiliary_Functions.h"
#include <cmath>
#include <complex>
#include <chrono>

MatrixXcd Kronecker_Product(MatrixXcd A, MatrixXcd B){

	MatrixXcd C(A.rows() * B.rows(), A.cols() * B.rows());

	for (int i = 0; i < A.rows(); i++){
		for (int j = 0; j < A.cols(); j++){
			C.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j)*B;
		}
	}

	return C;
}

MatrixXcd Implementing_Chiral_Symmetry_W(MatrixXcd W_aux, int N, int ress, int _chiral_deg){
	
	MatrixXcd Waux_A(ress, 1); MatrixXcd Waux_B(ress, 1);
	Waux_A.setZero(); Waux_B.setZero();
	
	MatrixXcd Waux_new(_chiral_deg * ress, 2);
	Waux_new.setZero();

	MatrixXcd W_new(_chiral_deg * ress, _chiral_deg * N);
	W_new.setZero();

	for (int k=1; k < W_aux.cols(); k += 2){
		Waux_A = W_aux.block(0, k-1, ress, 1);
		Waux_B = W_aux.block(0, k, ress, 1);	
	
		Waux_new.block(0,0, ress, 1) = Waux_A; Waux_new.block(0, 1, ress, 1) = Waux_B;
		Waux_new.block(ress, 0, ress, 1) = -Waux_B; Waux_new.block(ress, 1, ress, 1) = Waux_A;

		W_new.block(0, (k-1), _chiral_deg * ress, 2) = Waux_new;
	}

	// Uncomment the code below to verify the Chiral Symmetry Constraint //

	// cout << "\nThe new Matrix W1: \n" << W1_new << endl;

	// MatrixXcd I_spin = MatrixXcd::Identity(2, 2);
	// MatrixXcd I_ress = MatrixXcd::Identity(ress, ress);
	// MatrixXcd I_N = MatrixXcd::Identity(N, N);
	// MatrixXcd paulimatrix_y(2,2);
       
	// paulimatrix_y << 0, complex_identity,
	//	 	-complex_identity, 0;

	// cout << "\nExpressao1: \n" << W_new << endl;
	// cout << "\nExpressao2: \n" << Kronecker_Product(paulimatrix_y, I_ress)*W1_new*Kronecker_Product(I_N, paulimatrix_y) << endl;
	// cout << "\nDifference of Symmetry equation: " << (W_new-Kronecker_Product(paulimatrix_y, I_ress)*W_new*Kronecker_Product(I_N, paulimatrix_y)).cwiseAbs() << endl; 
	
	return W_new;
}

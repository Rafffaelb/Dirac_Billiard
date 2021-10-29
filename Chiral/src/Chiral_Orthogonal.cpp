#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../include/Chiral_Orthogonal.h"
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <random>
#include <fstream>
#include <string>

using namespace std;

Chiral_Orthogonal::Chiral_Orthogonal(double lambda, int num_steps, int spin_deg, int chiral_deg){

	this -> _lambda = lambda;
	this -> _num_steps = num_steps;
	this -> _spin_deg = spin_deg;
	this -> _chiral_deg = chiral_deg;
}

Chiral_Orthogonal::~Chiral_Orthogonal() {}

MatrixXcd Kronecker_Product(MatrixXcd A, MatrixXcd B);
MatrixXcd Implementing_Chiral_Symmetry_W(MatrixXcd W_aux, int N, int ress, int _chiral_deg);

void Chiral_Orthogonal::Create_W(MatrixXcd* W_pointer, int ress, int N1, int N2, double lambda, double y){

	complex<double> complex_identity(0,1);

	MatrixXcd W1_aux(ress, _chiral_deg * N1); 
	W1_aux.setZero();

	for (int j=1; j < ress+1; j++ ){
		for (int k=1; k < 2*N1+1; k++){
			std::complex<double> aux(y*(sqrt(((2.0*lambda))/(M_PI*(ress+1)))*sin(j*k*M_PI/(ress+1))), 0);
			W1_aux(j-1,k-1) = aux;
		}
	}

	MatrixXcd W2_aux(ress, _chiral_deg * N2);
	W2_aux.setZero();

	for (int j=1; j < ress+1; j++ ){
		for (int k=1; k < 2*N2+1; k++){
			std::complex<double> aux(y*(sqrt(((2.0*lambda))/(M_PI*(ress+1)))*sin(j*(k+(2*N1))*M_PI/(ress+1))), 0);
			W2_aux(j-1,k-1) = aux;
		}
	}

	// Uncomment the code below to verify the constraint Wp'*Wq = delta_pq //
	
	// cout << "\nW1'*W1: (It must be diagonal)\n" << W1_aux.adjoint()*W1_aux << endl;
	// cout << "\nW2'*W2: (It must be diagonal)\n" << W2_aux.adjoint()*W2_aux << endl;
	// cout << "\nW1'*W2: (It must be zero)\n" << W1_aux.adjoint()*W2_aux << endl;
	// cout << "\nW2'*W1: (It must be zero)\n" << W2_aux.adjoint()*W1_aux << endl;

	MatrixXcd W1_new(_chiral_deg * ress, _chiral_deg * N1);
	W1_new.setZero();
	W1_new << Implementing_Chiral_Symmetry_W(W1_aux, N1, ress, _chiral_deg);

	MatrixXcd W2_new(_chiral_deg * ress, _chiral_deg * N2);
	W2_new.setZero();
	W2_new << Implementing_Chiral_Symmetry_W(W2_aux, N2, ress, _chiral_deg);

	// Uncomment the code below to verify the constraint Wp'*Wq = delta_pq //
	//
	// cout << "\nW1_new'*W1_new: (It must be diagonal)\n" << W1_new.adjoint()*W1_new << endl;
	// cout << "\nW2_new'*W2_new: (It must be diagonal)\n" << W2_new.adjoint()*W2_new << endl;
	// cout << "\nW1_new'*W2_new: (It must be zero)\n" << W1_new.adjoint()*W2_new << endl;
	// cout << "\nW2_new'*W1_new: (It must be zero)\n" << W2_new.adjoint()*W1_new << endl;

	MatrixXcd W(_chiral_deg * ress, _chiral_deg * (N1 + N2));
        W << W1_new, W2_new;
	*W_pointer = W;

	// Uncomment the code below to verify the Chiral Symmetry Constraint //

	// MatrixXcd I_spin = MatrixXcd::Identity(2, 2);
	// MatrixXcd I_ress = MatrixXcd::Identity(ress, ress);
	// MatrixXcd I_n = MatrixXcd::Identity((N1+N2), (N1+N2));
	// MatrixXcd paulimatrix_y(2,2);
       
	// paulimatrix_y << 0, complex_identity,
	//	 	 -complex_identity, 0;
	 
	// cout << "\nDifference of Symmetry equation: " << (W-Kronecker_Product(paulimatrix_y, I_ress)*W*Kronecker_Product(I_n, paulimatrix_y)).cwiseAbs() << endl;
}

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

void Chiral_Orthogonal::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){

	MatrixXcd C1tio(2,2);
	MatrixXcd C2tio(2,2);

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);
	
	MatrixXcd C1_aux(2*N1, 2*N1);
	MatrixXcd C2_aux(2*N2, 2*N2);

	MatrixXcd C1(_chiral_deg * 2*N1, _chiral_deg * 2*N1);
	MatrixXcd C2(_chiral_deg * 2*N2, _chiral_deg * 2*N2);

	for (int i=1; i < 3; i++){
		for (int j=1; j < 3; j++){
			if (i == 1 && j == 1){
				std::complex<double> aux(1,0);
				C1tio(i-1,j-1) = aux;
			}
			else{
				std::complex<double> aux(0,0);
				C1tio(i-1,j-1) = aux;
			}
		}
	}

	for (int i=1; i < 3; i++){
		for (int j=1; j < 3; j++){
			if (i == 2 && j == 2){
				std::complex<double> aux(1,0);
				C2tio(i-1,j-1) = aux;
			}
			else{
				std::complex<double> aux(0,0);
				C2tio(i-1,j-1) = aux;
			}
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C1_aux.block((i-1)*identity1.rows(), (j-1)*identity1.cols(), identity1.rows(), identity1.cols()) = C1tio(i-1,j-1)*identity1;
			
		}
	}

	for (int i = 1; i < 3; i++){
		for (int j = 1; j < 3; j++){
			C2_aux.block((i-1)*identity2.rows(), (j-1)*identity2.cols(), identity2.rows(), identity2.cols()) = C2tio(i-1,j-1)*identity2;
		}
	}

	C1 << Kronecker_Product(C1_aux, MatrixXcd::Identity(2,2));
	C2 << Kronecker_Product(C2_aux, MatrixXcd::Identity(2,2));

	*C1_pointer << C1;
	*C2_pointer << C2;
}

void Chiral_Orthogonal::Create_H(MatrixXcd* H_pointer, int ress, double _lambda){

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::normal_distribution<double> distribution(0.0, 1.0);
	std::default_random_engine generator(seed);

	MatrixXcd A(ress, ress);
	MatrixXcd H1(ress, ress);
	MatrixXcd Symmetric(ress, ress);

	A.setZero();
	H1.setZero();
	Symmetric.setZero();

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		H1(i-1,i-1) = A(i-1,i-1)*((1/2)*_lambda*(1/sqrt(ress)));
		for (int j = 1 + 1; j < ress + 1; j++){
			H1(i-1,j-1) = A(i-1,j-1)*(_lambda*(1/sqrt(ress)));
		}
	}

	Symmetric << H1 + H1.transpose();

	MatrixXcd H(_chiral_deg * _spin_deg * ress, _chiral_deg * _spin_deg * ress);

	H.block(0, 0, ress, ress) =  MatrixXcd::Zero(ress, ress); H.block(0, ress, ress, ress) = Symmetric;
	H.block(ress, 0, ress, ress) = Symmetric.adjoint(); H.block(ress, ress, ress, ress) = MatrixXcd::Zero(ress, ress);

	*H_pointer = H;
}


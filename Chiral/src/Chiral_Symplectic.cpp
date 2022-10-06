#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../include/Chiral_Symplectic.h"
#include "../include/Auxiliary_Functions.h"
#include <cmath>
#include <complex>
#include <ctime>
#include <chrono>
#include <random>
#include <fstream>
#include <string>

using namespace std;

Chiral_Symplectic::Chiral_Symplectic(double lambda, int spin_deg, int chiral_deg){

	this -> _lambda = lambda;
	this -> _spin_deg = spin_deg;
	this -> _chiral_deg = chiral_deg;
}

Chiral_Symplectic::~Chiral_Symplectic() {}

void Chiral_Symplectic::Create_W(MatrixXcd* W_pointer, int ress, int N1, int N2, double lambda, double y){

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
	
	// cout << "\nW1_new'*W1_new: (It must be diagonal)\n" << W1_new.adjoint()*W1_new << endl;
	// cout << "\nW2_new'*W2_new: (It must be diagonal)\n" << W2_new.adjoint()*W2_new << endl;
	// cout << "\nW1_new'*W2_new: (It must be zero)\n" << W1_new.adjoint()*W2_new << endl;
	// cout << "\nW2_new'*W1_new: (It must be zero)\n" << W2_new.adjoint()*W1_new << endl;

	MatrixXcd W_without_spin(_chiral_deg * ress, _chiral_deg * (N1 + N2));
        W_without_spin << W1_new, W2_new;

	MatrixXcd W(_spin_deg * _chiral_deg * ress, _spin_deg * _chiral_deg * (N1 + N2));
	W << Kronecker_Product(W_without_spin, MatrixXcd::Identity(_spin_deg, _spin_deg));
	*W_pointer = W;

	// Uncomment the code below to verify the Chiral Symmetry Constraint //

	// MatrixXcd I_spin = MatrixXcd::Identity(2, 2);
	// MatrixXcd I_ress = MatrixXcd::Identity(ress, ress);
	// MatrixXcd I_n = MatrixXcd::Identity((N1+N2), (N1+N2));
	// MatrixXcd paulimatrix_z(2,2), paulimatrix_y(2,2);
       
	// paulimatrix_z << 1, 0,
	//   	 	 0, -1;
	
	// paulimatrix_y << 0, -complex_identity,
	//  		 complex_identity, 0; 

	// cout << "\nThe W constraint: \n" << (W-Kronecker_Product(Kronecker_Product(paulimatrix_z, I_ress), paulimatrix_y)*W*Kronecker_Product(Kronecker_Product(I_n, paulimatrix_z), paulimatrix_y)).cwiseAbs() << endl;
}

void Chiral_Symplectic::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);
	
	MatrixXcd C1_aux((N1+N2), (N1+N2));
	MatrixXcd C2_aux((N1+N2), (N1+N2));

	C1_aux.block(0, 0, N1, N1) << identity1; C1_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C1_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C1_aux.block(N1, N1, N2, N2) << MatrixXcd::Zero(N2, N2);

	C2_aux.block(0, 0, N1, N1) << MatrixXcd::Zero(N1, N1); C2_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C2_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C2_aux.block(N1, N1, N2, N2) << identity2;

	MatrixXcd C1_without_spin(_chiral_deg * (N1+N2), _chiral_deg * (N1+N2));
	MatrixXcd C2_without_spin(_chiral_deg * (N1+N2), _chiral_deg * (N1+N2));

	C1_without_spin << Kronecker_Product(C1_aux, MatrixXcd::Identity(_chiral_deg, _chiral_deg));
	C2_without_spin << Kronecker_Product(C2_aux, MatrixXcd::Identity(_chiral_deg, _chiral_deg));

	MatrixXcd C1(_chiral_deg * _spin_deg * (N1+N2), _chiral_deg * _spin_deg * (N1+N2));
	MatrixXcd C2(_chiral_deg * _spin_deg * (N1+N2), _chiral_deg * _spin_deg * (N1+N2));

	C1 << Kronecker_Product(C1_without_spin, MatrixXcd::Identity(_spin_deg, _spin_deg));
	C2 << Kronecker_Product(C2_without_spin, MatrixXcd::Identity(_spin_deg, _spin_deg));

	*C1_pointer << C1;
	*C2_pointer << C2;
}

void Chiral_Symplectic::Create_H(MatrixXcd* H_pointer, int ress, double _lambda){

	complex<double> complex_identity(0,1);

	MatrixXcd paulimatrix_x(2,2);
	MatrixXcd paulimatrix_y(2,2);
	MatrixXcd paulimatrix_z(2,2);

	paulimatrix_x << 0, 1,
		  	 1, 0;

	paulimatrix_y << 0, complex_identity,
		 	 -complex_identity, 0;

	paulimatrix_z << 1, 0,
			 0, -1;

	auto seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::normal_distribution<double> distribution(0.0, 1.0);
	std::default_random_engine generator(seed);

	MatrixXcd A(ress, ress);
	MatrixXcd B(ress, ress);
	MatrixXcd C(ress, ress);
	MatrixXcd D(ress, ress);

	MatrixXcd H0(ress, ress);
	MatrixXcd H1(ress, ress);
	MatrixXcd H2(ress, ress);
	MatrixXcd H3(ress, ress);

	MatrixXcd Ham(_spin_deg * ress, _spin_deg * ress);

	A.setZero();
	B.setZero();
	C.setZero();
	D.setZero();

	H0.setZero();
	H1.setZero();
	H2.setZero();
	H3.setZero();

	Ham.setZero();

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			B(i-1,j-1) = aux;
		}
	}
	
	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			C(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			D(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		H0(i-1,i-1) = A(i-1,i-1)*(_lambda*(1/sqrt(4*ress)));
		for (int j = i + 1; j < ress + 1; j++){
			H0(i-1,j-1) = A(i-1,j-1)*(_lambda*(1/sqrt(4*ress)));
			H0(j-1,i-1) = A(j-1,i-1)*(_lambda*(1/sqrt(4*ress)));
		}
	}

	H1 << (_lambda*(1/sqrt(4*ress)))*B;
	H2 << (_lambda*(1/sqrt(4*ress)))*C;
	H3 << (_lambda*(1/sqrt(4*ress)))*D;

	Ham << Kronecker_Product(H0, MatrixXcd::Identity(2,2)) + complex_identity*(Kronecker_Product(H1, paulimatrix_x) + Kronecker_Product(H2, paulimatrix_y) + Kronecker_Product(H3, paulimatrix_z));

	MatrixXcd H(_chiral_deg * _spin_deg * ress, _chiral_deg * _spin_deg * ress);

	H.block(0, 0, _spin_deg * ress, _spin_deg * ress) =  MatrixXcd::Zero(_spin_deg * ress, _spin_deg * ress); H.block(0, _spin_deg * ress, _spin_deg * ress, _spin_deg * ress) = Ham;
	H.block(_spin_deg * ress, 0, _spin_deg * ress, _spin_deg * ress) = Ham.adjoint(); H.block(_spin_deg * ress, _spin_deg * ress, _spin_deg * ress, _spin_deg * ress) = MatrixXcd::Zero(_spin_deg * ress, _spin_deg * ress);

	// Uncomment the code below to verify the Chiral Symmetry Constraint of Hamiltonian //

	// cout << "\nThe Hamiltonian Constraint: \n" << (H -(-Kronecker_Product(Kronecker_Product(paulimatrix_z, MatrixXcd::Identity(ress, ress)), MatrixXcd::Identity(_spin_deg, _spin_deg))*H*Kronecker_Product(Kronecker_Product(paulimatrix_z, MatrixXcd::Identity(ress, ress)), MatrixXcd::Identity(_spin_deg, _spin_deg)))) << endl;

	*H_pointer << H;
}

void Chiral_Symplectic::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps){
	std::ofstream output_G("Data_Analysis/Channel/Dirac_G_S_Channel.txt");
	std::ofstream output_P("Data_Analysis/Channel/Dirac_P_S_Channel.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 10; j++){
			if (j == 9){
				output_G << G(i,j).real() << std::endl;
				output_P << P(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
				output_P << P(i,j).real() << "\t";
			}
		}
	}	
}

void Chiral_Symplectic::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1){
	std::ofstream output_G("Data_Analysis/Gamma/Dirac_G_S_Gamma_N"+to_string(N1)+".txt");
	std::ofstream output_P("Data_Analysis/Gamma/Dirac_P_S_Gamma_N"+to_string(N1)+".txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 21; j++){
			if (j == 20){
				output_G << G(i,j).real() << std::endl;
				output_P << P(i,j).real() << std::endl;
			}
			else{
				output_G << G(i,j).real() << "\t";
				output_P << P(i,j).real() << "\t";
			}
		}
	}	
}

void Chiral_Symplectic::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps) {}

void Chiral_Symplectic::Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps) {}

void Chiral_Symplectic::Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps) {}

void Chiral_Symplectic::Save_txt_files_Bell_Parameter_Fixed_Base(MatrixXd Bell_Parameter_Fixed_Base, int num_steps) {}

void Chiral_Symplectic::Save_txt_files_Correlators_Bell_Inequality_Gamma(MatrixXd Correlator_C11, MatrixXd Correlator_C22, MatrixXd Correlator_C12, MatrixXd Correlator_C21, int num_steps) {}

void Chiral_Symplectic::Save_txt_files_Energy(MatrixXcd G, int num_steps, int N1){}

void Chiral_Symplectic::Save_txt_files_Energy_Gamma(MatrixXcd G, int num_steps, int N1, int gamma_idx) {}

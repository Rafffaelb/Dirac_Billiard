#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../include/Chiral_Orthogonal.h"
#include "../include/Auxiliary_Functions.h"
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
	// MatrixXcd paulimatrix_z(2,2);
       
	// paulimatrix_z << 1, 0,
	//  	 	 0, -1;
	 
	// cout << "\nThe W constraint: " << (W-Kronecker_Product(paulimatrix_z, I_ress)*W*Kronecker_Product(I_n, paulimatrix_z)).cwiseAbs() << endl;
}

void Chiral_Orthogonal::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2){

	MatrixXcd identity1 = MatrixXcd::Identity(N1,N1);
	MatrixXcd identity2 = MatrixXcd::Identity(N2,N2);
	
	MatrixXcd C1_aux((N1+N2), (N1+N2));
	MatrixXcd C2_aux((N1+N2), (N1+N2));

	C1_aux.block(0, 0, N1, N1) << identity1; C1_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C1_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C1_aux.block(N1, N1, N2, N2) << MatrixXcd::Zero(N2, N2);

	C2_aux.block(0, 0, N1, N1) << MatrixXcd::Zero(N1, N1); C2_aux.block(0, N1, N1, N2) << MatrixXcd::Zero(N1, N2);
	C2_aux.block(N1, 0, N2, N1) << MatrixXcd::Zero(N2, N1); C2_aux.block(N1, N1, N2, N2) << identity2;

	MatrixXcd C1(_chiral_deg * (N1+N2), _chiral_deg * (N1+N2));
	MatrixXcd C2(_chiral_deg * (N1+N2), _chiral_deg * (N1+N2));

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

	for (int i = 1; i < ress + 1; i++){
		for (int j = 1; j < ress + 1; j++){
			double aux = distribution(generator);
			A(i-1,j-1) = aux;
		}
	}

	for (int i = 1; i < ress + 1; i++){
		H1(i-1,i-1) = A(i-1,i-1)*(_lambda*(1/sqrt(ress)));
		for (int j = i + 1; j < ress + 1; j++){
			H1(i-1,j-1) = A(i-1,j-1)*(_lambda*(1/sqrt(2*ress)));
			H1(j-1,i-1) = A(j-1,i-1)*(_lambda*(1/sqrt(2*ress)));
		}
	}

	MatrixXcd Ham(ress, ress);
       	Ham << H1;

	MatrixXcd H(_chiral_deg * _spin_deg * ress, _chiral_deg * _spin_deg * ress);
	H.setZero();

	H.block(0, 0, ress, ress) =  MatrixXcd::Zero(ress, ress); H.block(0, ress, ress, ress) = Ham;
	H.block(ress, 0, ress, ress) = Ham.adjoint(); H.block(ress, ress, ress, ress) = MatrixXcd::Zero(ress, ress);

	// Uncomment the code below to verify the Chiral Symmetry Constraint of Hamiltonian //

	// MatrixXcd paulimatrix_z(2,2);
       
	// paulimatrix_z << 1, 0,
	//		 0, -1;

	// cout << "\nThe Hamiltonian Constraint: \n" << (H -(-Kronecker_Product(paulimatrix_z, MatrixXcd::Identity(ress, ress))*H*Kronecker_Product(paulimatrix_z, MatrixXcd::Identity(ress, ress)))) << endl;

	*H_pointer = H;
}

void Chiral_Orthogonal::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps){
	std::ofstream output_G("Data_Analysis/Channel/Dirac_G_O_Channel.txt");
	std::ofstream output_P("Data_Analysis/Channel/Dirac_P_O_Channel.txt");

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

void Chiral_Orthogonal::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1){

	std::ofstream output_G("Data_Analysis/Gamma/Dirac_G_O_Gamma_N"+to_string(N1)+".txt");
	std::ofstream output_P("Data_Analysis/Gamma/Dirac_P_O_Gamma_N"+to_string(N1)+".txt");

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

void Chiral_Orthogonal::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps){

	std::ofstream output_Concurrence("Data_Analysis/Concurrence/Dirac_Concurrence_O_Gamma.txt");
	std::ofstream output_Entanglement("Data_Analysis/Concurrence/Dirac_Entanglement_O_Gamma.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 21; j++){
			if (j == 20){
				output_Concurrence << Concurrence(i,j) << std::endl;
				output_Entanglement << Entanglement(i,j) << std::endl;
			}
			else{
				output_Concurrence << Concurrence(i,j) << "\t";
				output_Entanglement << Entanglement(i,j) << "\t";
			}
		}
	}	
}

void Chiral_Orthogonal::Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps){

	std::ofstream output_Bell_Parameter_Ress("Data_Analysis/Bell_Parameter/Bell_Ress/Dirac_Bell_Parameter_O_Ress.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 11; j++){
			if (j == 10){
				output_Bell_Parameter_Ress << Bell_Parameter_Ress(i,j) << std::endl;
			}
			else{
				output_Bell_Parameter_Ress << Bell_Parameter_Ress(i,j) << "\t";
			}
		}
	}
}

void Chiral_Orthogonal::Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps){

	std::ofstream output_Bell_Parameter_Gamma("Data_Analysis/Bell_Parameter/Bell_Gamma/Dirac_Bell_Parameter_O_Gamma.txt");
	std::ofstream output_Bell_Parameter_Dephase_Gamma("Data_Analysis/Bell_Parameter/Bell_Gamma/Dirac_Bell_Parameter_Dephase_O_Gamma.txt");

	for(int i = 0; i < num_steps; i++){
		for (int j = 0; j < 21; j++){
			if (j == 20){
				output_Bell_Parameter_Gamma << Bell_Parameter_Gamma(i,j) << std::endl;
				output_Bell_Parameter_Dephase_Gamma << Bell_Parameter_Dephase_Gamma(i,j) << std::endl;
			}
			else{
				output_Bell_Parameter_Gamma << Bell_Parameter_Gamma(i,j) << "\t";
				output_Bell_Parameter_Dephase_Gamma << Bell_Parameter_Dephase_Gamma(i,j) << "\t";
			}
		}
	}
}

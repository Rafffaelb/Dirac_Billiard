#include <iostream>
#include <chrono>
#include "../../include/Chiral_AbstractClass_h/Chiral.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include "../../include/Chiral_AbstractClass_h/Run_Simulation_Bell_Parameter_Fixed_Base.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

void Chiral::Run_Simulation_Bell_Parameter_Fixed_Base(){

	auto start = chrono::system_clock::now();

	double Gamma, y;
       	int ress, N1, N2, n, num_steps_aux;

	num_steps_aux = 4000000;
	ress = 100;

	Gamma = 1;

	y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));

	N1 = 1;
	N2 = N1;
	n = N1 + N2;

	// Create_ProjectionMatrices //

	MatrixXcd C1(_chiral_deg * _spin_deg * 2*N1, _chiral_deg * _spin_deg * 2*N1); MatrixXcd C2(_chiral_deg * _spin_deg * 2*N2, _chiral_deg * _spin_deg * 2*N2);
	C1.setZero(); C2.setZero();
	MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

	Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);
	
	MatrixXd Bell_Parameter_Fixed_Base(num_steps_aux, 2);

	Bell_Parameter_Fixed_Base.setZero();

	// Create W Matrices //

	MatrixXcd W(_chiral_deg * _spin_deg * ress, _chiral_deg * _spin_deg * n);
	W.setZero();
	MatrixXcd *W_pointer = &W;
	
	Create_W(W_pointer, ress, N1, N2, _lambda, y);
		
	#pragma omp parallel for shared(W, C1, C2)
	for (int step = 1; step < num_steps_aux + 1; step++){
		
		// Generate Hamiltonian Matrix //

		MatrixXcd H(_chiral_deg * _spin_deg * ress, _chiral_deg * _spin_deg * ress);
		H.setZero();
		MatrixXcd* H_pointer = &H;
		
		Create_H(H_pointer, ress, _lambda);

		// Create billiard setup //

		Quantum_chaotic_billiard billiard_setup(H, W, C1, C2, N1, N2);

		// Scattering Matrix //
		
		billiard_setup.Calculate_Smatrix();

		// Bell Parameter //
	
		billiard_setup.Calculate_Concurrence();	
		billiard_setup.Calculate_Bell_Parameter_Fixed_Base();

		Bell_Parameter_Fixed_Base(step-1, 0) = billiard_setup.getBell_Parameter();
		Bell_Parameter_Fixed_Base(step-1, 1) = billiard_setup.getConcurrence();

		if (step % _num_steps == 0){
			std::cout << "\nCurrent number of steps: " << step << std::endl;
		}

	}
	
	//Save Concurrence matrix as txt files //
	Save_txt_files_Bell_Parameter_Fixed_Base(Bell_Parameter_Fixed_Base, num_steps_aux);
		
	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}

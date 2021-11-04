#include <iostream>
#include <chrono>
#include "../../include/Chiral_AbstractClass_h/Chiral.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include "../../include/Chiral_AbstractClass_h/Run_Simulation_Concurrence_Gamma.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

void Chiral::Run_Simulation_Concurrence_Gamma(){

	auto start = chrono::system_clock::now();

	double Gamma, y;
       	int ress, N1, N2, n;

	ress = 100;

	N1 = 1;
	N2 = N1;
	n = N1 + N2;

	// Create_ProjectionMatrices //

	MatrixXcd C1(_chiral_deg * _spin_deg * 2*N1, _chiral_deg * _spin_deg * 2*N1); MatrixXcd C2(_chiral_deg *_spin_deg * 2*N2, _chiral_deg * _spin_deg * 2*N2);
	C1.setZero(); C2.setZero();
	MatrixXcd *C1_pointer = &C1; MatrixXcd *C2_pointer = &C2;

	Create_ProjectionMatrices(C1_pointer, C2_pointer, N1, N2);	
	
	MatrixXd Concurrence(_num_steps, 21);
	MatrixXd Entanglement(_num_steps,21);

	Concurrence.setZero();
	Entanglement.setZero();

	for (int gamma_idx = 1; gamma_idx < 22; gamma_idx++){
	
		if (gamma_idx == 1){
				
			Gamma = 0.0001;	
		}
		else{

			Gamma = double(gamma_idx - 1) / double(20);
		}

		double y = sqrt(double(1.0)/Gamma)*(1.0-sqrt(1.0-Gamma));

		// Create W Matrices //

		MatrixXcd W(_chiral_deg * _spin_deg * ress, _chiral_deg * _spin_deg * n);
		W.setZero();
		MatrixXcd *W_pointer = &W;

		Create_W(W_pointer, ress, N1, N2, _lambda, y);		
		
		#pragma omp parallel for shared(W, C1, C2)
		for (int step = 1; step < _num_steps + 1; step++){
		
			// Generate Hamiltonian Matrix //

			MatrixXcd H(_chiral_deg * _spin_deg * ress, _chiral_deg * _spin_deg * ress);
			H.setZero();
			MatrixXcd* H_pointer = &H;

			Create_H(H_pointer, ress, _lambda);
			
			// Create billiard setup //

			Quantum_chaotic_billiard billiard_setup(H, W, C1, C2, N1, N2);
			
			// Scattering Matrix //
		
			billiard_setup.Calculate_Smatrix();

			// Concurrence //
		
			billiard_setup.Calculate_Concurrence();

			Concurrence(step-1, gamma_idx-1) = billiard_setup.getConcurrence();
			Entanglement(step-1, gamma_idx-1) = billiard_setup.getEntanglement();


			if (step % _num_steps == 0){
				std::cout << "\nCurrent number of steps: " << step << "| Current index of Gamma: " << gamma_idx << std::endl;
			}

		}
		//Save Concurrence matrix as txt files //
		Save_txt_files_Concurrence_Gamma(Concurrence, Entanglement, _num_steps);
		
	}

	auto end = chrono::system_clock::now();
	auto elapsed =
		chrono::duration_cast<chrono::minutes>(end-start);
	cout << "\n Simulation time duration: " << elapsed.count() << "\n";
}
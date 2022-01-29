#include <iostream>
#include <cstring>
#include "include/Chiral_AbstractClass_h/Chiral.h"
#include <cmath>
#include "include/Chiral_Orthogonal.h"
#include "include/Chiral_Unitary.h"
#include "include/Chiral_Symplectic.h"

using namespace std;

int main(int argc, char **argv){

	double lambda;
	int num_steps, spin_deg, chiral_deg;

	lambda = 0.5;
	chiral_deg = 2;
	num_steps = 1000000;

	for (int i = 1; i < argc; i++){
		
		if (strcmp(argv[i],"Chiral_Orthogonal") == 0){
		
			spin_deg = 1;
			
			Chiral_Orthogonal chiral_orthogonal(lambda, num_steps, spin_deg, chiral_deg);

			for (int j = 1; j < argc; j++){
				
				if (strcmp(argv[j],"Channel") == 0){
					
					cout << "\n ####### Running Chiral Orthogonal (variable: Channel) ####### \n" << endl;
				
					chiral_orthogonal.Run_Simulation_Conductance_Channels();
				}

				if (strcmp(argv[j],"Gamma") == 0){

					cout << "\n ##### Running Chiral Orthogonal (variable: Gamma) ##### \n" << endl;

					chiral_orthogonal.Run_Simulation_Conductance_Gamma();
				}

				if (strcmp(argv[j],"Concurrence") == 0){

					cout << "\n ###### Running Orthogonal Concurrence (variable: Gamma) ###### \n" << endl;
				
					chiral_orthogonal.Run_Simulation_Concurrence_Gamma();
				}

				if (strcmp(argv[j],"Bell_Parameter_Ress") == 0){
					
					cout << "\n ###### Running Orthogonal Bell Parameter (variable: Ress) ##### \n" << endl;
				
					chiral_orthogonal.Run_Simulation_Bell_Parameter_Ress();
				}

				if (strcmp(argv[j],"Bell_Parameter_Gamma") == 0){
		
					cout << "\n ###### Running Orthogonal Bell Parameter (variable: Gamma) #### \n" << endl;
					
					chiral_orthogonal.Run_Simulation_Bell_Parameter_Gamma();
				}

				if (strcmp(argv[j],"Bell_Parameter_Fixed_Base") == 0){
				
					cout << "\n ###### Running Orthogonal Bell Parameter Fixed Base ##### \n" << endl;

					chiral_orthogonal.Run_Simulation_Bell_Parameter_Fixed_Base();
				}

				if (strcmp(argv[j],"Correlators_Gamma") == 0){
		
					cout << "\n ###### Running Orthogonal Correlators C_ij (variable: Gamma) #### \n" << endl;
					
					chiral_orthogonal.Run_Simulation_Correlators_Bell_Inequality_Gamma();
				}

			}
			chiral_orthogonal.~Chiral_Orthogonal();
		}
		else{
			if (strcmp(argv[i],"Chiral_Unitary") == 0){
				
				spin_deg = 1;

				Chiral_Unitary chiral_unitary(lambda, num_steps, spin_deg, chiral_deg);

				for (int j = 1; j < argc; j++){

					if (strcmp(argv[j],"Channel") == 0){
					
						cout << "\n ####### Running Chiral Unitary (variable: Channel) ####### \n" << endl;

						chiral_unitary.Run_Simulation_Conductance_Channels();
					}

					if (strcmp(argv[j],"Gamma") == 0){

						cout << "\n ###### Running Chiral Unitary (variable: Gamma) ##### \n" << endl;

						chiral_unitary.Run_Simulation_Conductance_Gamma();
					}
					
					if (strcmp(argv[j],"Concurrence") == 0){

						cout << "\n ###### Running Unitary Concurrence (variable: Gamma) ###### \n" << endl;

						chiral_unitary.Run_Simulation_Concurrence_Gamma();
					}

					if (strcmp(argv[j],"Bell_Parameter_Ress") == 0){
					
						cout << "\n ###### Running Unitary Bell Parameter (variable: Ress) ##### \n" << endl;
					
						chiral_unitary.Run_Simulation_Bell_Parameter_Ress();
					}

					if (strcmp(argv[j],"Bell_Parameter_Gamma") == 0){
						
						cout << "\n ##### Running Unitary Bell Parameter (variable: Gamma) ##### \n" << endl;

						chiral_unitary.Run_Simulation_Bell_Parameter_Gamma();
					}

					if (strcmp(argv[j],"Bell_Parameter_Fixed_Base") == 0){

						cout << "\n ##### Running Unitary Bell Parameter Fixed Base ##### \n" << endl;

						chiral_unitary.Run_Simulation_Bell_Parameter_Fixed_Base();
					}

					if (strcmp(argv[j],"Correlators_Gamma") == 0){
		
						cout << "\n ###### Running Unitary Correlators C_ij (variable: Gamma) #### \n" << endl;
					
						chiral_unitary.Run_Simulation_Correlators_Bell_Inequality_Gamma();
					}

				}
				chiral_unitary.~Chiral_Unitary();
			}
			else{
				if (strcmp(argv[i],"Chiral_Symplectic") == 0){
				
					spin_deg = 2;

					Chiral_Symplectic chiral_symplectic(lambda, num_steps, spin_deg, chiral_deg);

					for (int j = 1; j < argc; j++){

						if (strcmp(argv[j],"Channel") == 0){
					
							cout << "\n ####### Running Chiral Symplectic (variable: Channel) ####### \n" << endl;

							chiral_symplectic.Run_Simulation_Conductance_Channels();
						}

						if (strcmp(argv[j],"Gamma") == 0){
				
							cout << "\n ###### Running Chiral Symplectic (variable: Gamma) ####### \n" << endl;
							chiral_symplectic.Run_Simulation_Conductance_Gamma();
						}

					}
					chiral_symplectic.~Chiral_Symplectic();
				}
			}
		}
	}
	return 0;
}

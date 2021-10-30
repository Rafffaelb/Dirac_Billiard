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
	num_steps = 100000;

	for (int i = 1; i < argc; i++){
		
		if (strcmp(argv[i],"Chiral_Orthogonal") == 0){
		
			spin_deg = 1;
			
			Chiral_Orthogonal chiral_orthogonal(lambda, num_steps, spin_deg, chiral_deg);

			for (int j = 1; j < argc; j++){
				
				if (strcmp(argv[j],"Channel") == 0){
					
					cout << "\n ####### Running Chiral Orthogonal (variable: Channel) ####### \n" << endl;
				
					chiral_orthogonal.Run_Simulation_Conductance_Channels();
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
					}
					chiral_symplectic.~Chiral_Symplectic();
				}
			}
		}
	}
	return 0;
}

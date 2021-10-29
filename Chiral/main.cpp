#include <iostream>
#include <cstring>
#include "include/Chiral_AbstractClass_h/Chiral.h"
#include <cmath>
#include "include/Chiral_Orthogonal.h"

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
	}
	return 0;
}

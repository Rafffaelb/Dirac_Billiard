#ifndef CHIRAL_SYMPLECTIC_H
#define CHIRAL_SYMPLECTIC_H

#include "Chiral_AbstractClass_h/Chiral.h"

class Chiral_Symplectic: public Chiral{

	public:
		Chiral_Symplectic(double lambda, int num_steps, int spin_deg, int chiral_deg);

		~Chiral_Symplectic();

		void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double _y);

		void Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2);

		void Create_H(MatrixXcd* H_pointer, int _ress, double _lambda);
		
		void Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps);
		void Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1);
};
#endif

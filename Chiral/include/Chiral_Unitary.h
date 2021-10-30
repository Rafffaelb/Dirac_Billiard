#ifndef CHIRAL_UNITARY_H
#define CHIRAL_UNITARY_H

#include "Chiral_AbstractClass_h/Chiral.h"

class Chiral_Unitary: public Chiral{

	public:
		Chiral_Unitary(double lambda, int num_steps, int spin_deg, int chiral_deg);

		~Chiral_Unitary();

		void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double _y);

		void Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2);

		void Create_H(MatrixXcd* H_pointer, int _ress, double _lambda);
		
		void Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps);
};
#endif
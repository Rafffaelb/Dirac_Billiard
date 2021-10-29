#ifndef CHIRAL_ORTHOGONAL_H
#define CHIRAL_ORTHOGONAL_H

#include "Chiral_AbstractClass_h/Chiral.h"

class Chiral_Orthogonal: public Chiral{

	public:
		Chiral_Orthogonal(double lambda, int num_steps, int spin_deg, int chiral_deg);

		~Chiral_Orthogonal();

		void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double _y);

		void Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2);

		void Create_H(MatrixXcd* H_pointer, int _ress, double _lambda);
};
#endif

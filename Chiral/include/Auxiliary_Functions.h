#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

using namespace std;
using namespace Eigen;

MatrixXcd Kronecker_Product(MatrixXcd A, MatrixXcd B);

MatrixXcd Implementing_Chiral_Symmetry_W(MatrixXcd W_aux, int N, int ress, int _chiral_deg);
	
#endif

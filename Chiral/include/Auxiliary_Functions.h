#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

using namespace std;
using namespace Eigen;

MatrixXcd Kronecker_Product(MatrixXcd A, MatrixXcd B);

MatrixXcd Implementing_Chiral_Symmetry_W(MatrixXcd W_aux, int N, int ress, int _chiral_deg);

MatrixXcd product_vector_matrix(Vector3d vector, MatrixXcd paulimatrix_x, MatrixXcd paulimatrix_y, MatrixXcd paulimatrix_z);

#endif

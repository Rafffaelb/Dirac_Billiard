#include <iostream>
#include <chrono>
#include "../../include/Chiral_AbstractClass_h/Chiral.h"
#include "../../include/Quantum_chaotic_billiard.h"
#include <eigen3/Eigen/Dense>
#include <omp.h>

using namespace std;
using namespace Eigen;

Chiral::Chiral() {};

Chiral::~Chiral() {};

void Chiral::Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double y) {};

void Chiral::Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2) {};

void Chiral::Create_H(MatrixXcd* H_pointer, int _ress, double _lambda) {};

void Chiral::Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps) {};

void Chiral::Save_txt_files_Gamma(MatrixXcd G, MatrixXcd P, int num_steps, int N1) {};

void Chiral::Save_txt_files_Concurrence_Gamma(MatrixXd Concurrence, MatrixXd Entanglement, int num_steps) {};

void Chiral::Save_txt_files_Bell_Parameter_Ress(MatrixXd Bell_Parameter_Ress, int num_steps) {};

void Chiral::Save_txt_files_Bell_Parameter_Gamma(MatrixXd Bell_Parameter_Gamma, MatrixXd Bell_Parameter_Dephase_Gamma, int num_steps) {};

void Chiral::Save_txt_files_Bell_Parameter_Fixed_Base(MatrixXd Bell_Parameter_Fixed_Base, int num_steps_aux) {};

void Chiral::Save_txt_files_Correlators_Bell_Inequality_Gamma(MatrixXd Correlator_C11, MatrixXd Correlator_C22, MatrixXd Correlator_C12, MatrixXd Correlator_C21, int num_steps) {};

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

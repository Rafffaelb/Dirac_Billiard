#include <iostream>
#ifndef QUANTUM_CHAOTIC_BILLIARD_H
#define QUANTUM_CHAOTIC_BILLIARD_H

using namespace std;
using namespace Eigen;

class Quantum_chaotic_billiard
{
	private:
		int _N1;
		int _N2;
		
		complex<double> _G;
		complex<double> _P;
		double _Concurrence;
		double _Entanglement;
		double _Bell_Parameter;
		double _Bell_Parameter_Dephase;
		double _Correlator_C11;
		double _Correlator_C22;
		double _Correlator_C12;
		double _Correlator_C21;

		MatrixXcd _H;
		MatrixXcd _W;
		MatrixXcd _C1;
		MatrixXcd _C2;
		MatrixXcd _S;
	
	public:
		Quantum_chaotic_billiard(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2, int N1, int N2);

		void Set_Setup(MatrixXcd H, MatrixXcd W, MatrixXcd C1, MatrixXcd C2, int N1, int N2);

		void Calculate_Smatrix(double Energy);

		void Calculate_G_and_P();

		complex<double> getG();

		complex<double> getP();

		void Calculate_Concurrence();

		double getConcurrence();
		
		double getEntanglement();

		void Calculate_Bell_Parameter();

		void Calculate_Bell_Parameter_Fixed_Base();

		double getBell_Parameter();

		double getBell_Parameter_Dephase();

		void Calculate_Correlators();

		double getCorrelator_C11();
		double getCorrelator_C22();
		double getCorrelator_C12();
		double getCorrelator_C21();
};

#endif

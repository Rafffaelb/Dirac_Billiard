#include <eigen3/Eigen/Dense>
#ifndef CHIRAL_H
#define CHIRAL_H

using namespace Eigen;

class Chiral{

	protected:

		int _num_steps;
		int _spin_deg;
		int _chiral_deg;
		double _lambda;

	public:

		Chiral();

		virtual ~Chiral();

		virtual void Create_W(MatrixXcd* W_pointer, int _ress, int N1, int N2, double _lambda, double y) = 0;
		
		void Run_Simulation_Conductance_Channels();

		virtual void Create_H(MatrixXcd* H_pointer, int _ress, double _lambda) = 0;

		virtual void Create_ProjectionMatrices(MatrixXcd* C1_pointer, MatrixXcd* C2_pointer, int N1, int N2) = 0;

		virtual void Save_txt_files_Channels(MatrixXcd G, MatrixXcd P, int num_steps) = 0;
};

#endif

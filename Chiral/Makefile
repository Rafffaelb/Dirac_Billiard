CC = g++

all:
	$(CC) -mfma -fopenmp -O3 main.cpp src/Chiral_AbstractClass_cpp/Chiral.cpp src/Chiral_Orthogonal.cpp src/Chiral_Unitary.cpp src/Chiral_Symplectic.cpp src/Quantum_chaotic_billiard.cpp src/Auxiliary_Functions.cpp src/Chiral_AbstractClass_cpp/Run_Simulation_Conductance_Channels.cpp src/Chiral_AbstractClass_cpp/Run_Simulation_Conductance_Gamma.cpp src/Chiral_AbstractClass_cpp/Run_Simulation_Concurrence_Gamma.cpp src/Chiral_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Ress.cpp src/Chiral_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Gamma.cpp src/Chiral_AbstractClass_cpp/Run_Simulation_Bell_Parameter_Fixed_Base.cpp -o Chiral

clean:
	rm Chiral

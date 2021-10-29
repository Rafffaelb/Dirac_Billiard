CC = g++

all:
	$(CC) -mfma -fopenmp -O3 main.cpp src/Chiral_AbstractClass_cpp/Chiral.cpp src/Chiral_Orthogonal.cpp src/Chiral_Unitary.cpp src/Chiral_Symplectic.cpp src/Quantum_chaotic_billiard.cpp src/Auxiliary_Functions.cpp src/Chiral_AbstractClass_cpp/Run_Simulation_Conductance_Channels.cpp -o Chiral

clean:
	rm Chiral

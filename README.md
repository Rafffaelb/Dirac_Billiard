# Dirac's Billiard

Hello! My name is Rafael and in this repository I perform numerical simulations to calculate properties of electrical transport in Dirac's billiard such as conductance, power shot noise, concurrence, the entanglement of formation, and violations of Bell's inequality. All of the results of this repository can be found in the jupyter notebook's files of the directory "/Chiral/Data_Analysis". Dirac's billiard is a type of quantum billiard used to describe a system with chiral/sub-lattice symmetry. Systems with chiral symmetry can be found in some categories of topological insulators and in single-layer graphene quantum billiard, illustrated in the Figure 1 from the thesis [[1]](#1) that produces the reference [[2]](#2).

<p align="center" width="100%">
  Figure 1<br>
    <img width="40%" src="https://user-images.githubusercontent.com/53999015/147810879-3242820a-749b-4b06-9f81-92ba5c8cebe9.png">
</p>

The difference between Schrödinger and Dirac billiards remains mainly in the electron behavior described by the Schrödinger and Dirac equations, respectively. Moreover, the dimension of hamiltonian and coupling matrix between lead and cavity are double due to the chiral/sub-lattice degree of freedom. This difference can be viewed by the peculiar form of Dirac billiard hamiltonian:

<p align="center" width="100%">
    <img width="30%" src="https://user-images.githubusercontent.com/53999015/149149596-7ca5f772-7e74-4896-931c-605de358dac9.png">
</p>

Here, the T-matrix block of the H-matrix has dimension M x M, where M is the number of resonances.

## Conductance

**In this repository**, I analyze the influence of the three fundamental symmetries (chGOE, chGUE, chGSE) in the conductance of electrons inside the Dirac billiard. I obtain the distributions of conductance, as well as the average and variance in terms of the number of open channels for each class of symmetry. I also perform simulations to see how the conductance of each class of symmetry depends on the barrier between the lead and chaotic cavity, which Γ is called the transparency of the barrier. Moreover, I analyze the weak localization correction of conductance as a function of the barrier's transparency Γ and compare with the Schrodinger billiard.

As the Schrödinger billiard, the conductance in a Dirac billliard also depends on the number of open channels in each lead, the potential barrier between the lead and the chaotic resonance cavity and especially the class of symmetry of the hamiltonian. Each open channels has two subchannels due to the chirality of the system. In graphene material, this subchannels can be viewed as the two sublattices in the honeycomb lattice.

<p align="center" width="100%">
  Figure 2<br>
    <img width="60%" src="https://user-images.githubusercontent.com/53999015/148553888-1240bffa-9381-4286-a088-9cc90bc7f55c.png">
</p>

As we can see in Figure 2, the distribution of realizations of conductance depends on the β which indicates a specific symmetry of the hamiltonian. These symmetries are fundamentals to describe the physics of quantum billiards. In a Dirac billiard, there are three fundamental symmetries which are listed below:

#### Chiral Gaussian Orthogonal Ensemble (chGOE)
The index β = 1 (blue), also indicated in the literature by "BDI", represents a system that has time-reversal symmetry (TRS) and spin rotational symmetry (SRS). The Hamiltonian of this system has to be a real matrix.
#### Chiral Gaussian Unitary Ensemble (chGUE)
The β = 2 (red), also indicated by literature as "AIII", represents a system that the TRS is broken. This happens when there is a magnetic field applied in the quantum dot, the magnetic field breaks the time-reversal symmetry. When the TRS is broken, it doesn't matter if SRS is preserved or broken. In this case, each element of the Hamiltonian is a complex number.
#### Chiral Gaussian Symplectic Ensemble (chGSE)
The β = 4 (green), indicated by literature as "CII", represents a system that SRS is broken but TRS is preserved. In this case, the Hamiltonian has a quaternionic structure. In these systems, there is a special interaction called spin-orbit coupling that arises due to the relation between the apparent magnetic field seen from the electron perspective and the magnetic moment of the electron associated with its intrinsic spin. 

## Shot Noise

**In this repository**, I perform numerical simulations of shot-noise power, also known as the spectral density of the noise. In the same way as the conductance, I analyze the influence of the chGOE, chGUE, and chGSE ensembles on the distributions of shot-noise power, as well as the average and variance. I also explore the dependence of shot-noise power with the number of open channels and the barrier's transparency Γ, recovering some results of the literature. The shot-noise produced by scattered electrons in a chaotic cavity can be used to indicate the presence of entangled electron pairs by violations of Bell's inequality.

## Concurrence

**In this repository**, I analyze the influence of the time-reversal symmetry in the concurrence and entanglement of formation for the Dirac billiard. So I study the concurrence and entanglement of formation produced in the chiral/sub-lattice degree of freedom and compare this production of entangled states when the time-reversal symmetry is broken or preserved. I also obtain the distributions of concurrence and entanglement of formation, as well as the average and variance in terms of the barrier's transparency Γ, recovering numerical and analytical results from the literature [[3]](#3). Comparing Schrödinger and Dirac repositories, the Dirac billiard produces more entanglement pairs than Schrödinger one.

The concurrence is a property that quantifies the degree of entanglement of a system [[4]](#4) and can be viewed as how much a state of the system can be separated in two states from different subsystems. In this repository we set the number of open channels in each lead as N = 1, the main goal is to study the production of entanglement in chiral/sub-lattice degrees of freedom. 

## References

<a id="1">[1]</a>
Barros, M. S. M. (2018). Teoria do ruído e fenômenos de interferência quântica em nano estruturas quirais. Institutional Repository of UFPB. https://repositorio.ufpb.br/jspui/handle/123456789/14825

<a id="2">[2]</a>
Barros, M. S. M., Júnior, A. J. N., Macedo-Junior, A. F., Ramos, J. G. G. S., & Barbosa, A. L. R. (2013). Open chaotic Dirac billiards: Weak (anti)localization, conductance fluctuations, and decoherence. Physical Review B, 88(24). doi:10.1103/physrevb.88.245133.

<a id="3">[3]</a>
Anomalous entanglement in chaotic Dirac billiards
J. G. G. S. Ramos, I. M. L. da Silva, and A. L. R. Barbosa
Phys. Rev. B 90, 245107 – Published 1 December 2014

<a id="4">[4]</a>
W. K. Wootters, Entanglement of formation of an arbitrary state of two qubits, Physical Review Letters 80, 2245 (1998).

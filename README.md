# Dirac's Billiard

Hello! My name is Rafael and in this repository I perform numerical simulations to calculate properties of electrical transport in Dirac's billiard such as conductance, power shot noise, concurrence, the entanglement of formation, and violations of Bell's inequality. All of the results of this repository can be found in the jupyter notebook's files of the directory "/Chiral/Data_Analysis". Dirac's billiard is a type of quantum billiard used to describe a system with chiral/sub-lattice symmetry. Systems with chiral symmetry can be found in some categories of topological insulators and in single-layer graphene quantum billiard, illustrated in the Figure 1 from the thesis [[1]](#1) that produces the reference [[2]](#2).

<p align="center" width="100%">
  Figure 1<br>
    <img width="40%" src="https://user-images.githubusercontent.com/53999015/147810879-3242820a-749b-4b06-9f81-92ba5c8cebe9.png">
</p>

The difference between Schrödinger and Dirac billiards remains mainly in the electron behavior described by the Schrödinger and Dirac equations, respectively. Moreover, the dimension of hamiltonian and coupling matrix between lead and cavity are double due to the chiral/sub-lattice degree of freedom.

## Conductance

**In this repository**, I analyze the influence of the three fundamental symmetries (chGOE, chGUE, chGSE) in the conductance of electrons inside the Dirac billiard. I obtain the distributions of conductance, as well as the average and variance in terms of the number of open channels for each class of symmetry. I also perform simulations to see how the conductance of each class of symmetry depends on the barrier between the lead and chaotic cavity, which Γ is called the transparency of the barrier. Moreover, I analyze the weak localization correction of conductance as a function of the barrier's transparency Γ and compare with the Schrodinger billiard. 




## References

<a id="1">[1]</a>
Barros, M. S. M. (2018). Teoria do ruído e fenômenos de interferência quântica em nano estruturas quirais. Institutional Repository of UFPB. https://repositorio.ufpb.br/jspui/handle/123456789/14825

<a id="2">[2]</a>
Barros, M. S. M., Júnior, A. J. N., Macedo-Junior, A. F., Ramos, J. G. G. S., & Barbosa, A. L. R. (2013). Open chaotic Dirac billiards: Weak (anti)localization, conductance fluctuations, and decoherence. Physical Review B, 88(24). doi:10.1103/physrevb.88.245133.

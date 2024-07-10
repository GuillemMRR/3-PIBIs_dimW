Supplementary code for:
Guillem Müller-Rigat, Albert Aloy, Maciej Lewenstein, Matteo Fadel and Jordi Tura, Three-outcome multipartite Bell inequalities: applications to dimension witnessing and spin-nematic squeezing in many-body systems, 
[arXiv.2406.12823](https://doi.org/10.48550/arXiv.2406.12823).

### dimW
The folder contains the routines for the computation of the dimension-constrainted bound. 
- [bounds_qubit_Q.py](dimW/bounds_qubit_Q.py)  variational bound and quantum bound. For the last, we employ [ncpol2sdpa](https://ncpol2sdpa.readthedocs.io/en/stable/). 
- [PIBIout3qubits.m](/dimW/PIBIout3qubits.m), [PIBIout3qubits.m](/dimW/PIBINV_test.m): applies the SDP relaxation according to the NV hierarchy supported in Matlab's package [qdimsum](https://denisrosset.github.io/qdimsum/).   
- [bounds_comparison.nb](/dimW/bounds_comparison.nb) compares the bound with the maximal violation to construct Fig. 2 of the paper.   

### appW

The folder contains the verification of the witnesses and its application to the spin-1 Bose-Einstein condensate experiment. 

- [type2.nb](appW/type2.nb) addresses effective two-level witness Eq. 23  and Fig. 4 of the paper.
- [type1.nb](appW/type1.nb) addresses three-level witness Eq. 28  and Fig. 5 of the paper.

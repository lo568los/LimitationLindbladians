## Introduction
This repository contains code used for writing the paper "Rigorous results on quantum master equations preserving complete positivity and
obeying local conservation laws" [arxivlink]. The naming convention in the code closely matches the notation used in the paper. For any questions / comments, please feel free to contact me. I will try my best to get back to you. 

## Requirements
1. Matlab with CVX http://cvxr.com/cvx/ installed correctly. Mosek solver is preferred.
2. Python with Qutip, Numpy, Scipy installed.

## Code
1. **redfield_gamma_construction** : This python script computes the matrix $\Gamma$ for a setup of a two-site open XXZ spin chain, with the first site coupled to the bath (Appendix B of the paper).
2. **minimize_tauopt**, **minimize_tauopt_function** : These MATLAB scripts use CVX to compute the minimum value of $\tau$, subject to certain constraints (Section IV of the paper)
3. **thermalization_matlab_interface** : This python script reads from the .mat file that stores the optimal $\Gamma^{(L)}$ and $H^{(L)}_{LS}$ computed by CVX, and then constructs the quantum master equation (QME) from these values. It then computes the steady-state of the QME and eigenvalues of its Liouvillian operator. The trace distance of the computed steady-state with the Gibbs-state is a measure of thermalization, while the eigenvalues of the Liouvillian can be used to check whether the uniqueness of the steady-state. This can be used for sanity checks / other cross checks / or to study the dynamics of the QME obtained from CVX. 
4.  **compute_tau**, and **compute_tau_singepoint** : These MATLAB scripts read from .mat files, and compute tau for a modified set of parameters(Section IV of the paper)


The rest of the code is self-explanatory, and consists mainly of plotting functions / helper functions. **/python/helper_code_qutip.py** is a python file consists of helper functions written using QuTiP. **/matlab/helper_functions** is a folder storing helper functions for MATLAB (this folder must be added to the MATLAB PATH). 
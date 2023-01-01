## Introduction
This repository contains code used for writing the paper "Searching Lindbladians that obey local conservation laws and show thermalization" [arxivlink]. The naming convention in the code closely matches the notation used in the paper. For any questions / comments, please feel free to contact me. I will try my best to get back to you. 

## Requirements
1. Matlab with CVX (http://cvxr.com/cvx/) installed correctly. Mosek solver is preferred.
2. Python with Qutip, Numpy, Scipy installed.

## Code
1. **redfield_gamma_construction.py** : This python script computes the matrix $\Gamma$ for a setup of a two-site open XXZ spin chain for a set of parameters, with the first site coupled to the bath (Appendix B of the paper).
2. **minimize_tauopt.m**, : This MATLAB script use CVX to compute the minimum value of $\tau$, subject to certain constraints (Section IV of the paper). The entire workspace at the end of the computation can be stored in a file named "thermal_data.mat". This .mat file can then be used to recover the $ \Gamma^{(L)}_{LS} $ and $ H^{(L)}_{LS} $ to be used in computing plots of $\tau$ vs some parameter that is varied in the bulk of the system (Fig. 6 of the paper). 
3.  **minimize_tauopt_function.m** : MATLAB function that uses CVX to compute the minimum value ot $\tau$, subject to certain constraints. Can be called for plotting values of $\tau_{opt}$ vs any  parameter.
4.  **compute_tau.m**, and **compute_tau_singepoint.m** : These MATLAB scripts read from .mat files, and compute tau for a modified set of parameters in the bulk. (Fig 6 of the paper). 
5.  **plotting_tauopt.m** computes and plots $\tau_{opt}$ vs $g$, and saves the data in a specified folder. The folder contains various plotting scripts, which can used to plot the data (Fig. 2, Fig. 4)
6.   **plotting_tauopt_vs_beta.m** computes and plots $\tau_{opt}$ vs $\beta$, and saves the data in a specified folder. The folder contains various plotting scripts, which can used to plot the data (Fig. 3, Fig. 5)
7.  **thermalization_matlab_interface.py** : This python script reads from the .mat file that stores the optimal $\Gamma^{(L)}$ and $H^{(L)}_{LS}$ computed by CVX, and then constructs the quantum master equation (QME) from these values. It then computes the steady-state of the QME and eigenvalues of its Liouvillian operator. The trace distance of the computed steady-state with the Gibbs-state is a measure of thermalization, while the eigenvalues of the Liouvillian can be used to check whether the uniqueness of the steady-state. This can be used for sanity checks / other cross checks / or to study the dynamics of the QME obtained from CVX. 


The rest of the code is self-explanatory, and consists mainly of plotting functions and helper functions.  

 **/python/helper_code_qutip.py** is a python file consists of helper functions written using QuTiP.  
  **/matlab/helper_functions** is a folder storing helper functions for MATLAB (this folder must be added to the MATLAB PATH).  
All plotting scripts are located in the same folder as the data they plot.

## Introduction
This repository contains code used for writing the paper "Semidefinite Programming for understanding limitations of Lindblad Equations" [arxivlink](https://arxiv.org/abs/2301.02146). The naming convention in the code closely matches the notation used in the paper. For any questions / comments, please feel free to contact me. I will try my best to get back to you. 

## Requirements
1. Matlab with CVX (http://cvxr.com/cvx/) installed correctly. Mosek solver is preferred.
2. Python with Qutip, Numpy, Scipy installed.

## Code
1. **final_data_beta.py, final_data_g.py** : This python script computes the exact zeroth-order non-equilibrium steady state (NESS) diagonal elements from the Redfield formalism by solving linear equations [cf. Eq. (12) in the paper] for a setup of a two-site open XXZ/XX spin chain for a given set of parameters ($e,~g,~\beta_L,~\beta_R$), with the first $N_L$ and last $N_R$ sites coupled to baths. The additional .sub and .sh bash files with the same names are for running the codes in parallel on a cluster.
2. **minimize_tauopt_ness_g.m, minimize_tauopt_ness_beta.m**, : These MATLAB scripts use CVX to compute the minimum value of $\tau$, subject to certain constraints [cf. Eq. (25) of the paper]. The entire workspace at the end of the computation can be stored in a file named "diag_data...". This .mat file can then be used to recover the $\Gamma^{(L,R)}_{LS}$ and $H^{(L,R)}\_{LS}$  to be used in computing plots of $\tau_{\rm opt}$ vs $g$ or $\beta_L$ that is varied in the bulk of the system (Panels (a),(c) of Figs. 2, 3, A1, A2 of the paper) 
3. **minimize_taucoh_ness_g.m, minimize_taucoh_ness_beta.m** : These MATLAB scripts use CVX to compute the minimum value of $\tau^{\rm coh}$, subject to certain constraints [cf. Eq. (28) of the paper]. The entire workspace at the end of the computation can be stored in a file named "coh_data...". used in computing plots of $\tau$ vs $g$ or $\beta_L$ that is varied in the bulk of the system (Panels (b),(d) of Figs. 2, 3, A1, A2 of the paper). By setting $N_R = 0$, we can recover Fig. 4 of the paper as well.
4. **minimize_tauopt_function.m, minimize_tauopt_function2.m** : MATLAB functions that use CVX to compute the minimum value ot $\tau$ and $\tau^{\rm coh}$ respectively, subject to certain constraints. Can be called for plotting values of $\tau_{opt}$ and $\tau^{\rm coh}_{\rm opt}$ vs any parameter.
5. **ness_data**: This folder contains all the data used for plotting the optimal $\tau$ values for populations and coherences. Obtained from the python files highlighted in 1.


The rest of the code is self-explanatory, and consists mainly of plotting functions and helper functions.  

 **/python/helper_code_qutip.py** is a python file consists of helper functions written using QuTiP.  
  **/matlab/helper_functions** is a folder storing helper functions for MATLAB (this folder must be added to the MATLAB PATH).  
All plotting scripts are located in the same folder as the data they plot.

from qutip import *
import numpy as np
from scipy import integrate
from helper_code_qutip import * 
import sys
import scipy.io

beta_r = float(sys.argv[1])
#g = float(sys.argv[2])
s = float(sys.argv[2])
beta_l = float(sys.argv[3])
e = float(sys.argv[4])

ham_type = 1
g = 0.01


def optimized_L2_red(eigstates, C1, D1, C2, D2, create_sm_list_left, create_sm_list_right, dims):
    """
    Optimized function to construct the L2_red superoperator, corrected for summation over all operators.
    
    Args:
        eigstates (list): List of eigenstate Qobj vectors.
        C1, C2, D1, D2 (np.ndarray): Matrices of constants.
        create_sm_list_left (list): List of left-side creation/destruction operators.
        create_sm_list_right (list): List of right-side creation/destruction operators.
        dims (list): Dimensions of the system.

    Returns:
        Qobj: The L2_red superoperator.
    """
    superop_total = 0
    N = len(eigstates)

    # Combine all operator lists for easier looping
    all_op_lists = [(create_sm_list_left, C1, D1), (create_sm_list_right, C2, D2)]

    # Loop through all eigenstates for the alpha and gamma indices
    for alpha in range(N):
        for gamma in range(N):
            E_alpha = eigstates[alpha] * eigstates[alpha].dag()
            E_gamma = eigstates[gamma] * eigstates[gamma].dag()
            
            # Loop through both the left and right operator sets
            for op_list, C_matrix, D_matrix in all_op_lists:  #left or right baths
                for op in op_list:
                    # Terms corresponding to C and C*
                    term1 = (spre(op.dag())*spost(E_alpha * op * E_gamma) - spost(E_alpha * op * E_gamma * op.dag())) * C_matrix[alpha, gamma] #since keeping epsilon^2 positive
                    term2 = (spre(E_alpha * op * E_gamma) * spost(op.dag()) - spre(op.dag()*E_alpha * op * E_gamma)) * D_matrix[alpha, gamma]
                    
                    # Add terms and their Hermitian conjugates
                    superop_total += term1 + term1.dag() + term2 + term2.dag()

    # Re-normalize dimensions and return
    superop_total.dims = [[dims, dims], [dims, dims]]
    return superop_total


#matlab_data_g = scipy.io.loadmat(f'../matlab/data_plotting_vss/coh_data_new2_NL1={NL1},e={e:.2f},beta_r={beta_r:.1f},g={g:.4f},ham_type=1.mat',mat_dtype=False)



NL1 = 1
NL = NL1
NR = NL1
NM = 2
N = NL + NR + NM
dL = 2**NL
dR = 2**NR
dM = 2**NM
d = dL * dR * dM
dims = [2]*N

create_sm_list_left = [create_sm(N,i + 1) for i in range(NL)]
create_sm_list_right = [create_sm(N,NM + NL + i + 1) for i in range(NR)]

epsilon = 0.01 #system bath coupling strength


delta = 1
beta1 = beta_r  #right baths
beta2 = beta_l

#gvals = np.logspace(-2,1,12)
w0list = [2.0] * N  #system frequencies
glist = [-g] * (N - 1)  #coupling strengths


gamma1 = 1
gamma2 = 1
limit_value = 700
b_val = 50
mu1 = -1e-10
mu2 = -1e-10
epsilon = 0.01
tb = 0.01

H_S = create_hamiltonian(w0list,glist,delta,N)
eigenenergies, eigenstates = H_S.eigenstates()
number = len(eigenenergies) # should be 2^N


print("s = ", s, "NL = ",NL1)

integral11=np.empty((number,number),dtype=np.cdouble) #stores J * N integral for left bath
integral12=np.empty((number,number),dtype=np.cdouble) # stores J integral (just to check) for the left bath
integral21=np.empty((number,number),dtype=np.cdouble) #stores J*N integral for right bath
integral22=np.empty((number,number),dtype=np.cdouble)

        #print("Integral calculations at beta2 = {} and e = {} are : ".format(beta2,e))

for i in range(number):
    for k in range(number):
                freq=eigenenergies[k]-eigenenergies[i]
                #print(f"Absolute frequency  for i = {i}, k = {k} is ",np.absolute(freq))
                #print(i,k,freq)
                if( np.absolute(freq) >= 1/10**10):
                    integral11[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func1,0,b_val,args=(s,tb,beta2,mu2,gamma1),limit=limit_value,weight='cauchy',wvar=eigenenergies[k]-eigenenergies[i])[0] #func 1
                    integral12[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath,0,b_val,args=(s,tb,gamma1),limit=limit_value,weight='cauchy',wvar=eigenenergies[k]-eigenenergies[i])[0]  #left bath done
                    integral21[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func1,0,b_val,args=(s,tb,beta1,mu1,gamma2),limit=limit_value,weight='cauchy',wvar=eigenenergies[k]-eigenenergies[i])[0] #func 1
                    integral22[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath,0,b_val,args=(s,tb,gamma2),limit=limit_value,weight='cauchy',wvar=eigenenergies[k]-eigenenergies[i])[0]  #right bath
        
                if (np.absolute(freq)<=1/10**10):  #The problem is arising here....
                    integral11[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func2,0,b_val,args=(s,tb,beta2,mu2,gamma1),limit=limit_value)[0]
                    integral12[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath_2,0,b_val,args=(s,tb,gamma1),limit=limit_value)[0]
                    integral21[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func2,0,b_val,args=(s,tb,beta1,mu1,gamma2),limit=limit_value)[0]
                    integral22[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath_2,0,b_val,args=(s,tb,gamma2),limit=limit_value)[0]
                
            
            #expected=1.0j*(eigenenergies[k]-eigenenergies[i])/(2*tb*tb)
        #        print(i,k,integral2[i,k],expected)
    
    
        # PAY ATTENTION TO THE WAY THESE COEFFICIENTS ARE BEING COMPUTED
    
constant12=np.empty((number,number),dtype=np.cdouble)
constant11=np.empty((number,number),dtype=np.cdouble)
constant21=np.empty((number,number),dtype=np.cdouble)
constant22=np.empty((number,number),dtype=np.cdouble)
        
        
        
for i in range(number):
        for k in range(number):
                freq = eigenenergies[k]-eigenenergies[i]
                if np.abs(freq) >= 1/10**10:
                    constant12[i,k]=(integral12[i,k]+integral11[i,k])+(0.5*(spectral_bath(eigenenergies[k]-eigenenergies[i],s,tb,gamma1))+func1(eigenenergies[k]-eigenenergies[i],s,tb,beta2,mu2,gamma1))/(np.abs(freq)**(s-1))  #full coefficient created this is nbar+1
                    constant11[i,k]=integral11[i,k]+(0.5*func1(eigenenergies[k]-eigenenergies[i],s,tb,beta2,mu2,gamma1))/(np.abs(freq)**(s-1))                                  # the full coefficient is created
                    
                    constant22[i,k]=(integral22[i,k]+integral21[i,k])+(0.5*(spectral_bath(eigenenergies[k]-eigenenergies[i],s,tb,gamma2)+func1(eigenenergies[k]-eigenenergies[i],s,tb,beta1,mu1,gamma2)))/(np.abs(freq)**(s-1))   #full coefficient created this is nbar+1
                    constant21[i,k]=integral21[i,k]+0.5*func1(eigenenergies[k]-eigenenergies[i],s,tb,beta1,mu1,gamma2)/(np.abs(freq)**(s-1))

"""l2_red = optimized_L2_red(eigenstates, constant11,constant12, constant21, constant22, create_sm_list_left, create_sm_list_right, dims)

l0 = liouvillian(H_S)
l_total = l0 + epsilon**2*l2_red

rho_red = steadystate(l_total)
c_1 = create_sm(N,1)
c_N = create_sm(N,N)

"""


pre=-1.0j*H_S
post=1.0j*H_S

L=spre(pre)+spost(post)

for i in range(number):
    for k in range(number):
        vi=eigenstates[i]
        vk=eigenstates[k]

        print(constant11[i,k],constant12[i,k])
        for c_1 in create_sm_list_left:
        
            op1=epsilon*epsilon*constant11[i,k]*vi*vi.dag()*c_1*vk*vk.dag()*c_1.dag()
            op2=epsilon*epsilon*constant12[i,k]*c_1.dag()*vi*vi.dag()*c_1*vk*vk.dag()
            
            op3=epsilon*epsilon*constant11[i,k]*c_1.dag()
            op4=vi*vi.dag()*c_1*vk*vk.dag()
            op5=epsilon*epsilon*constant12[i,k]*c_1.dag()
            
            
            L=L+spre(-op2-op1.dag())+spost(-op1-op2.dag())
            L=L+spre(op3)*spost(op4)+spre(op4)*spost(op5)+spre(op4.dag())*spost(op3.dag()) +spre(op5.dag())*spost(op4.dag())

        for c_N in create_sm_list_right:
        
            op1=epsilon*epsilon*constant21[i,k]*vi*vi.dag()*c_N*vk*vk.dag()*c_N.dag()
            op2=epsilon*epsilon*constant22[i,k]*c_N.dag()*vi*vi.dag()*c_N*vk*vk.dag()
            
            op3=epsilon*epsilon*constant21[i,k]*c_N.dag()
            op4=vi*vi.dag()*c_N*vk*vk.dag()
            op5=epsilon*epsilon*constant22[i,k]*c_N.dag()
            
            
            L=L+spre(-op2-op1.dag())+spost(-op1-op2.dag())
            L=L+spre(op3)*spost(op4)+spre(op4)*spost(op5)+spre(op4.dag())*spost(op3.dag()) +spre(op5.dag())*spost(op4.dag())
        
        
        
#Variables needed for for iterative-lgmres to work. 
return_info=True
#print('Redfield Liouvillian constructed, Computing steady-state ...')
ss_redfield = steadystate(L,return_info=return_info)
L_eigen = L.eigenenergies()
print("Smallest eigenvalues of L_red are ", L_eigen[-3:])
#list_red.append(ss_redfield)
np.savetxt(f"rho_red_s={s:.2f}_NL1={NL1}_e={e:.2f}_beta_r={beta_r:.1f}_g={g:.4f}.txt", ss_redfield.full())

#print("Smallest eigenvalues of Ltotal are ", l_total.eigenenergies()[-3:])

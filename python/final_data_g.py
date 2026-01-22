from qutip import *
import numpy as np
from scipy import integrate
from helper_code_qutip import * 
import scipy.io
import sys

beta_r = float(sys.argv[1])
g = float(sys.argv[2])
ham_type = int(sys.argv[3])
beta_l = float(sys.argv[4])
e = float(sys.argv[5])

NL1 = 3
NL2 = 0
NM = 0

N = NL1 + NL2 + NM
dL1 = 2**NL1
dL2 = 2**NL2
dM = 2**NM
d = 2**N
dims = [2]*N

create_sm_list_left = [create_sm(N,i + 1) for i in range(NL1)]
create_sm_list_right = [create_sm(N,NM + NL1 + i + 1) for i in range(NL2)]

## Firstly, we have to define the integral function C and D

def integral1(i,k,tb,beta,mu,gamma,eigenergies,limit_value = 700,b=50):
    freq=eigenergies[k]-eigenergies[i]
    if( np.absolute(freq) >= 1/10**10):
        integral = (-1.0j/(2*np.pi))*integrate.quad(func1,0,b,args=(tb,beta,mu,gamma),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0]
    else:
        integral = (-1.0j/(2*np.pi))*integrate.quad(func2,0,b,args=(tb,beta,mu,gamma),limit=limit_value)[0]
    return integral

def integral2(i,k,tb,gamma,eigenergies,limit_value = 700,b=50):
    freq=eigenergies[k]-eigenergies[i]
    if( np.absolute(freq) >= 1/10**10):
        integral = (-1.0j/(2*np.pi))*integrate.quad(spectral_bath,0,b,args=(tb,gamma),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0]
    else:
        integral = (-1.0j/(2*np.pi))*integrate.quad(spectral_bath_2,0,b,args=(tb,gamma),limit=limit_value)[0]
    return integral

def C(i,k,tb,beta,mu,gamma,eigenergies):
    val = integral1(i,k,tb,beta,mu,gamma,eigenergies) + 0.5*(func1(eigenergies[k]-eigenergies[i],tb,beta,mu,gamma))

    return val

def D(i,k,tb,beta,mu,gamma,eigenergies):
    val = integral1(i,k,tb,beta,mu,gamma,eigenergies) + integral2(i,k,tb,gamma,eigenergies) + 0.5*(spectral_bath(eigenergies[k]-eigenergies[i],tb,gamma)+func1(eigenergies[k]-eigenergies[i],tb,beta,mu,gamma))
    return val

def L2_red(rho,eigstates,number, constant11,constant12,constant21,constant22):
    sum = 0
    rho = Qobj(rho)
    rho.dims = [dims,dims]
    for i in range(number):
        for k in range(number):
            vi = eigstates[i]
            vk = eigstates[k]

            proj_i = vi*vi.dag()
            proj_k = vk*vk.dag()

            for l in range(NL1):
                op = commutator(rho*proj_i*create_sm_list_left[l]*proj_k,create_sm_list_left[l].dag())*constant11[i,k] + commutator(create_sm_list_left[l].dag(),proj_i*create_sm_list_left[l]*proj_k*rho)*constant12[i,k]
                sum += op
                sum += op.dag()

            for l in range(NL2):
                op = commutator(rho*proj_i*create_sm_list_right[l]*proj_k,create_sm_list_right[l].dag())*constant21[i,k] + commutator(create_sm_list_right[l].dag(),proj_i*create_sm_list_right[l]*proj_k*rho)*constant22[i,k]
                sum += op
                sum += op.dag()
    data = sum.full()
    sum = np.array(data,dtype=complex)
    return sum

def re_ness_g(beta_r,beta_l,g,ham_type,e):
    #Define the parameters
    print("Beta_r is ",beta_r, "and Beta_l is ",beta_l)
    w0list = np.linspace(1,1,N)
    for i in range(int(N/2),N):
        w0list[i] = 1 + e
   
    delta = 1
    beta1 = beta_r  #right baths
    beta2 = beta_l

    #gvals = np.logspace(-2,1,12)
    

    gamma1 = 1
    gamma2 = 1
    limit_value = 700
    b_val = 50
    mu1 = -1e-10
    mu2 = -1e-10
    epsilon = 0.01
    tb = 0.01

    s = 1

    gamma_list = [1,1,1]  #for 3 sites


  
    print("g:",g)
    glist = np.linspace(g,g,N-1)
    #print("g is ",g)
        #Define the Hamiltonian
    if ham_type == 1:
        H_S = create_hamiltonian3(w0list,glist,delta,N)
    elif ham_type == 2:
        H_S = create_hamiltonian3(w0list,glist,0,N)
    eigenergies,eigstates=H_S.eigenstates()
    number = len(eigenergies)
    #Define the change of basis unitary matrix
    U = np.zeros((number,number),dtype=complex)
    for i in range(number):
        U[:,i] = eigstates[i].full().flatten()


    #rho_th = (-beta2*H_S).expm()/((-beta2*H_S).expm()).tr() #left thermal density matrix
    #print(rho_th)
    number = len(eigenergies)
    integral11=np.empty((number,number),dtype=np.cdouble) #stores J * N integral for left bath
    integral12=np.empty((number,number),dtype=np.cdouble) # stores J integral (just to check) for the left bath
    integral21=np.empty((number,number),dtype=np.cdouble) #stores J*N integral for right bath
    integral22=np.empty((number,number),dtype=np.cdouble)

            #print("Integral calculations at beta2 = {} and e = {} are : ".format(beta2,e))

    for i in range(number):
        for k in range(number):
                    freq=eigenergies[k]-eigenergies[i]
                    #print(f"Absolute frequency  for i = {i}, k = {k} is ",np.absolute(freq))
                    #print(i,k,freq)
                    if( np.absolute(freq) >= 1/10**10):
                        integral11[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func1,0,b_val,args=(s,tb,beta2,mu2,gamma1),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0] #func 1
                        integral12[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath,0,b_val,args=(s,tb,gamma1),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0]  #left bath done
                        integral21[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func1,0,b_val,args=(s,tb,beta1,mu1,gamma2),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0] #func 1
                        integral22[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath,0,b_val,args=(s,tb,gamma2),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0]  #right bath
            
                    if (np.absolute(freq)<=1/10**10):  #The problem is arising here....
                        integral11[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func2,0,b_val,args=(s,tb,beta2,mu2,gamma1),limit=limit_value,points=[0])[0]
                        integral12[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath_2,0,b_val,args=(s,tb,gamma1),limit=limit_value,points=[0])[0]
                        integral21[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func2,0,b_val,args=(s,tb,beta1,mu1,gamma2),limit=limit_value,points=[0])[0]
                        integral22[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath_2,0,b_val,args=(s,tb,gamma2),limit=limit_value,points=[0])[0]
                    
                
                #expected=1.0j*(eigenergies[k]-eigenergies[i])/(2*tb*tb)
            #        print(i,k,integral2[i,k],expected)
        
        
            # PAY ATTENTION TO THE WAY THESE COEFFICIENTS ARE BEING COMPUTED
        
    constant12=np.empty((number,number),dtype=np.cdouble)
    constant11=np.empty((number,number),dtype=np.cdouble)
    constant21=np.empty((number,number),dtype=np.cdouble)
    constant22=np.empty((number,number),dtype=np.cdouble)
            
            
            
    for i in range(number):
            for k in range(number):
                    constant12[i,k]=integral12[i,k]+integral11[i,k]+0.5*(spectral_bath(eigenergies[k]-eigenergies[i],tb,gamma1)+func1(eigenergies[k]-eigenergies[i],tb,beta2,mu2,gamma1))    #full coefficient created this is nbar+1
                    constant11[i,k]=integral11[i,k]+0.5*func1(eigenergies[k]-eigenergies[i],s,tb,beta2,mu2,gamma1)                                       # the full coefficient is created
                    
                    constant22[i,k]=integral22[i,k]+integral21[i,k]+0.5*(spectral_bath(eigenergies[k]-eigenergies[i],tb,gamma2)+func1(eigenergies[k]-eigenergies[i],tb,beta1,mu1,gamma2))    #full coefficient created this is nbar+1
                    constant21[i,k]=integral21[i,k]+0.5*func1(eigenergies[k]-eigenergies[i],s,tb,beta1,mu1,gamma2)   # the full coefficient is created
                    #print(i,k,constant11[i,k],constant12[i,k],constant21[i,k],constant22[i,k])



    ## Now we will write out the matrix elements

    A = np.zeros((number,number),dtype=complex)

    for i in range(number):
        for k in range(number):
            sum1 = 0
            vi = eigstates[i]
            vk = eigstates[k]
            #proj_i = vi*vi.dag()
            proj_k = vk*vk.dag()
            for y in range(number):
                for l in range(NL1):
                    proj_y = eigstates[y]*eigstates[y].dag()
                    op1 = commutator(proj_k*create_sm_list_left[l]*proj_y,create_sm_list_left[l].dag())*constant11[k,y]
                    sum1 += gamma_list[l]*epsilon*epsilon*vi.dag()*(op1 + op1.dag())*vi

                    op2 = commutator(create_sm_list_left[l].dag(),proj_y*create_sm_list_left[l]*proj_k)*constant12[y,k]
                    sum1 += gamma_list[l]*epsilon*epsilon*vi.dag()*(op2 + op2.dag())*vi

                for l in range(NL2):
                    proj_y = eigstates[y]*eigstates[y].dag()
                    op1 = commutator(proj_k*create_sm_list_right[l]*proj_y,create_sm_list_right[l].dag())*constant21[k,y]
                    sum1 += gamma_list[l]*epsilon*epsilon*vi.dag()*(op1 + op1.dag())*vi

                    op2 = commutator(create_sm_list_right[l].dag(),proj_y*create_sm_list_right[l]*proj_k)*constant22[y,k]
                    sum1 += gamma_list[l]*epsilon*epsilon*vi.dag()*(op2 + op2.dag())*vi
            
    #print(np.array(sum1))
            #print(np.array(sum1)[0])
            A[i,k] = np.array(sum1)[0][0]

    b = np.zeros((number),dtype=complex)
    A_new = A[:-1]
    A_new = np.vstack([A_new,np.ones((1,number))])
    b[-1] = 1  ## Last element of b is 1

    x = np.linalg.solve(A_new,b)

    print("Correctness check:",np.dot(A_new,x))
    print(np.dot(A[-1],x))

    x_real = [np.real(x[i]) for i in range(number)]

    rho = np.diag(x_real)

    #set U matrix whose columns are the eigenvectors of the Hamiltonian

    

    rho_comp2 = np.dot(U,np.dot(rho,U.T.conjugate()))
    #rho_ness_arr.append(rho_comp2)
    l2_red = L2_red(rho_comp2,eigstates,number,constant11,constant12,constant21,constant22)
    print("L2_red shape:", l2_red.shape)

    #l2_red_arr.append(L2_redfield)

    data_dict = {"dm_ness":rho_comp2,"L2_red":l2_red, "beta2":beta2, "g":g, "e":e, "ham_type":ham_type}

    scipy.io.savemat(f'ness_data_neqall3_NL1={NL1},NL2={NL2},NM={NM},e={e:.2f},beta_r={beta_r:.1f},beta_l={beta_l:.1f},g={g:.4f}_{ham_type}.mat',data_dict)

re_ness_g(beta_r,beta_l,g,ham_type,e)



        

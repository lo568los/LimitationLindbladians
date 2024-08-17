## Here, we will implement the Redfield and Local Lindblad forms of the master equation

from qutip import *
import matplotlib.pyplot as plt
import numpy as np
from helper_code_qutip import *
from scipy import integrate
import scipy.io

def evolve(state,H_S,list1,list2,C11,C12,C21,C22,number,epsilon,c_1,c_N,indices1,indices2):
    term1=1.0j*commutator(state,H_S)
    
    term2=0

    for index in indices1:
        i=index[0]
        k=index[1]
        term2=term2+commutator(state*list1[i][k],c_1.dag())*C11[i,k]
        term2=term2+commutator(c_1.dag(),list1[i][k]*state)*C12[i,k]
        
    for index in indices2:
        i=index[0]
        k=index[1]
        term2=term2+commutator(state*list2[i][k],c_N.dag())*C21[i,k]
        term2=term2+commutator(c_N.dag(),list2[i][k]*state)*C22[i,k]


    return term1-epsilon*epsilon*(term2+term2.dag())



def convergencecheck(array):
    length=len(array)
    max=np.max(array)
    min=np.min(array)
    mean=np.mean(array)
    

    #print(max,min,"mean =",mean)
    if ((max-min) < 0.05*abs(mean) and (max-min) > -0.05*abs(mean)):
        return 1
    
    return 0


def gmatrix(omega, gamma, beta, mu, tb, i , j): #computes the i,j the element of the submatrix of the lth bath. 
    # l is determined by beta, mu and tb etc
    
    submatrix=np.zeros((2,2),dtype=np.complex128)
    submatrix[0,0]=1
    submatrix[0,1]=1.0j
    submatrix[1,0]=-1.0j
    submatrix[1,1]=1                                  
    if (omega <= 0):
        factor=np.sqrt(func1(-omega,tb,beta,mu,gamma)*2)/(8*np.pi)    
        return factor*submatrix[i-1,j-1]
    if (omega > 0):
        factor=np.sqrt(2*(func1(omega,tb,beta,mu,gamma)+spectral_bath(omega,tb,gamma)))/(8*np.pi)
        return factor*submatrix.conj()[i-1,j-1]

#declaring parameters
#

def lamb_integrand(omega,E1,E2,alpha,beta,gamma1,gamma2,beta1,beta2,mu1,mu2,tb,flag): #this function computes the integrand of the lamb shift formula.
    term=0
    
    if (alpha <0 or beta<0):
        print('indices are wrong. ERror')
        return 0
    elif (alpha>4 or beta>4):
        print('indices are wrong. err0r')
        return 0
    #alpha, beta have to be between 1 and 4.
    elif (alpha <=2 and beta>=3): #we output zero
        return 0
    elif (alpha>=3 and beta <=2):
        return 0
    elif (alpha<=2 and beta<=2): #we are in the first bath setup
        term=term+gmatrix(omega-E1,gamma1,beta1,mu1,tb,alpha,1)*gmatrix(omega+E2,gamma1,beta1,mu1,tb,1,beta)
        term=term+gmatrix(omega-E1,gamma1,beta1,mu1,tb,alpha,2)*gmatrix(omega+E2,gamma1,beta1,mu1,tb,2,beta)
        if (flag==0):
            return term.real
        elif (flag==1):
            return term.imag
        else:
            print('flag invalid')
            return term
    elif (alpha >=3 and beta>=3):
        term=term+gmatrix(omega-E1,gamma2,beta2,mu2,tb,alpha-2,1)*gmatrix(omega+E2,gamma2,beta2,mu2,tb,1,beta-2)
        term=term+gmatrix(omega-E1,gamma2,beta2,mu2,tb,alpha-2,2)*gmatrix(omega+E2,gamma2,beta2,mu2,tb,2,beta-2)
        if (flag==0):
            return term.real
        elif (flag==1):
            return term.imag
        else:
            print('flag invalid')
            return term
        
Tc_list = [0.001,0.18,0.20]
Th_list = np.linspace(0.001,8,15)


betalist1 = [1/Tc for Tc in Tc_list]
betalist2 = [1/Th for Th in Th_list]

print(len(betalist1),len(betalist2))

redfield_ss = []
lle_ss = []
g = 0
for beta1 in betalist1:
    list_red = []
    list_lle = []
    for beta2 in betalist2:
        limit_value=700
        b=50
        N=2
        w0min=1
        w0max=1
        delta=1
        gmin=16
        gmax=16
            
        w0list=np.linspace(w0min,w0max,N)
        glist=np.linspace(gmin,gmax,N-1)

        g = glist[0]
            
        tb=0.01
        epsilon=0.1
        gamma1=1 #gamma1 is the coupling to left bath. It shows up in spectral bath function
        gamma2=1    #gamma2 is the coupling to the right bath.    
            
            

        mu=-0.5
            
            
        delta=1
        mu1=mu
        mu2=mu

        H_S=create_hamiltonian2(w0list,glist,N)
    
    
        c_N=create_sm(N,N)  # we couple the Nth spin to the bath
        c_1=create_sm(N,1)
        
        
        
        eigenergies,eigstates=H_S.eigenstates()
        
        #print("eigenenergies are : ",eigenergies)
        
        spectrum=max(eigenergies)-min(eigenergies)
        
        
        
        number=len(eigenergies)
        
        integral11=np.empty((number,number),dtype=np.cdouble) #stores J * N integral for left bath
        integral12=np.empty((number,number),dtype=np.cdouble) # stores J integral (just to check) for the left bath
        integral21=np.empty((number,number),dtype=np.cdouble) #stores J*N integral for right bath
        integral22=np.empty((number,number),dtype=np.cdouble)

        for i in range(number):
            for k in range(number):
                freq=eigenergies[k]-eigenergies[i]
                #print(i,k,freq)
                if( np.absolute(freq) >= 1/10**10):
                    integral11[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func1,0,b,args=(tb,beta1,mu1,gamma1),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0] #func 1
                    integral12[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath,0,b,args=(tb,gamma1),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0]  #left bath done
                    integral21[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func1,0,b,args=(tb,beta2,mu2,gamma2),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0] #func 1
                    integral22[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath,0,b,args=(tb,gamma2),limit=limit_value,weight='cauchy',wvar=eigenergies[k]-eigenergies[i])[0]  #right bath
        
                if (np.absolute(freq)<=1/10**10):
                    integral11[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func2,0,b,args=(tb,beta1,mu1,gamma1),limit=limit_value)[0]
                    integral12[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath_2,0,b,args=(tb,gamma1),limit=limit_value)[0]
                    integral21[i,k]=(-1.0j/(2*np.pi))*integrate.quad(func2,0,b,args=(tb,beta2,mu2,gamma2),limit=limit_value)[0]
                    integral22[i,k]=(-1.0j/(2*np.pi))*integrate.quad(spectral_bath_2,0,b,args=(tb,gamma2),limit=limit_value)[0]
                
            
            #expected=1.0j*(eigenergies[k]-eigenergies[i])/(2*tb*tb)
        #        print(i,k,integral2[i,k],expected)
    
    
        # PAY ATTENTION TO THE WAY THESE COEFFICIENTS ARE BEING COMPUTED
    
        constant12=np.empty((number,number),dtype=np.cdouble)
        constant11=np.empty((number,number),dtype=np.cdouble)
        constant21=np.empty((number,number),dtype=np.cdouble)
        constant22=np.empty((number,number),dtype=np.cdouble)
        
        
        
        for i in range(number):
            for k in range(number):
                constant12[i,k]=integral12[i,k]+integral11[i,k]+0.5*(spectral_bath(eigenergies[k]-eigenergies[i],tb,gamma1)+func1(eigenergies[k]-eigenergies[i],tb,beta1,mu1,gamma1))    #full coefficient created this is nbar+1
                constant11[i,k]=integral11[i,k]+0.5*func1(eigenergies[k]-eigenergies[i],tb,beta1,mu1,gamma1)                                       # the full coefficient is created
                
                constant22[i,k]=integral22[i,k]+integral21[i,k]+0.5*(spectral_bath(eigenergies[k]-eigenergies[i],tb,gamma2)+func1(eigenergies[k]-eigenergies[i],tb,beta2,mu2,gamma2))    #full coefficient created this is nbar+1
                constant21[i,k]=integral21[i,k]+0.5*func1(eigenergies[k]-eigenergies[i],tb,beta2,mu2,gamma2)   # the full coefficient is created
                #print(i,k,constant11[i,k],constant12[i,k],constant21[i,k],constant22[i,k])
        list1=[]
        list2=[]
        
        
        for i in range(number):
            list1.append([])
            list2.append([])
        
        
        
        matrix=np.zeros((number,number))
        
        dim=[]
        for k in range(N):
            dim.append(2)    
        
        zeromatrix=Qobj(matrix,dims=[dim,dim])
        
        
        indices1=[]
        indices2=[]
        
        
        
        
        for i in range(number):
            for k in range(number):
                list1[i].append(eigstates[i]*eigstates[i].dag()*c_1*eigstates[k]*eigstates[k].dag())
                list2[i].append(eigstates[i]*eigstates[i].dag()*c_N*eigstates[k]*eigstates[k].dag())
                
                if(tracedist(eigstates[i]*eigstates[i].dag()*c_1*eigstates[k]*eigstates[k].dag(),zeromatrix)!=0):
                    indices1.append((i,k))
                if(tracedist(eigstates[i]*eigstates[i].dag()*c_N*eigstates[k]*eigstates[k].dag(),zeromatrix)!=0):
                    indices2.append((i,k))

        pre=-1.0j*H_S
        post=1.0j*H_S
        
        L=spre(pre)+spost(post)
        
        for i in range(number):
            for k in range(number):
                vi=eigstates[i]
                vk=eigstates[k]
                
                op1=epsilon*epsilon*constant11[i,k]*vi*vi.dag()*c_1*vk*vk.dag()*c_1.dag()
                op2=epsilon*epsilon*constant12[i,k]*c_1.dag()*vi*vi.dag()*c_1*vk*vk.dag()
                
                op3=epsilon*epsilon*constant11[i,k]*c_1.dag()
                op4=vi*vi.dag()*c_1*vk*vk.dag()
                op5=epsilon*epsilon*constant12[i,k]*c_1.dag()
                
                
                L=L+spre(-op2-op1.dag())+spost(-op1-op2.dag())
                L=L+spre(op3)*spost(op4)+spre(op4)*spost(op5)+spre(op4.dag())*spost(op3.dag()) +spre(op5.dag())*spost(op4.dag())
                
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
        list_red.append(ss_redfield)
        #print("List1 length:", len(list_red))

        #print('Redfield steady-state computed for Tc = ',1/beta1,' and Th = ',1/beta2)


        ################## Local lindbald shit #############################
    
    
        Delta1=(-1.0*epsilon*epsilon/(2*np.pi))*integrate.quad(spectral_bath,0,b,args=(tb,gamma1),weight='cauchy',wvar=w0list[0])[0] #Delta
        Deltadash1=(-1.0*epsilon*epsilon/(2*np.pi))*integrate.quad(func1,0,b,args=(tb,beta1,mu1,gamma1),weight='cauchy',wvar=w0list[0])[0] #Delta
        
        
        DeltaN=(-1.0*epsilon*epsilon/(2*np.pi))*integrate.quad(spectral_bath,0,b,args=(tb,gamma2),weight='cauchy',wvar=w0list[N-1])[0] #Delta
        DeltadashN=(-1.0*epsilon*epsilon/(2*np.pi))*integrate.quad(func1,0,b,args=(tb,beta2,mu2,gamma2),weight='cauchy',wvar=w0list[N-1])[0] #Delta
        
        
        H=H_S+(Deltadash1+0.5*Delta1)*create_sigmaz(N,1)+(DeltadashN+0.5*DeltaN)*create_sigmaz(N,N)
        
        
        Cops=[]
        
        Cops.append(epsilon*np.sqrt(func1(w0list[0],tb,beta1,mu1,gamma1)+spectral_bath(w0list[0], tb, gamma1))*create_sm(N,1))
        Cops.append(epsilon*np.sqrt(func1(w0list[0],tb,beta1,mu1,gamma1))*create_sm(N, 1).dag())
        Cops.append(epsilon*np.sqrt(func1(w0list[N-1],tb,beta2,mu2,gamma2)+spectral_bath(w0list[N-1],tb,gamma2))*create_sm(N,N));
        Cops.append(epsilon*np.sqrt(func1(w0list[N-1],tb,beta2,mu2,gamma2))*create_sm(N,N).dag())
        
        
        #print('Local-Lindblad Liouvillian constructed, Computing steady-state ..')
        ss_lindblad=steadystate(H,Cops,return_info=return_info)

        list_lle.append(ss_lindblad)
        #print("List2 length:", len(list_lle))

        #print('Local-Lindblad steady-state computed for Tc = ',1/beta1,' and Th = ',1/beta2)

    
    #print("List1 length:", len(list_red))
    #print("List1: ",list_red)

    redfield_ss.append(list_red)
    lle_ss.append(list_lle)

#print('All steady-states computed!')

#print(lle_ss[0][0])

def concurrence_plot(Th_list,Tc_list,reduced_dm_list):
    for i in range(len(Tc_list)):
        concurrence_list = []
        for j in range(len(Th_list)):
            concurrence_list.append(concurrence(reduced_dm_list[i][j]))
            print(f'Concurrence Qutip for Tc = {Tc_list[i]} and Th = {round(Th_list[j],2)}: ',concurrence(reduced_dm_list[i][j]))
        plt.plot(Th_list,concurrence_list,label='Tc/E = '+str(Tc_list[i]))

    plt.xlabel('Thermal bath temperature (T_h/E)')
    plt.ylabel('Concurrence')
    plt.title(f'Concurrence vs Th for different Tc at g = {g}, pc = 1, ph = 1')
    plt.legend()

    plt.show()
print("Concurrence plot for Redfield ")
concurrence_plot(Th_list,Tc_list,redfield_ss)
print()
print("Concurrence plot for Local Lindblad ")
concurrence_plot(Th_list,Tc_list,lle_ss)

print("density matrix for Redfield: ",redfield_ss[0])
print("density matrix for LLe: ",lle_ss[0])







from qutip import *
import matplotlib.pyplot as plt
import numpy as np




def create_vacuum(N):
    psi=basis(2,1)
    for k in range(2,N+1):
        psi=tensor(psi,basis(2,1))

    return psi

def create_allup(N):
    psi=basis(2,0)
    for k in range(2,N+1):
        psi=tensor(psi,basis(2,0))    

    return psi



def create_magnetization(N):
    op=0

    for k in range(1,N+1):
        op=op+create_sigmaz(N,k)

    return op



def create_sm(N,pos): #creates the sigma_minus operator for the given position. N>=2
    if pos==1:
        op=fock(2,1)*fock(2,0).dag()
        for k in range(2,N+1):
            op=tensor(op,qeye(2))

    else:
        op=qeye(2)
        for k in range(2,N+1):
            if k==pos:
                op=tensor(op,fock(2,1)*fock(2,0).dag())
            else:
                op=tensor(op,qeye(2))

    return op


def create_sigmax(N,pos):
    if pos==1:
        op=sigmax()
        for k in range(2,N+1):
            op=tensor(op,qeye(2))

    else:
        op=qeye(2)
        for k in range(2,N+1):
            if k==pos:
                op=tensor(op,sigmax())
            else:
                op=tensor(op,qeye(2))
    return op


def create_sigmay(N,pos):
    if pos==1:
        op=sigmay()
        for k in range(2,N+1):
            op=tensor(op,qeye(2))

    else:
        op=qeye(2)
        for k in range(2,N+1):
            if k==pos:
                op=tensor(op,sigmay())
            else:
                op=tensor(op,qeye(2))
    return op


def create_sigmaz(N,pos):
    if pos==1:
        op=sigmaz()
        for k in range(2,N+1):
            op=tensor(op,qeye(2))

    else:
        op=qeye(2)
        for k in range(2,N+1):
            if k==pos:
                op=tensor(op,sigmaz())
            else:
                op=tensor(op,qeye(2))

    return op

def create_projector0(N,pos):
    if pos==1:
        op=fock_dm(2,0)
        for k in range(2,N+1):
            op=tensor(op,qeye(2))

    else:
        op=qeye(2)
        for k in range(2,N+1):
            if k==pos:
                op=tensor(op,fock_dm(2,0))
            else:
                op=tensor(op,qeye(2))

    return op

def create_projector1(N,pos):
    if pos==1:
        op=fock_dm(2,1)
        for k in range(2,N+1):
            op=tensor(op,qeye(2))

    else:
        op=qeye(2)
        for k in range(2,N+1):
            if k==pos:
                op=tensor(op,fock_dm(2,1))
            else:
                op=tensor(op,qeye(2))

    return op

def create_projector01(N,pos):
    if pos==1:
        op=tensor(fock(2,0),fock(2,1))*(tensor(fock(2,1),fock(2,0)).dag())
        for k in range(2,N):
            op=tensor(op,qeye(2))

    else:
        op=qeye(2)
        for k in range(2,N):
            if k==pos:
                op=tensor(op,tensor(fock(2,0),fock(2,1))*(tensor(fock(2,1),fock(2,0)).dag()))
            else:
                op=tensor(op,qeye(2))

    return op



def create_hamiltonian(w0list,glist,delta,N):
    
    H=(w0list[N-1]/2)*create_sigmaz(N,N)

    for k in range(1,N):
        H=H+(w0list[k-1]/2)*create_sigmaz(N,k) - glist[k-1]*(create_sigmax(N,k)*create_sigmax(N,k+1) + create_sigmay(N,k)*create_sigmay(N,k+1) + delta*create_sigmaz(N,k)*create_sigmaz(N,k+1))

    return H

def create_hamiltonian2(w0list,glist,N):
    
    H=(w0list[N-1])*create_projector0(N,N)

    for k in range(1,N):
        H=H+(w0list[k-1])*create_projector0(N,k) + (2*glist[k-1])*(create_projector01(N,k) + create_projector01(N,k).dag())

    return H

def create_hamiltonian3(w0list,glist,delta,N):
    
    H=(w0list[N-1])*create_projector0(N,N)

    for k in range(1,N):
        H=H+(w0list[k-1])*create_projector0(N,k) + glist[k-1]*(create_sigmax(N,k)*create_sigmax(N,k+1) + create_sigmay(N,k)*create_sigmay(N,k+1) + delta*create_sigmaz(N,k)*create_sigmaz(N,k+1))

    return H


def create_hamiltonian_v2(w0list,glist,deltalist,N): #includes deltalist
    
    H=(w0list[N-1]/2)*create_sigmaz(N,N)

    for k in range(1,N):
        H=H+(w0list[k-1]/2)*create_sigmaz(N,k) - glist[k-1]*(create_sigmax(N,k)*create_sigmax(N,k+1) + create_sigmay(N,k)*create_sigmay(N,k+1) + deltalist[k-1]*create_sigmaz(N,k)*create_sigmaz(N,k+1))

    return H



def spectral_bath(omega,s,tb,gamma=1):
    if(omega <=0):
        return 0

    return gamma*(omega**s)*np.exp(-omega*omega*tb)



def spectral_bath_2(omega,s,tb,gamma=1):
    if(omega <=0):
        return 0
    return gamma*(omega**(s-1))*np.exp(-omega*omega*tb)
    



def nbar(omega,beta,mu):
    x = beta * (omega - mu)
    
    # Use absolute value to check for smallness. 
    # 1e-6 is usually sufficient for switching to Taylor expansion.
    if np.abs(x) < 1e-6:
        # Handle the exact singularity to prevent ZeroDivisionError
        if x == 0:
            # Return a large finite number (cutoff) representing divergence
            # Adjust this magnitude based on your physics needs/units
            return 1e15 
        
        # Taylor expansion of 1/(exp(x)-1) around 0:
        # 1/x - 1/2 + x/12 ...
        # Adding the -0.5 term ensures better continuity at the switch point
        return (1.0 / x) - 0.5

    # Standard calculation for values away from singularity
    try:
        return 1.0 / (np.exp(x) - 1)
    except OverflowError:
        # If exp(x) overflows (very large x), the occupation is effectively 0
        return 0.0


def func1(omega,s,tb,beta,mu,gamma=1):
    if(omega <=0):
        return 0

    return spectral_bath(omega,s,tb,gamma)*nbar(omega,beta,mu)



def func2(omega,s,tb,beta,mu,gamma=1):
    if(omega<=0):
        return 0

    return spectral_bath_2(omega,s,tb,gamma)*nbar(omega,beta,mu)




def GramSchmidt(matrix_list): 
    #accepts as input a list of matrices, and returns an orthonormal basis obtained via gram schmidting
    #we assume input has identity at the start, and therefore the output list has identity at the end for convenience.

    
    length = len(matrix_list)
    output_list = []

    for index1 in range(length):
        oper = matrix_list[index1]

        for index2 in range(index1):
            numer = hilbert_schmidt_innerproduct(output_list[index2],oper)
            denom = hilbert_schmidt_innerproduct(output_list[index2],output_list[index2])

            oper = oper - (numer/denom)*output_list[index2]
        
        output_list.append(oper)

    for index in range(len(output_list)): 
        output_list[index] = output_list[index] / ( np.sqrt(hilbert_schmidt_innerproduct(output_list[index],output_list[index]) ) ) 

    if not basis_is_orthonormal(output_list):
        print("WARNING : basis is NOT orthonormal")
    output_list.reverse()
    return output_list


def hilbert_schmidt_innerproduct(oper1,oper2):
    return (oper1.dag()*oper2).tr()
        
def basis_is_orthonormal(list): #accepts a list of operators and outputs true if they are orthonormal
    tol = 1e-8
    length = len(list)
    
    for i in range(length):
        for k in range(length):
            operi = list[i]
            operk = list[k]
            innerproduct = hilbert_schmidt_innerproduct(operi,operk)
            if i==k:
                if np.absolute(innerproduct-1) > tol:
                    return False
            else:
                if np.absolute(innerproduct) > tol:
                    return False
    
    return True

def create_basis_startingwith_identity(dim):
    #creates a list of basis operators starting with identity.
    output_list = []

    for i in range(dim):
        for k in range(dim):
            matrix = np.zeros((dim,dim),dtype = np.cdouble)
            matrix[i,k] = 1
            output_list.append(Qobj(matrix))
    
    output_list[0] = Qobj(np.identity(dim)) #make sure 1st is identity.

    return output_list

                

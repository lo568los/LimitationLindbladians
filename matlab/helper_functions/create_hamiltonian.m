% creates a hamiltonian 

function H=create_hamiltonian(w0list,glist,delta,N)
    H=0.5*w0list(N)*create_sigmaz(N,N);
    
    for k=1:N-1
        H=H+0.5*w0list(k)*create_sigmaz(N,k)-glist(k)*(create_sigmax(N,k)*create_sigmax(N,k+1)+create_sigmay(N,k)*create_sigmay(N,k+1)+delta(k)*create_sigmaz(N,k)*create_sigmaz(N,k+1));
        
    end
end


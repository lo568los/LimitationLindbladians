% creates a hamiltonian 2

function H=create_hamiltonian2(w0list,glist,N)
    H= w0list(N)*create_projector0(N,N); %0 
    
    for k=1:N-1
        H=H+ w0list(k)*create_projector0(N,k) + glist(k)*(create_projector01(N,k) + create_projector10(N,k));
        
    end
end
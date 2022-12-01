% creates magnetization operator

function output=create_magnetization(N)
    output=0;
    
    for k=1:N
        output=output+create_sigmaz(N,k);
    end
    
    
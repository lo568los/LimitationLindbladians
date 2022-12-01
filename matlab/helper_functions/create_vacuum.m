function output=create_vacuum(N)
    vac=[0,0;0,1];
    output=vac;
    for k=2:N
        output=kron(output,vac);
    end
    
end

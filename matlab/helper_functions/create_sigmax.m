function output=create_sigmax(N,pos)
    if pos==1
        output=sigmax();
        for k=2:N
            output=kron(output,eye(2));
        end
  
    
    else
        output=eye(2);
        for k=2:N
            if k==pos
                output=kron(output,sigmax());
            else
                output=kron(output,eye(2));
            end
        end
    end
end

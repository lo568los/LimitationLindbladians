function output=create_sigmay(N,pos)
    if pos==1
        output=sigmay();
        for k=2:N
            output=kron(output,eye(2));
        end
  
    
    else
        output=eye(2);
        for k=2:N
            if k==pos
                output=kron(output,sigmay());
            else
                output=kron(output,eye(2));
            end
        end
    end
end

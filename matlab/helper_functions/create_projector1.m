%create projector 1 at correct position

function output=create_projector1(N,pos)
    if pos==1
        output=projector_1();
        for k=2:N
            output=kron(output,eye(2));
        end
  
    
    else
        output=eye(2);
        for k=2:N
            if k==pos
                output=kron(output,projector_1());
            else
                output=kron(output,eye(2));
            end
        end
    end
end
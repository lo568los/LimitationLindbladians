% create projector_01 at specific position

function output=create_projector01(N,pos)
    if pos==1
        output=projector_01();
        for k=2:N-1
            output=kron(output,eye(2));
        end
  
    
    else
        output=eye(2);
        for k=2:N-1
            if k==pos
                output=kron(output,projector_01());
            else
                output=kron(output,eye(2));
            end
        end
    end
end

%% function to check if basis is orthonormal

function output = basis_is_orthonormal(F)

    num = length(F);
    
    for i  = 1:num
        for j = 1:num
            ip = trace(F{i}'*F{j});
    
            if (i == j)
                if (abs(ip-1) > 1e-9)
                    output = false;
                    return;
                end
            else
                if (abs(ip) > 1e-9)
                    output = false;
                    return;
                end
            end
        end
    end

    output = true;
    return;















end
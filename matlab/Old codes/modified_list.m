function [output] = modified_list(N,w0list,delta)
    w0list2 = zeros(N,1)+1;
    for index = length(w0list)/2:length(w0list)
        w0list2(index) = w0list(index)+delta;
    end

    output.w0list2 = w0list2;
    return
end
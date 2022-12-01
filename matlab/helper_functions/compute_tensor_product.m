% returns tensor prodduct of sets

function output = compute_tensor_product(list1, list2)
    
    output = {};
    l1 = length(list1);
    l2 = length(list2);

    for i =1:l1
        for j = 1:l2
            output{l2*i + j - l2} = kron( list1{i}, list2{j} );
        end
    end

    return;
end

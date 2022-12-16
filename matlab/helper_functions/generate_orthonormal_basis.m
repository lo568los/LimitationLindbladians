% generate_orthonormal_basis on input on number of qubits

function output = generate_orthonormal_basis(N)
    
    f = {-create_sigmaz(1,1) /sqrt(2),create_sigmaminus(1,1),create_sigmaplus(1,1),eye(2)/sqrt(2)};

    if N == 1
        output = f;
        return ;
    end

    output = compute_tensor_product(f,generate_orthonormal_basis(N-1));
end

function output=create_sigmaplus(N,pos)
    output = 0.5*(create_sigmax(N,pos)+1.0j*create_sigmay(N,pos));
end
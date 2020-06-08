function Sigma = get_sigma(Z,tau,l,nugget)
    K = length(Z);
    if K > 1
        Sigma = tau^2 .* exp(- squareform(pdist(Z,'euclidean')) ./ l);
        if nugget > 0
            Sigma = Sigma + diag(ones(1,K)*(nugget.^2));
        end
    else
        Sigma = tau^2 + nugget^2;
    end
end
    
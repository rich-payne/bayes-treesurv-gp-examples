function theint = montecarlo_int(thetree,res,nsamps)
    % approximate true marginal likelihood given hyperparameters
    tau = res.tau;
    l = res.l;
    mu = res.mu;
    Z = res.Z; 
    nugget = thetree.nugget;
    K = length(res.f);
    Sigma = tau^2 * exp(-(1/(2*l^2)) * squareform(pdist(Z,'squaredeuclidean')) ) + diag(ones(1,K)*nugget);

    % Draw Fs
    
    % nsamps = 1000;
    thechol = chol(Sigma,'lower');
    fs = mu + thechol * normrnd(0,1,K,nsamps);
    % The estimated marginal
    theint = log(mean(exp(sum(repmat(res.ns,1,nsamps) .* fs - repmat(res.a + res.b,1,nsamps) .* exp(fs),1))));
end
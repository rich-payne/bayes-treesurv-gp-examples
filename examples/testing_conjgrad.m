% Conclusion:  conjgrad is much slower in this case, but faster in
% simple cases

% The sparse matrix thing seems to be the best option at the moment

nreps = 1000;
maxK = 100;
allK = round(linspace(1,maxK,10));
TIMES = zeros(length(allK),2);
cntr = 1;
for ii=allK
    K = ii
    Sigma = get_sigma((1:K)',1,1,1e-10);
    D = diag(rand(K,1)*10);
    %rng(50)
    f = rand(K,1);
    tic
    for jj=1:nreps
        tmp = (inv(Sigma) + D) \ f;
    end
    time1 = toc;
    tmp1 = tmp;
    Ds = sqrt(D);
    ds = sqrt(diag(D))';
    tic
    for jj=1:nreps
        Sf = Sigma * f;
        SDs = Sigma .* ds;
        tmp = Sf - SDs * conjgrad(eye(K) + ds' .* SDs,SDs' * f);
    end
    time2 = toc;
    tmp2 = tmp;
    TIMES(cntr,:) = [time1, time2];    
    cntr = cntr + 1;
end

[tmp1,tmp2]

plot(TIMES)


% Try a banded matrix with sparseness
nreps = 1000;
maxK = 500;
allK = round(linspace(3,maxK,20));
TIMES2 = zeros(length(allK),3);
cntr = 1;
rho = 3;
for ii=allK
    K = ii
    %Sigma = get_sigma((1:K)',1,1,1e-10);
    Sigmainv = diag([1, (1 + rho^2) * ones(1,K-2), 1]) + ...
        diag(-rho * ones(K-1,1),1) + ...
        diag(-rho * ones(K-1,1),-1);
    B = [
        [1; (1 + rho^2) * ones(K-2,1); 1],[0;-rho * ones(K-1,1)],...
            [-rho * ones(K-1,1);0]];
    Sigmainv = spdiags(B,[0,1,-1],K,K);
    Sigmainv = 1./(1 - rho^2) .* Sigmainv;
    ds = rand(K,1)*10;
    D = diag(ds);
    SigmainvDsparse = spdiags(ds,0,Sigmainv);
    %rng(50)
    SigmainvD = full(SigmainvDsparse);
    f = rand(K,1);
    tic
    for jj=1:nreps
        tmp = (SigmainvD) \ f;
    end
    time11 = toc;
    tmp11 = tmp;
    Ds = sqrt(D);
    ds = sqrt(diag(D))';
    tic
    for jj=1:nreps
        tmp2 = SigmainvDsparse \ f;
    end
    time22 = toc;
    tmp22 = tmp;
      
    
    % A dense matrix...
    Sigma = get_sigma((1:K)',1,1,1e-10);
    tic
    for jj=1:nreps
        tmp = (inv(Sigma) + D) \ f;
    end
    time33 = toc;
    
    TIMES2(cntr,:) = [time11, time22,time33];  
    
    cntr = cntr + 1;
end

% [tmp11,tmp22]

plot(allK,TIMES2)


    
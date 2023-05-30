function [Y, X, Y_uncensored] = gen_data_cox(seed)
    rng(seed)
    n = 1000; % sample sizes
    censparms = 3.5; % censoring parameters
    p = 9;
    X = rand(n,p);
    % Add in a categorical variable
    X = [X, binornd(1,.5,n,1)];
    beta = 2; % Hazard rate over time (exponential hazard, scale parameter)
    betas = [-1 1 2 0 0 -2, -1, 1, 1.5 -1.5]';
    nu = X * betas;
    thehaz = (1/beta) * exp(nu); % Hazard function
    % Therefore each individual comes from a different exponential failure
    % distribution.
    % Failure times
    y = exprnd(1./thehaz,n,1); % the new beta is 1/hazard
    Y_uncensored = y;
    % Add some censoring
    cens = exprnd(censparms,n,1);
    ind = y <= cens;
    disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (n))));
    y(~ind) = cens(~ind);
    ds = double(ind);
    Y = [y,ds];
    X = array2table(X);
end
function [Y, X] = gen_data_true_tree(seed)
% Create custom survival function
    times = [0,1,2,3,4,5,6,7,8,9,10,11];
    surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
    % plot(times,surv)
    ht = -log(surv);
    theta1 = 1;
    surv1 = exp(-ht * theta1);

    % CART example (Figure 2)
    N = 1000; % sample sizes
    censparm = 5; % censoring parameters

    rng(seed)
    x1 = unifrnd(0,10,N,1);
    grps = {'A','B','C','D'};
    x2 = randsample(grps,N,true)';
    dat = table(x1,x2);
    y = zeros(N,1);
    I1 = ismember(x2,{'A','B'}) & x1 > 5;
    y(I1) = wblrnd(5,2,sum(I1),1);
    I2 = ismember(x2,{'A'}) & x1 <= 5;
    y(I2) = wblrnd(1,5,sum(I2),1);
    I3 = ismember(x2,{'B'}) & x1 <= 5;
    y(I3) = wblrnd(.5,.9,sum(I3),1);
    I4 = ismember(x2,{'C','D'}) & x1 <= 3;
    y(I4) = wblrnd(5,5,sum(I4),1);
    I5 = ismember(x2,{'C','D'}) & x1 > 3 & x1 <= 7;
    y(I5) = wblrnd(.5,.5,sum(I5),1);
    I6 = ismember(x2,{'C','D'}) & x1 >  7;
    y(I6) = [interp1(surv1,times,rand(sum(I6),1))];
    % Do some censoring
    cens = exprnd(censparm,N,1);
    ind = y <= cens;
    y(~ind) = cens(~ind);
    Y = [y,double(ind)];
    X = dat;
    disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (N))));
end
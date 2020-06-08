addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
% 
% a = .1;
% a0 = 1.1;
% b0 = 1.3;
% theta_shape = 1;
% theta_rate = 1;

% Try simulating from the model...
% Generate Rs on a grid
% rng(288)
% omega = .5;
% ngrid = 1000;
% xx = linspace(0,10,ngrid);
% rs = zeros(ngrid-1,1);
% c = .1;
% for ii=2:length(xx)
%     rs(ii-1) = gamrnd(c*(omega*(xx(ii) - xx(ii-1))), 1/c);
% end
% rsum = cumsum(rs);
% plot(xx(2:length(xx)),rsum)
% surv_discrete = exp(-rsum);
% plot(xx(2:length(xx)),surv_discrete)
% plot(surv_discrete,xx(2:length(xx)))
% 
% u = rand;
% 
% [newsurv,I] = sort(surv_discrete);
% newx = xx(2:length(xx));
% newx = newx(I);
% plot(surv_discrete,newx)
% 
% [usurv,ia,~] = unique(surv_discrete);
% newx = xx(2:length(xx));
% newx = newx(ia);
% plot(usurv,newx)

% Create my own survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
%plot(times,surv)
ht = -log(surv);

% Data from first function
theta1 = 1;
surv1 = exp(-ht*theta1);
theta2 = 2;
surv2 = exp(-ht*theta2);
% plot(times,surv1);
% hold on;
% plot(times,surv2);
% hold off;

% Simulate survival times from each
rng(33366)
n1 = 250;
n2 = 250;
X = table([rand(n1,1)/2; rand(n2,1)/2 + .5]);
y = [interp1(surv1,times,rand(n1,1)); interp1(surv2,times,rand(n2,1))];
% Random censoring function.
%cens = rand(n1 + n2,1) + 5; % Original Censoring...
cens = exprnd(50,n1+n2,1);
ind = y < cens;
ds = double(ind);
y(~ind) = cens(~ind);
Y = [y,ds];
[~,yind] = sort(y);
Y = Y(yind,:);
X = X(yind,:);

% Subset data for testing
Ysub1 = Y(1:250,:);
Xsub1 = X(1:250,:);
Ysub2 = Y(251:350,:);
Xsub2 = X(251:350,:);

Ystd = Y;
Ystd(:,1) = Ystd(:,1)/max(Ystd(:,1));
K = 100;
nugget = 0;

[marg_y,res] = get_marginal(Ystd,K,[],1e-9,1,.01,nugget,1);
[marg_y1,res1] = get_marginal(Ystd(X{:,1} < .5,:),K,[],1e-9,1,.1,nugget,1);
[marg_y2,res2] = get_marginal(Ystd(X{:,1} >= .5,:),K,[],1e-9,1,.1,nugget,1);
marg_y
marg_y1 + marg_y2
figure()
out = get_surv(Y,res1,1000,1,[],.05);
hold on
    plot(times,surv1,'-r')
hold off
figure()
out = get_surv(Y,res2,1000,1,[],.05);
hold on
    plot(times,surv2,'-r')
hold off





parpool(8);
tic
Tree_Surv(Y,X,'nmcmc',100,'burn',0,...
    'filepath','../output/sim2_fast_K20_rep1/','seed',1,'saveall',1,'bigK',20,...
    'swapfreq',10,'nprint',10);
time1 = toc;

tic
Tree_Surv(Y,X,'nmcmc',100,'burn',0,...
    'filepath','../output/sim2_fast_K100_rep1/','seed',1,'saveall',1,'bigK',100,...
    'swapfreq',10,'nprint',10);
time2 = toc;

time2 / time1


% Tree_Surv(Y,X,'nmcmc',10000,'burn',10000,'filepath','../output/sim2_rep1/',...
%     'seed',1001,'bigK',20,'nprint',1000,'saveall',1);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim2_rep2/',...
%     'seed',2222,'bigK',20,'nprint',1000,'saveall',1,...
%     'resume','../output/sim2_rep1/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim2_1overlprior_rep3ish/',...
%     'seed',252,'bigK',20,'nprint',1000,'saveall',1,...
%     'resume','../output/sim2_rep2/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim2_intercept_rep1/',...
%     'seed',1001,'bigK',20,'nprint',1000,'saveall',1);

% Full prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim2_fullprior_rep1/',...
%     'seed',1001,'nprint',1000,'saveall',1);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,'filepath','../output/sim2_fullprior_rep2/',...
%     'seed',191,'nprint',1000,'saveall',1,...
%     'resume','../output/sim2_fullprior_rep1/');

% % Full prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim2_fullprior_mumarg_rep1/',...
%     'seed',1001,'nprint',1000,'saveall',1);
% second iter not run!!!
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,'filepath','../output/sim2_fullprior_mumarg_rep2/',...
%     'seed',92001,'nprint',1000,'saveall',1,'resume','../output/sim2_fullprior_mumarg_rep1/');

% Full prior, mu marginalized, prior reversed
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim2_fullprior_mumarg_revprior_rep1/',...
%     'seed',1001,'nprint',1000,'saveall',1);

% Full prior, mu marginalized, prior reversed, large var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_rep1/',...
%     'seed',1001,'nprint',1000,'saveall',1);

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',1001,'saveall',1,...
%     'nprint',1000);


% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 250 points
% Tree_Surv(Ysub1,Xsub1,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/',...
%     'seed',453121,'saveall',1,'bigK',20,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 100 points
% Tree_Surv(Ysub2,Xsub2,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/',...
%     'seed',450091,'saveall',1,'bigK',20,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
%     'seed',123117,'saveall',1,'bigK',20,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/',...
%     'seed',1217,'saveall',1,'bigK',20,...
%     'nprint',1000,'resume','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/');

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
%     'seed',12311,'saveall',1,...
%     'nprint',1000,'resume','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/',...
%     'seed',112233,'saveall',1,...
%     'nprint',1000,'resume','../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/');


if 0 
    load('../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/mcmc_id1.mat')
    load('../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/mcmc_id1.mat')
    load('../output/sim2_fast_K20_rep1/mcmc_id1.mat')
    
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    % [~,I] = max(output.llike + output.prior);    
    
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
%     Ystd = Y;
%     Ystd(:,1) = Y(:,1) ./ max(Y(:,1));
%     thetree = llike_termnodes(thetree,Ystd);
    
    
    Treeplot(thetree)
    %tabulate(output.treesize)

    figure()
    get_surv_tree(thetree,Y,X,10000,1,.25,[])
    %y_max = max(Y(:,1));
    hold on
    plot(times,surv1,'r')
    hold off
    
    figure()
    get_surv_tree(thetree,Y,X,10000,1,.75,[])
    hold on
    plot(times,surv2,'r');
    hold off
    
    
    
    % Subset of the data wth 250 points
    load('../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    %tabulate(output.treesize)
    figure()
    get_surv_tree(thetree,Ysub1,Xsub1,10000,1,.25,[])
    hold on
    plot(times,surv1,'r')
    hold off
    
    % Subset of the data wth 100 points
    load('../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub2);
    Treeplot(thetree)
    figure()
    get_surv_tree(thetree,Ysub2,Xsub2,10000,1,.25,[])
    %y_max = max(Y(:,1));
    hold on
    plot(times,surv1,'r')
    hold off
    
    
    
    
    
    
    % Test out a fixed grid!
    
    
    truetree = thetree;
    truetree.Allnodes{2}.Rule(2) = {.5};
    truetree = fatten_tree(truetree,X);
    truetree = llike_termnodes(truetree,Y);
    
    
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1) ./ max(Ystd(:,1));
    ysub1 = Ystd(X{:,1} < .5,:);
    ysub2 = Ystd(X{:,1} > .5,:);
    ysub3 = Ystd(X{:,1} < .31786,:);
    ysub4 = Ystd(X{:,1} > .31786,:);
    
    % ysub(:,1) = ysub(:,1) ./ max(Y(:,1));
    K = 100;
    s = linspace(0,1.01,K); % Fixed grid
    eps = 1e-10;
    tau = 1;
    l = .01;
    nugget = 1e-10;
    EB = 1;
    [marg_y1,res1] = get_marginal(ysub1,K,s,eps,tau,l,nugget,EB);
    [marg_y2,res2] = get_marginal(ysub2,K,s,eps,tau,l,nugget,EB);
    [marg_y3,res3] = get_marginal(ysub3,K,s,eps,tau,l,nugget,EB);
    [marg_y4,res4] = get_marginal(ysub4,K,s,eps,tau,l,nugget,EB);
    marg_y1 + marg_y2
    marg_y3 + marg_y4
    
    
    get_surv(Y,res,1000,1,[],.05)
    hold on
      plot(times,surv1,'r')
    hold off
    
    
    ysub = Y(X{:,1} > .5,:);
    ysub(:,1) = ysub(:,1) ./ max(Y(:,1));
    K = 20;
    s = [];
    eps = 1e-10;
    tau = 1;
    l = 1;
    mu = 0;
    nugget = 1e-10;
    EB = 1;
    [marg_y,res2] = get_marginal(ysub,K,s,eps,tau,l,mu,nugget,EB);
    get_surv(max(Y(:,1)),res2,1000,1,[])
    hold on
      plot(times,surv2,'r')
    hold off
    


    plot(output.As)
    plot(output.Omegas)

    [GAM,grid,thetas] = plot_surv(Y,X,thetree,10,50,[],[],[],1000);
    [GAM2,grid2,thetas2] = plot_surv(Y,X,thetree,10,50,[],10,[],1000);
    [GAM3,grid3,thetas3] = plot_surv(Y,X,thetree,10,50,[],1,[],1);

    ii = 2;
    cumhaz = squeeze(GAM(ii,:,:));
    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,min(exp(-cumhaz)))
    plot(grid,max(exp(-cumhaz)))
    plot(times,surv1)
    plot(times,surv2)
    hold off


    ii = 2;
    cumhaz2 = squeeze(GAM2(ii,:,:));
    plot(grid,mean(exp(-cumhaz2)))
    hold on
    plot(grid,min(exp(-cumhaz2)))
    plot(grid,max(exp(-cumhaz2)))
    plot(times,surv1)
    plot(times,surv2)
    hold off



    plot(grid,exp(-cumhaz(2,:)))

    options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
    [~,lmax] = fminunc(@(theta) log_ptheta_t(theta,Y,thetree.a,thetree.omega,thetree.theta_shape,thetree.theta_rate),1,options);
    afunc = @(a) ptheta_t(30,Y,a,1,1,.001,lmax+1000);
    fplot(afunc,[0 25])
    fplot(@(theta) ptheta_t(theta,Y,1000,1,1,.001,lmax + 1000),[.01 5])

    afunc(10)

end            
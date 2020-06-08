addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

% Create my own survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
%plot(times,surv)
ht = -log(surv);
% Data from first function
theta1 = 1;
surv1 = exp(-ht*theta1);

% Simulate survival times from each
%rng(33366)
%n1 = 1000;
%X = table([rand(n1,1)]);
% y = [interp1(surv1,times,rand(n1,1))];
% % Random censoring function.
% cens = rand(n1,1)*12;
% ind = y < cens;
% ds = double(ind);
% y(~ind) = cens(~ind);
% Y = [y,ds];
% [~,yind] = sort(y);
% Y = Y(yind,:);
% X = X(yind,:);


% CART example (Figure 2)
rng(3328)
N = 1200;
x1 = unifrnd(0,10,N,1);
grps = {'A','B','C','D'};
x2 = randsample(grps,N,true)';
dat = table(x1,x2);
y = zeros(N,1);
I1 = ismember(x2,{'A','B'}) & x1 > 5;
y(I1) = exprnd(1,sum(I1),1);
I2 = ismember(x2,{'A','B'}) & x1 <= 5;
y(I2) = exprnd(10,sum(I2),1);
I3 = ismember(x2,{'C','D'}) & x1 <= 3;
y(I3) = exprnd(.5,sum(I3),1);
I4 = ismember(x2,{'C','D'}) & x1 > 3 & x1 <= 7;
y(I4) = exprnd(3,sum(I4),1);
I5 = ismember(x2,{'C','D'}) & x1 >  7;
%y(I5) = exprnd(10,sum(I5),1);
y(I5) = [interp1(surv1,times,rand(sum(I5),1))];
% Do some censoring
cens = exprnd(50,N,1);
ind = y <= cens;
y(~ind) = cens(~ind);
Y = [y,double(ind)];
X = dat;

% Run MCMC
parpool(8);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',10000,...
%     'filepath','../output/sim6_rep1/','seed',5250,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_rep2/','seed',52503521,'saveall',1,...
%     'nprint',1000,'resume','../output/sim6_rep1/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_1overlprior_rep3ish/','seed',12,'saveall',1,...
%     'nprint',1000,'resume','../output/sim6_rep2/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_intercept_rep1/','seed',5250,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
%     'filSubgroup identification from randomized clinical trial dataepath','../output/sim6_intercept_rep2_test/','seed',3325,'saveall',1,...
%     'nprint',100,'resume','../output/sim6_intercept_rep1_test/');
% Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
%     'filepath','../output/sim6_intercept_rep3_test/','seed',9999,'saveall',1,...
%     'nprint',100,'resume','../output/sim6_intercept_rep2_test/');
% Tree_Surv(Y,X,'nmcmc',2000,'burn',0,...
%     'filepath','../output/sim6_intercept_rep4_test/','seed',9339,'saveall',1,...
%     'nprint',100,'resume','../output/sim6_intercept_rep3_test/');
% Tree_Surv(Y,X,'nmcmc',5000,'burn',0,...
%     'filepath','../output/sim6_intercept_rep5_test/','seed',7339,'saveall',1,...
%     'nprint',100,'resume','../output/sim6_intercept_rep4_test/');
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim6_intercept_rep6_test/','seed',7332359,'saveall',1,...
%     'nprint',100,'resume','../output/sim6_intercept_rep5_test/');

% Full prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_fullprior_rep1/','seed',5250,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim6_fullprior_rep2/','seed',525990,'saveall',1,...
%     'nprint',1000,'resume','../output/sim6_fullprior_rep1/');

% Full prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_rep1/','seed',5250,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_rep2/','seed',3325,'saveall',1,...
%     'nprint',1000,'resume','../output/sim6_fullprior_mumarg_rep1/');

% Full prior, mu marginalized, reversed prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_revprior_rep1/','seed',5250,'saveall',1,...
%     'nprint',1000);

% Full prior, mu marginalized, reversed prior, large var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_revprior_bigvartau_rep1/','seed',5250,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_revprior_bigvartau_rep2/','seed',250,'saveall',1,...
%     'nprint',1000,'resume','../output/sim6_fullprior_mumarg_revprior_bigvartau_rep1/');

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',250,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim6_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
    'seed',772821,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim6_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim6_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
    'seed',5311,'saveall',1,...
    'nprint',1000,'resume','../output/sim6_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');


if 0 
    load('../output/sim6_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/mcmc_id1.mat')
    
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree)
    
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(6,{'A'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - expcdf(xx,1);
        plot(xx,yy)
        xlim([0,10])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(4,{'A'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - expcdf(xx,10);
        plot(xx,yy)
        %xlim([0,10])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(1,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,.5);
        plot(xx,yy)
        xlim([0,5])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(5,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,3);
        plot(xx,yy)
        xlim([0,15])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,10000,1,table(8,{'C'}),[])
    hold on
        plot(times,surv1)
        xlim([0,15])
    hold off
    
    % This code is designed for rep6
    newnode = thetree.Allnodes{8};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 20;
    s = [];
    eps = 1e-10;
    tau = .24;
    l = .00079;
    mu = 4.5;
    nugget= 1e-10;
    EB = 1;
    
    
    [marg_y,res] = get_marginal(Ystd(newnode.Xind,:),K,s,eps,tau,l,mu,nugget,EB)
    get_surv(Y,res,1000,1,[])
    hold on;
      xx = linspace(.01,3,100);
      plot(xx,1 - expcdf(xx,.5));
    hold off;
    
    
    get_surv_tree(thetree,Y,dat,1000,1,[],[])
    
    newnode = thetree.Allnodes{8};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 20;
    s = [];
    eps = 1e-10;
    tau = .24;
    l = .00079;
    mu = 4.5;
    nugget= 1e-10;
    EB = 1;
    lchild = thetree.Allnodes{nodeind(thetree,newnode.Lchild)};
    rchild = thetree.Allnodes{nodeind(thetree,newnode.Rchild)};
    [marg_y,res] = get_marginal(Ystd(newnode.Xind,:),K,s,eps,tau,l,mu,nugget,EB);
    [marg_y_lc,res_lc] = get_marginal(Ystd(lchild.Xind,:),K,s,eps,tau,l,mu,nugget,EB);
    [marg_y_rc,res_rc] = get_marginal(Ystd(rchild.Xind,:),K,s,eps,tau,l,mu,nugget,EB);
    montecarlo_int(thetree,res,100000)
    marg_y
    montecarlo_int(thetree,res_lc,100000)
    marg_y_lc
    montecarlo_int(thetree,res_rc,100000)
    marg_y_rc
    
    
    
    
    
    
    inds = termnodes(thetree);
    alllikes = [];
    for ii=1:length(inds)
        thelike = thetree.Allnodes{inds(ii)}.Llike;
        if ~isreal(thelike)
            ii
            thelike
        end            
        alllikes(ii) = thelike;
    end
    
    thenode = thetree.Allnodes{inds(5)};
    Ystd = Y; Ystd(:,1) = Ystd(:,1) / max(Ystd(:,1));
    K = 20;
    s = [];
    eps = 1e-10;
    tau = thenode.tau;
    l = thenode.l;
    mu = thenode.mu;
    nugget = 1e-10;
    EB = 0;
    [marg_y,res] = get_marginal(Ystd(thetree.Allnodes{inds(5)}.Xind,:),K,s,eps,tau,l,mu,nugget,EB);
    
    thetree.Allnodes{inds(5)}.Llike
    
    
    
    
        
    newtree = llike_termnodes(thetree,Ystd);
    

    get_surv_tree(thetree,Y,X,1000,1,[],[])
    get_surv_tree(thetree,Y,X,1000,1,table(3,{'D'}),[])
    hold on
        % plot(times,surv)
        xx = linspace(.01,max(Y(:,1)),100);
        plot(xx,1-expcdf(xx,.5));
    hold off

    
    [~,ids] = termnodes(thetree);
    I = ids(1);
    thetree.Allnodes{I}
    
    ystd = Y;
    ystd(:,1) = ystd(:,1) / max(Y(:,1));
    [marg_y,res] = get_marginal(ystd,20,[],1e-10,thetree.Allnodes{I}.tau,...
        thetree.Allnodes{I}.l,thetree.Allnodes{I}.mu,1e-10,0);
    
    
    
    
    
    
%     get_surv_tree(thetree,Y,X,1000,1,.7,[])
%     hold on
%         xx = linspace(.01,max(Y(:,1)),100);
%         yy = 1 - gamcdf(xx,1,1/10);
%         plot(xx,yy)
%     hold off
    
    
    
    
    
    
    
    % We have very sensitive results based on a and/or the prior for a
    tabulate(output1.treesize)
    tabulate(output2.treesize)
    tabulate(output3.treesize)
    tabulate(output4.treesize)

    [~,I1] = max(output1.llike);
    [~,I2] = max(output2.llike);
    [~,I3] = max(output3.llike);
    [~,I4] = max(output4.llike);

    tree1 = output1.Trees{I1};
    tree2 = output2.Trees{I2};
    tree3 = output3.Trees{I3};
    tree4 = output4.Trees{I4};

    [GAM1,grid1,thetas1] = plot_surv(Y,X,tree1,10,50,[],[],[],[],1000);
    [GAM2,grid2,thetas2] = plot_surv(Y,X,tree2,10,50,[],[],[],[],1000);
    [GAM3,grid3,thetas3] = plot_surv(Y,X,tree3,10,50,[],[],[],1e-12,1000);
    [GAM4,grid4,thetas4] = plot_surv(Y,X,tree4,10,50,[],[],[],1e-12,1000);

    % save('sim6-all.mat','-v7.3');

    ii = 1;
    cumhaz = squeeze(GAM4(ii,:,:));
    grid = grid4;
    qtiles = quantile(exp(-cumhaz),[.025,.975]);
    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,qtiles(1,:))
    plot(grid,qtiles(2,:));
    hold off


    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,min(exp(-cumhaz)))
    plot(grid,max(exp(-cumhaz)))
    plot(times,surv1)
    hold off




    a=1;
    omega=1;
    theta_shape = 1;
    theta_rate = .001;

    [ahat,val] = get_a(Y,omega,theta_shape,theta_rate);

    options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
    [~,lmax] = fminunc(@(theta) log_ptheta_t(theta,Y,a,omega,theta_shape,theta_rate),1,options);
    afunc = @(a) ptheta_t(theta1,Y,a,omega,theta_shape,theta_rate,lmax-1000);
    fplot(afunc,[0 500])
    fplot(@(theta) ptheta_t(theta,Y,.5,omega,theta_shape,theta_rate,lmax + 1000),[.01 5])

    % Find optimum a Value through marginalization
    [res1,C] = mt(Y,1,omega,theta_shape,theta_rate);
    log(res1) + C

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








    afunc(5)
end
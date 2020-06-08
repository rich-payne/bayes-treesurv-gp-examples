% Same as sim8 but we have changed one of the Weibull distributions

addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

% Create my own survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
%plot(times,surv)
ht = -log(surv);
% Data from first function
theta1 = 1;
surv1 = exp(-ht*theta1);

% CART example (Figure 2)
rng(424777)
%N = 5000;
N = 1200;
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
cens = exprnd(50,N,1);
ind = y <= cens;
y(~ind) = cens(~ind);
Y = [y,double(ind)];
X = dat;


Ystd = Y;
Ystd(:,1) = Ystd(:,1) ./ max(Ystd(:,1));
nugget = 0;
l = .01;
tau =.1;
EB = 1;
K = 20;
tic 
for ii=1:100
    [marg_y,res] = get_marginal(Ystd,K,[],1e-9,tau,l,nugget,EB);
end




% Subset data for testing
Ysub1 = Y(1:600,:);
Xsub1 = X(1:600,:);
Ysub2 = Y(601:700,:);
Xsub2 = X(601:700,:);

% dattosave = table(Y(:,1),Y(:,2),x1,x2);
% dattosave.Properties.VariableNames{'Var1'} = 'y';
% dattosave.Properties.VariableNames{'Var2'} = 'cens';
% writetable(dattosave,'../data/sim9.csv')

% Plot Survival functions
xx = linspace(.01,15,400);
plot(xx,1-wblcdf(xx,5,2))
hold on
    plot(xx,1-wblcdf(xx,1,5))
    plot(xx,1-wblcdf(xx,.5,.9))
    plot(xx,1-wblcdf(xx,5,5))
    plot(xx,1-wblcdf(xx,.5,.5))
    plot(times,surv)
hold off

% A & B
figure()
plot(xx,1-wblcdf(xx,5,2))
title('A and B')
hold on
    plot(xx,1-wblcdf(xx,1,5))
    plot(xx,1-wblcdf(xx,.5,.9))
hold off

% C & D 
figure()
plot(xx,1-wblcdf(xx,5,5))
title('C and D')
hold on
    plot(xx,1-wblcdf(xx,.5,.5))
    plot(times,surv)
hold off

parpool(8);
tic
Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
    'filepath','../output/sim9_fast_K20_rep1/','seed',1133,'saveall',1,...
    'nprint',100,'swapfreq',10,'bigK',20);
time1 = toc;

Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
    'filepath','../output/sim9_fast_K20_rep2/','seed',339,'saveall',1,...
    'nprint',100,'swapfreq',10,'bigK',20,...
    'resume','../output/sim9_fast_K20_rep1/');

Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim9_fast_K20_rep3/','seed',811923,'saveall',1,...
    'nprint',100,'swapfreq',10,'bigK',20,...
    'resume','../output/sim9_fast_K20_rep2/');





Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim9_fast_K100_rep1/','seed',1133,'saveall',1,...
    'nprint',100,'swapfreq',10,'bigK',100);

% Full prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_rep1/','seed',18211,'saveall',1,...
%     'nprint',100);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_rep2/','seed',11,'saveall',1,...
%     'nprint',100,'resume','../output/sim9_fullprior_mumarg_rep1/');

% Full prior, mu marginalized, reversed prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_rep1/','seed',18211,'saveall',1,...
%     'nprint',100);

% Full prior, mu marginalized, reversed prior, large var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_rep1/','seed',18211,'saveall',1,...
%     'nprint',100);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_rep2/','seed',18211,'saveall',1,...
%     'nprint',100,'resume','../output/sim9_fullprior_mumarg_revprior_bigvartau_rep1/');

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',18211,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 600 points
Tree_Surv(Ysub1,Xsub1,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/',...
    'seed',53121,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 100 points
Tree_Surv(Ysub2,Xsub2,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/',...
    'seed',50091,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
%     'seed',5551,'saveall',1,'bigK',20,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/',...
%     'seed',51,'saveall',1,'bigK',20,...
%     'nprint',1000,'resume','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/');
Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep3/',...
    'seed',1515,'saveall',1,'bigK',20,...
    'nprint',1000,'resume','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/');

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
%     'seed',2351,'saveall',1,...
%     'nprint',1000,'resume','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/',...
%     'seed',1909,'saveall',1,...
%     'nprint',1000,'resume','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/');
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep4/',...
    'seed',199,'saveall',1,...
    'nprint',1000,'resume','../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/');

if 0
    % You ought to check to see if there are any gaps or clusters of data 
    %   which may be fixed by adjusting the number of bins, etc.
    load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/mcmc_id1.mat');
    load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/mcmc_id1.mat');
    load('../output/sim9_fast_K20_rep3/mcmc_id1.mat')
    
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree)
    
    get_surv_tree(thetree,Y,X,1000,1,[],[])
    
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(6,{'A'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - wblcdf(xx,5,2);
        plot(xx,yy)
        xlim([0,15])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(4,{'A'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - wblcdf(xx,1,5);
        plot(xx,yy)
        xlim([0,3])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(4,{'B'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - wblcdf(xx,.5,.9);
        plot(xx,yy)
        xlim([0,5])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(1,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - wblcdf(xx,5,5);
        plot(xx,yy)
        xlim([0,5])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(5,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - wblcdf(xx,.5,.5);
        plot(xx,yy)
        xlim([0,15])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(8,{'C'}),[])
    hold on
        plot(times,surv1)
        xlim([0,15])
    hold off
    
    
    % Subset of the data wth 500 points
    load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    
    % Subset of the data wth 100 points
    load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub2);
    Treeplot(thetree)
    
    
    % 9 i s 7.7085
    % 10 is "c"
    newnode = thetree.Allnodes{6};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 100;
    %s = linspace(0,1.01,K+1);
    s = [];
    eps = 1e-10;
    tau = 1;
    l = .01;
    %mu = 4.5;
    nugget= 1e-10;
    EB = 1;
    lchild = thetree.Allnodes{nodeind(thetree,newnode.Lchild)};
    rchild = thetree.Allnodes{nodeind(thetree,newnode.Rchild)};
    [marg_y,res] = get_marginal(Ystd(newnode.Xind,:),K,s,eps,tau,l,nugget,EB);
    [marg_y_lc,res_lc] = get_marginal(Ystd(lchild.Xind,:),K,s,eps,tau,l,nugget,EB);
    [marg_y_rc,res_rc] = get_marginal(Ystd(rchild.Xind,:),K,s,eps,tau,l,nugget,EB);
    marg_y
    marg_y_lc + marg_y_rc
    
    
    % Find the true marginal
    I1 = ismember(x2,{'A','B'}) & x1 > 5;
    [marg_y1,res1] = get_marginal(Ystd(I1,:),K,s,eps,tau,l,nugget,EB);
    I2 = ismember(x2,{'A'}) & x1 <= 5;
    [marg_y2,res2] = get_marginal(Ystd(I2,:),K,s,eps,tau,l,nugget,EB);
    I3 = ismember(x2,{'B'}) & x1 <= 5;
    [marg_y3,res3] = get_marginal(Ystd(I3,:),K,s,eps,tau,l,nugget,EB);
    I4 = ismember(x2,{'C','D'}) & x1 <= 3;
    [marg_y4,res4] = get_marginal(Ystd(I4,:),K,s,eps,tau,l,nugget,EB);
    I5 = ismember(x2,{'C','D'}) & x1 > 3 & x1 <= 7;
    [marg_y5,res5] = get_marginal(Ystd(I5,:),K,s,eps,tau,l,nugget,EB);
    I6 = ismember(x2,{'C','D'}) & x1 >  7;
    [marg_y6,res6] = get_marginal(Ystd(I6,:),K,s,eps,tau,l,nugget,EB);
    
    thetree.Lliketree
    marg_y1 + marg_y2 + marg_y3 + marg_y4 + marg_y5 + marg_y6
   
    
    
    
   
    
    newnode
    lchild
    rchild
    
    montecarlo_int(thetree,res,100000)
    marg_y - get_lprior(newnode.tau,newnode.l,newnode.mu)
    montecarlo_int(thetree,res_lc,100000)
    marg_y_lc - get_lprior(lchild.tau,lchild.l,lchild.mu)
    montecarlo_int(thetree,res_rc,100000)
    marg_y_rc - get_lprior(rchild.tau,rchild.l,rchild.mu)
    
    %get_surv(Y_orig,res,ndraw,graph,ystar)
    get_surv(Y,res,1000,1,[],.05)
    hold on 
        plot(times,surv)
    hold off 
    figure()
    get_surv(Y,res_lc,1000,1,[],.05)
    hold on 
        plot(times,surv)
    hold off 
    figure()
    get_surv(Y,res_rc,1000,1,[],.05)
    hold on 
        plot(times,surv)
    hold off 
    
    plot(res.f .* res.ns)
    hold on
        plot(res_lc.f .* res_lc.ns)
        plot(res_rc.f .* res_rc.ns)
    hold off
    
    
    figure()
    plot(-exp(res.f) .* (res.a + res.b))
    hold on
        plot(-exp(res_lc.f) .* (res_lc.a + res_lc.b))
        plot(-exp(res_rc.f) .* (res_rc.a + res_rc.b))
    hold off
    
    figure()
    plot(res.f .* res.ns -exp(res.f) .* (res.a + res.b))
    hold on
        plot(res_lc.f.*res_lc.ns -exp(res_lc.f) .* (res_lc.a + res_lc.b))
        plot(res_rc.f.*res_rc.ns-exp(res_rc.f) .* (res_rc.a + res_rc.b))
    hold off
    
    llike = sum(res.f .* res.ns -exp(res.f) .* (res.a + res.b));
    llike_lc = sum(res_lc.f.*res_lc.ns -exp(res_lc.f) .* (res_lc.a + res_lc.b));
    llike_rc = sum(res_rc.f.*res_rc.ns-exp(res_rc.f) .* (res_rc.a + res_rc.b));
    
    Sigma = get_sigma(res.Z,res.tau,res.l,sqrt(nugget));
    Sigma_lc = get_sigma(res_lc.Z,res_lc.tau,res_lc.l,sqrt(nugget));
    Sigma_rc = get_sigma(res_rc.Z,res_rc.tau,res_rc.l,sqrt(nugget));
    eprior = -.5 * (res.f - res.mu)' * (Sigma \ (res.f - res.mu));
    eprior_lc = -.5 * (res_lc.f - res_lc.mu)' * (Sigma_lc \ (res_lc.f - res_lc.mu));
    eprior_rc = -.5 * (res_rc.f - res_rc.mu)' * (Sigma_rc \ (res_rc.f - res_rc.mu));
    
    thedet = -.5*ldet(res.Omegainv) - .5*ldet(Sigma);
    thedet_lc = -.5*ldet(res_lc.Omegainv) - .5*ldet(Sigma_lc);
    thedet_rc = -.5*ldet(res_rc.Omegainv) - .5*ldet(Sigma_rc);
    
    plot([llike,eprior,thedet,marg_y],'b')
    hold on
        plot([llike_lc,eprior_lc,thedet_lc,marg_y_lc],'r')
        plot([llike_rc,eprior_rc,thedet_rc,marg_y_rc],'g')
        plot([llike_lc + llike_rc,...
              eprior_lc + eprior_rc,...
              thedet_lc + thedet_rc,marg_y_lc + marg_y_rc],'m')
    hold off
    
    marg_y - get_lprior(res.tau,res.l,res.mu)
    llike + eprior + thedet
    
end


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
y(I4) = wblrnd(3,3,sum(I4),1);
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

% Subset data for testing
Ysub1 = Y(1:600,:);
Xsub1 = X(1:600,:);
Ysub2 = Y(601:700,:);
Xsub2 = X(601:700,:);


% Plot Survival functions
xx = linspace(.01,15,400);
plot(xx,1-wblcdf(xx,5,2))
hold on
    plot(xx,1-wblcdf(xx,1,5))
    plot(xx,1-wblcdf(xx,.5,.9))
    plot(xx,1-wblcdf(xx,3,3))
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
plot(xx,1-wblcdf(xx,3,3))
title('C and D')
hold on
    plot(xx,1-wblcdf(xx,.5,.5))
    plot(times,surv)
hold off

parpool(8);

% Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
%     'filepath','../output/sim8_fullprior_rep1/','seed',18211,'saveall',1,...
%     'nprint',100);
% Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
%     'filepath','../output/sim8_fullprior_rep2/','seed',118,'saveall',1,...
%     'nprint',100,'resume','../output/sim8_fullprior_rep1/');
% Tree_Surv(Y,X,'nmcmc',5000,'burn',0,...
%     'filepath','../output/sim8_fullprior_rep3/','seed',53211,'saveall',1,...
%     'nprint',100,'resume','../output/sim8_fullprior_rep2/');


% Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
%     'filepath','../output/sim8_fullprior_bigN_rep1/','seed',55309,'saveall',1,...
%     'nprint',100);
% Tree_Surv(Y,X,'nmcmc',2000,'burn',0,...
%     'filepath','../output/sim8_fullprior_bigN_rep2/','seed',519,'saveall',1,...
%     'nprint',100,'resume','../output/sim8_fullprior_bigN_rep1/');
% Tree_Surv(Y,X,'nmcmc',5000,'burn',0,...
%     'filepath','../output/sim8_fullprior_bigN_rep3/','seed',915,'saveall',1,...
%     'nprint',100,'resume','../output/sim8_fullprior_bigN_rep2/');


% Tree_Surv(Y,X,'nmcmc',5000,'burn',0,...
%     'filepath','../output/sim8_fullprior_K40_rep1/','seed',915,'saveall',1,...
%     'nprint',100,'bigK',40);

% Tree_Surv(Y,X,'nmcmc',5000,'burn',0,...
%     'filepath','../output/sim8_fullprior_K40_rep2/','seed',915,'saveall',1,...
%     'nprint',100,'bigK',40,'resume','../output/sim8_fullprior_K40_rep1/');

% Full prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_rep1/','seed',18211,'saveall',1,...
%     'nprint',100);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_rep2/','seed',3532,'saveall',1,...
%     'nprint',100,'resume','../output/sim8_fullprior_mumarg_rep1/');

% Full prior, mu marginalized, reversed prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_rep1/','seed',18211,'saveall',1,...
%     'nprint',100);

% Full prior, mu marginalized, reversed prior, large var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_rep1/','seed',18211,'saveall',1,...
%     'nprint',100);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_rep2/','seed',11,'saveall',1,...
%     'nprint',100,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_rep1/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_rep3/','seed',11,'saveall',1,...
%     'nprint',100,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_rep2/');

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_rep2/','seed',12,'saveall',1,...
%     'nprint',1000,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_rep1/');

% Reversed prior, larger var for tau, l has t-prior, fixed hessian and
%   gradient
% Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/','seed',12,'saveall',1,...
%     'nprint',100);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',1233,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 600 points
Tree_Surv(Ysub1,Xsub1,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/',...
    'seed',253121,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 100 points
Tree_Surv(Ysub2,Xsub2,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/',...
    'seed',250091,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
%     'seed',887321,'saveall',1,'bigK',20,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/',...
%     'seed',81,'saveall',1,'bigK',20,...
%     'nprint',1000,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/');
Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep3/',...
    'seed',901,'saveall',1,'bigK',20,...
    'nprint',1000,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/');

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
%     'seed',9901,'saveall',1,...
%     'nprint',1000,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/',...
%     'seed',3511,'saveall',1,...
%     'nprint',1000,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/');
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep4/',...
    'seed',1135,'saveall',1,...
    'nprint',1000,'resume','../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/');

if 0
    % You ought to check to see if there are any gaps or clusters of data 
    %   which may be fixed by adjusting the number of bins, etc.
    load('../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/mcmc_id1.mat');
    load('../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/mcmc_id1.mat');
    
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
        yy = 1 - wblcdf(xx,3,3);
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
    load('../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    
    % Subset of the data wth 100 points
    load('../output/sim8_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub2);
    Treeplot(thetree)
    
   
    % 9 i s 7.7085
    % 10 is "c"
    newnode = thetree.Allnodes{1};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 100;
    s = [];
    eps = 1e-10;
    tau = 1;
    l = .01;
    mu = 4.5;
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
    
    
    
    
    plot(res.f)
    
    
    newnode
    lchild
    rchild
    
    length(lchild.Xind)
    sum(res_lc.ns)
    length(rchild.Xind)
    sum(res_rc.ns)
    length(newnode.Xind)
    sum(res.ns)
    
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
    
    plot(res.f)
    hold on;
        plot(res_lc.f)
        plot(res_rc.f)
    hold off;
    
    plot(res.ns)
    hold on;
        plot(res_lc.ns)
        plot(res_rc.ns)
    hold off;
    
    plot(res.f .* res.ns)
    hold on
        plot(res_lc.f .* res_lc.ns)
        plot(res_rc.f .* res_rc.ns)
        plot(res_rc.f .* res_rc.ns + ...
             res_lc.f .* res_lc.ns)
    hold off
    sum(res.f .* res.ns)
    sum(res_rc.f .* res_rc.ns + res_lc.f .* res_lc.ns)
    plot(res.f .* res.ns - (res_rc.f .* res_rc.ns + res_lc.f .* res_lc.ns))
    
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
    
    
    
    % FULL MARGINALIZATION
    % 9 i s 7.7085
    % 10 is "c"
    newnode = thetree.Allnodes{1};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 20;
    s = [];
    eps = 1e-10;
    tau = .01;
    l = .01;
    mu = 4.5;
    nugget= 1e-10;
    EB = 1;
    lchild = thetree.Allnodes{nodeind(thetree,newnode.Lchild)};
    rchild = thetree.Allnodes{nodeind(thetree,newnode.Rchild)};
    [marg_y,res] = get_marginal(Ystd(newnode.Xind,:),K,s,eps,tau,l,nugget,EB);
    [marg_y_lc,res_lc] = get_marginal(Ystd(lchild.Xind,:),K,s,eps,tau,l,nugget,EB);
    [marg_y_rc,res_rc] = get_marginal(Ystd(rchild.Xind,:),K,s,eps,tau,l,nugget,EB);
    marg_y
    marg_y_lc + marg_y_rc
    
      
    newtree = llike_termnodes(thetree,Ystd);
    newnode = newtree.Allnodes{1};
    newtree = llike(newtree,newnode.Id,Ystd);
    newnode = newtree.Allnodes{1};
    newlchild = newtree.Allnodes{nodeind(newtree,newnode.Lchild)};
    newrchild = newtree.Allnodes{nodeind(newtree,newnode.Rchild)};
    
    %ypart = Ystd(newnode.Xind,:);
    %[fullmarg,ZZ] = get_fullmarginal(Ystd,thenode,thetree,delta_z,delta_pi)
    delta_z = 1;
    delta_pi = 2.5;
    [fullmarg,ZZ] = get_fullmarginal(Ystd,newnode,newtree,delta_z,delta_pi)
    [fullmarg_lc,ZZ] = get_fullmarginal(Ystd,newnode,newtree,delta_z,delta_pi)
    [fullmarg_rc,ZZ] = get_fullmarginal(Ystd,newnode,newtree,delta_z,delta_pi)
    
    scatter3(ZZ(:,1),ZZ(:,2),ZZ(:,3))
    
    
    loglikefunc(newnode,newtree,Ystd)
    
    % Explore the posterior 
    dl = 1;
    dtau = 5;
    nd = 20;
    %lgrid = linspace(thenode.l - nd*dl,thenode.l + nd*dl,nd + 1);
    %taugrid = linspace(thenode.tau - nd*dtau,thenode.tau + nd*dtau,nd+1)
    lgrid = linspace(.001,.05,50);
    taugrid = linspace(.01,5,50);
    negind = lgrid <= 0;
    lgrid(negind) = [];
    negind = taugrid <= 0;
    taugrid(negind) = [];
    
    
    thenode = newnode;
    Z = zeros(length(taugrid),length(lgrid));
    for ii=1:length(taugrid)
        for jj=1:length(lgrid)
            Z(ii,jj) = get_marginal(Ystd(thenode.Xind,:),K,s,eps,taugrid(ii),lgrid(jj),nugget,0);
        end
    end
    
    thenode = newlchild;
    Z_lc = zeros(length(taugrid),length(lgrid));
    for ii=1:length(taugrid)
        for jj=1:length(lgrid)
            Z_lc(ii,jj) = get_marginal(Ystd(thenode.Xind,:),K,s,eps,taugrid(ii),lgrid(jj),nugget,0);
        end
    end
    
    thenode = newrchild;
    Z_rc = zeros(length(taugrid),length(lgrid));
    for ii=1:length(taugrid)
        for jj=1:length(lgrid)
            Z_rc(ii,jj) = get_marginal(Ystd(thenode.Xind,:),K,s,eps,taugrid(ii),lgrid(jj),nugget,0);
        end
    end
    
    [myX,myY] = meshgrid(taugrid,lgrid);
    subplot(3,1,1)
    surf(myX,myY,exp(Z'))
    subplot(3,1,2)
    surf(myX,myY,exp(Z_lc'))
    subplot(3,1,3)
    surf(myX,myY,exp(Z_rc'))
    
    max(max(exp(Z)))
    max(max(exp(Z_lc))) * max(max(exp(Z_rc)))
    
    % Approximate marginal thing
    v = (taugrid(2) - taugrid(1)) * (lgrid(2) - lgrid(1));
    v * sum(sum(exp(Z)))
    (v * sum(sum(exp(Z_lc)))) * (v* sum(sum(exp(Z_rc))))
    
    
end


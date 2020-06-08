addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
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
y(I5) = exprnd(10,sum(I5),1);
% Do some censoring
cens = exprnd(50,N,1);
ind = y <= cens;
y(~ind) = cens(~ind);
Y = [y,double(ind)];
[~,I] = sort(Y(:,1));
Y = Y(I,:);
dat = dat(I,:);
X = dat;

% Subset data for testing
Ysub1 = Y(1:600,:);
Xsub1 = X(1:600,:);
Ysub2 = Y(601:700,:);
Xsub2 = X(601:700,:);

parpool(8);
Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
    'filepath','../output/sim3_fast_K20_rep1/','seed',3003,'nprint',100,...
    'saveall',1,'bigK',20,'swapfreq',10);

Tree_Surv(Y,X,'nmcmc',1000,'burn',0,...
    'filepath','../output/sim3_fast_K20_rep2/','seed',4004,'nprint',100,...
    'saveall',1,'bigK',20,'swapfreq',10,...
    'resume','../output/sim3_fast_K20_rep1/');


% Tree_Surv(Y,dat,'nmcmc',10000,'burn',10000,...
%     'filepath','../output/sim3_rep1/','seed',3003,'nprint',1000,...
%     'saveall',1);
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_rep2/','seed',35221,'nprint',1000,...
%     'saveall',1,'resume','../output/sim3_rep1/');
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_1overlprior_rep3ish/','seed',53918,'nprint',1000,...
%     'saveall',1,'resume','../output/sim3_rep2/');
% Tree_Surv(Y,dat,'nmcmc',1000,'burn',0,...
%     'filepath','../output/sim3_intercept_rep1/','seed',3003,'nprint',1000,...
%     'saveall',1);
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_intercept_rep2/','seed',4004,'nprint',1000,...
%     'saveall',1,'resume','../output/sim3_intercept_rep1/');

% Full Prior
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_fullprior_rep1/','seed',3003,'nprint',1000,...
%     'saveall',1);
% Tree_Surv(Y,dat,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim3_fullprior_rep2/','seed',5903,'nprint',1000,...
%     'saveall',1,'resume','../output/sim3_fullprior_rep1/');


% Full prior, mu marganalized
% % Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
% %     'filepath','../output/sim3_fullprior_mumarg_rep1/','seed',3003,'nprint',1000,...
% %     'saveall',1);
% Tree_Surv(Y,dat,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_rep2/','seed',382291,'nprint',1000,...
%     'saveall',1,'resume','../output/sim3_fullprior_mumarg_rep1/');

% Full prior, mu marginalized, reversed prior
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_revprior_rep1/','seed',3003,'nprint',1000,...
%     'saveall',1);

% Full prior, mu marginalized, reversed prior, large var for tau
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_rep1/','seed',3003,'nprint',1000,...
%     'saveall',1);

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',3003,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 600 points
Tree_Surv(Ysub1,Xsub1,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/',...
    'seed',853121,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 100 points
Tree_Surv(Ysub2,Xsub2,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/',...
    'seed',850091,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
% Tree_Surv(Y,dat,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
%     'seed',311,'saveall',1,'bigK',20,...
%     'nprint',1000);
Tree_Surv(Y,dat,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep2/',...
    'seed',1188,'saveall',1,'bigK',20,...
    'nprint',1000,'resume','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/');

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,dat,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,dat,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
%     'seed',8632,'saveall',1,...
%     'nprint',1000,'resume','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');
Tree_Surv(Y,dat,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep3/',...
    'seed',2991,'saveall',1,...
    'nprint',1000,'resume','../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/');

if 0
    load('../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/mcmc_id1.mat')
    load('../output/sim3_fast_K20_rep1/mcmc_id1.mat')
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike);
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,dat);
    Treeplot(thetree)
    %plot(output.treesize)
    %tabulate(output.treesize)
    get_surv_tree(thetree,Y,dat,1000,1,[],[])
    
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
        % Kaplan Meijer Curve
%         thenode = thetree.Allnodes{9};
%         [f,x,flo,fup] = ecdf(Y(thenode.Xind,1),'censoring',Y(thenode.Xind,2) == 0,'bounds','on');
%         plot(x,1-f,x,1-flo,'--',x,1-fup,'--')
    hold off
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(8,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,10);
        plot(xx,yy)
    hold off
    
    
    % Subset of the data wth 600 points
    load('../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    
    
    
    %ind = ismember(dat{:,2},{'C','D'}) & dat{:,1} > 3.0274 & dat{:,1} <= 6.9197;
    % sum(Y(ind,:))
    thenode = thetree.Allnodes{9};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 20;
    eps=1e-10;
    tau=1e-10;
    l = thenode.l;
    tau = thenode.tau;
    nugget = 1e-10;
    EB = 0;
    %(Y,K,s,eps,tau,l,nugget,EB)
    [marg_y,res] = get_marginal(Ystd(thenode.Xind,:),K,[],eps,tau,l,nugget,EB);
    
    
    figure()
    get_surv_tree(thetree,Y,dat,1000,1,table(8,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,10);
        plot(xx,yy)
        % xlim([0,15])
    hold off
    
    

    get_surv_tree(thetree,Y,dat,1000,1,[],[])

    tind = output.treesize == 2;
    [~,I] = max(output.llike(tind));
    thetree = output.Trees(tind);
    thetree = thetree{I};
    Treeplot(thetree);
    
    
    plot(output.As)
    plot(output.Omegas)

    [GAM,grid,thetas] = plot_surv(Y,dat,thetree,10,50,[],[],[],100);

    ii = 5;
    cumhaz = squeeze(GAM(ii,:,:));
    qtiles = quantile(exp(-cumhaz),[.025,.975]);
    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,qtiles(1,:))
    plot(grid,qtiles(2,:));
    fplot(@(x) exp(-x/10),[0,60])
    hold off
end

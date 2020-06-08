% See how we do with spurious Covariates
addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

rng(33366)
n1 = 50;
n2 = 50;
X = table(rand(n1+n2,1),rand(n1+n2,1),normrnd(0,1,n1+n2,1));
ind = X{:,1} < .25;
y = zeros(n1+n1,1);
y(ind) = exprnd(5,sum(ind),1);
y(~ind) = exprnd(1,sum(~ind),1);
% Random censoring function.
cens = exprnd(10,n1+n2,1);
ind = y < cens;
ds = double(ind);
y(~ind) = cens(~ind);
Y = [y,ds];
[~,yind] = sort(y);
Y = Y(yind,:);
X = X(yind,:);


parpool(8);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',10000,'filepath','../output/sim4_rep1/','seed',4004,...
%     'nprint',1000,'saveall',1);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim4_rep2/','seed',35599,...
%     'nprint',1000,'saveall',1,'resume','../output/sim4_rep1/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim4_1overlprior_rep3ish/','seed',99,...
%     'nprint',1000,'saveall',1,'resume','../output/sim4_rep2/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim4_intercept_rep1/','seed',4004,...
%     'nprint',1000,'saveall',1);

% Full prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim4_fullprior_rep1/','seed',4004,...
%     'nprint',1000,'saveall',1);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,'filepath','../output/sim4_fullprior_rep2/','seed',4111204,...
%     'nprint',1000,'saveall',1,'resume','../output/sim4_fullprior_rep1/');

% Full prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim4_fullprior_mumarg_rep1/','seed',4004,...
%     'nprint',1000,'saveall',1);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,'filepath','../output/sim4_fullprior_mumarg_rep2/','seed',82219,...
%     'nprint',1000,'saveall',1,'resume','../output/sim4_fullprior_mumarg_rep1/');

% Full prior, mu marginalized, reversed prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim4_fullprior_mumarg_revprior_rep1/','seed',4004,...
%     'nprint',1000,'saveall',1);

% Full prior, mu marginalized, reversed prior, large var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,'filepath','../output/sim4_fullprior_mumarg_revprior_bigvartau_rep1/','seed',4004,...
%     'nprint',1000,'saveall',1);

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim4_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim4_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',4004,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim4_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
    'seed',9721,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim4_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
    'seed',12,'saveall',1,...
    'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim4_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
%     'seed',12991,'saveall',1,...
%     'nprint',1000,'resume','../output/sim4_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');


if 0 
    load('../output/sim4_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/mcmc_id1.mat')

    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree)
    %plot(output.treesize)
    %tabulate(output.treesize)
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,[.1 0 0 ],[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - expcdf(xx,5);
        plot(xx,yy)
    hold off

    figure()
    get_surv_tree(thetree,Y,X,1000,1,[.7 0 0],[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - expcdf(xx,1);
        plot(xx,yy)
    hold off

    
    
    

    tind = output.treesize == 2;
    [~,I] = max(output.llike(tind));
    thetree = output.Trees(tind);
    thetree = thetree{I};
    Treeplot(thetree);

    plot(output.As)
    plot(output.Omegas)

    [GAM,grid,thetas] = plot_surv(Y,X,thetree,10,50,[],[],[],1000);

    ii = 2;
    cumhaz = squeeze(GAM(ii,:,:));
    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,min(exp(-cumhaz)))
    plot(grid,max(exp(-cumhaz)))
    plot(times,surv1)
    plot(times,surv2)
    hold off
end






            





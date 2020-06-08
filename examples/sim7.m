addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

% Generate some data from exponential proportional hazards model
rng(32551)
n = 1000;
p = 5;
X = rand(n,p);
% Add in a categorical variable
X = [X, binornd(1,.5,n,1)];
beta = 2; % Hazard rate over time (exponential hazard, scale parameter)
betas = [-1 1 2 0 0 -2]';
nu = X * betas;
thehaz = (1/beta) * exp(nu); % Hazard function
% Therefore each individual comes from a different exponential failure
% distribution.
% Failure times
y = exprnd(1./thehaz,n,1); % the new beta is 1/hazard
% Add some censoring
cens = exprnd(50,n,1);
ind = y <= cens;
y(~ind) = cens(~ind);
ds = double(ind);

Y = [y,ds];
X = array2table(X);

parpool(8);
% leafmin = 3;
% Tree_Surv(Y,X,'nmcmc',10000,'burn',10000,...
%     'filepath','../output/sim7_rep1/','seed',1,'saveall',1,...
%     'leafmin',leafmin,'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',10000,...
%     'filepath','../output/sim7_rep1/','seed',1,'saveall',1,...
%     'leafmin',25,'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_rep2/','seed',10023,'saveall',1,...
%     'leafmin',25,'nprint',1000,'resume','../output/sim7_rep1/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_1overlprior_rep3ish/','seed',301,'saveall',1,...
%     'leafmin',25,'nprint',1000,'resume','../output/sim7_rep2/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_intercept_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);

% Full Prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_fullprior_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim7_fullprior_rep2/','seed',1871,'saveall',1,...
%     'nprint',1000,'resume','../output/sim7_fullprior_rep1/');

% Full Prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_rep2/','seed',266,'saveall',1,...
%     'nprint',1000,'resume','../output/sim7_fullprior_mumarg_rep1/');

% Full Prior, mu marginalized, reversed prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_revprior_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);

% Full Prior, mu marginalized, reversed prior, large var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',1,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
%     'seed',351942,'saveall',1,'bigK',20,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);

% Subset data for testing
Ysub1 = Y(1:500,:);
Xsub1 = X(1:500,:);
Ysub2 = Y(501:600,:);
Xsub2 = X(501:600,:);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 500 points
Tree_Surv(Ysub1,Xsub1,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/',...
    'seed',353121,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20, subset of 100 points
Tree_Surv(Ysub2,Xsub2,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/',...
    'seed',350091,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
    'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
    'seed',8321,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
    'seed',1289,'saveall',1,...
    'nprint',1000);

if 0
    load('../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/mcmc_id1.mat')
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    Treeplot(thetree)
    get_surv_tree(thetree,Y,X,1000,1,[],[])
    
    x0 = [.5 .1 .8 .5 .5];
    res = get_surv_tree(thetree,Y,X,1000,1,x0,[]);
    x0hazard = 1/beta * exp(x0 * betas);
    orig_grid = res.ystar*max(Y(:,1));
    hold on;
    plot(orig_grid,1-expcdf(orig_grid,1/x0hazard))
    hold off;
    
    
    
    % Subset of the data wth 500 points
    load('../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    
    % Subset of the data wth 100 points
    load('../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub2);
    Treeplot(thetree)
    
    
    
    
    
    
    
    % Quick check on something...
    Xtmp = X{(X{:,2} <= .22947) & (X{:,3} > .76316),:};
    allinvhazards = 1./(1/beta * exp(Xtmp *  betas));
    allsurvs = zeros(size(Xtmp,1),length(orig_grid));
    for ii=1:size(Xtmp,1)
        allsurvs(ii,:) = 1 - expcdf(orig_grid,allinvhazards(ii));
    end
    plot(allsurvs')
    avgsurv = mean(allsurvs);
    hold on
    plot(res.ystar,avgsurv)
    hold off
    
    

    % Ensemble method...
    llikes = output.llike;
    [ul,ui] = unique(llikes);
    %[ul_sorted,uis] = sort(ul);
    a = max(ul);
    tree_probs = exp(ul - (a + log(sum(exp(ul - a))))); % in smallest to largest order
    nmc = 10000; % number of monte carlo runs;
    x0 = array2table([.1 .1 .9 .9 .5]); % Value of the covariate
    % sample the trees
    treeind = randsample(length(tree_probs),nmc,true,tree_probs);
    %treeind_u = unique(treeind);
    tograph = 0;
    ystar = linspace(.01,.99,102);
    pdraws = zeros(nmc,length(ystar));
    theind = 1;
    for ii=1:length(treeind)
        ndraw = sum(treeind == ii);
        if ndraw > 0
            tmp = get_surv_tree(output.Trees{ui(ii)},Y,X,ndraw,tograph,x0,ystar);
            pdraws(theind:(theind + ndraw -1),:) = tmp.surv;
            theind = theind + ndraw;
        end
    end
    % Now plot posterior draws...
    pmean = mean(pdraws);
    qtiles = quantile(pdraws,[.025,.975],1);
    orig_grid = ystar*max(Y(:,1));
    plot(orig_grid,pmean,'-k',orig_grid,qtiles(1,:),'--k',orig_grid,qtiles(2,:),'--k')
    % True survival function
    x0hazard = 1/beta * exp(x0{1,:} * betas);
    hold on;
    plot(orig_grid,1-expcdf(orig_grid,1/x0hazard))
    hold off;
end









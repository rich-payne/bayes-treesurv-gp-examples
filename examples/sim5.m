% Compares two weibull survival functions
addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
% Simulate data
rng(38230424);
n = 100;
x = rand(n,1);
ind = x < .5;
y = zeros(n,1);
y(ind) = wblrnd(5,2,sum(ind),1);
y(~ind) = wblrnd(1,5,sum(~ind),1);
% Do some censoring
cens = gamrnd(1,5,n,1);
%cens = 1000 * ones(n,1);
ind = y < cens;
ds = double(ind); % 1 if not censored, 0 if censored;
y(~ind) = cens(~ind);

X = table(x);
[~,I] = sort(y);
y = y(I);
ds = ds(I);
X = X(I,:);

Y = [y, ds];

parpool(8);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',10000,...
%     'filepath','../output/sim5_rep1/','seed',550,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_rep2/','seed',550551,'saveall',1,...
%     'nprint',1000,'resume','../output/sim5_rep1/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_1overlprior_rep3ish/','seed',3332,'saveall',1,...
%     'nprint',1000,'resume','../output/sim5_rep2/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_intercept_rep1/','seed',550,'saveall',1,...
%     'nprint',1000);

% Full Prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_fullprior_rep1/','seed',550,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim5_fullprior_rep2/','seed',559930,'saveall',1,...
%     'nprint',1000,'resume','../output/sim5_fullprior_rep1/');

% Full Prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_fullprior_mumarg_rep1/','seed',550,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim5_fullprior_mumarg_rep2/','seed',909,'saveall',1,...
%     'nprint',1000,'resume','../output/sim5_fullprior_mumarg_rep1/');

% Full Prior, mu marginalized, reversed prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_fullprior_mumarg_revprior_rep1/','seed',550,'saveall',1,...
%     'nprint',1000);

% Full Prior, mu marginalized, reversed prior, large var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_fullprior_mumarg_revprior_bigvartau_rep1/','seed',550,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim5_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',909,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim5_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
    'seed',2317,'saveall',1,'bigK',20,...
    'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim5_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);
Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
    'filepath','../output/sim5_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
    'seed',9013,'saveall',1,...
    'nprint',1000,'resume','../output/sim5_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');

if 0
    load('../output/sim5_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/mcmc_id1.mat')
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
%     Ystd = Y;
%     Ystd(:,1) = Y(:,1) ./ max(Y(:,1));
%     thetree = llike_termnodes(thetree,Ystd);
    Treeplot(thetree)
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,.25,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - wblcdf(xx,5,2);
        plot(xx,yy)
    hold off

    figure()
    get_surv_tree(thetree,Y,X,1000,1,.7,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - wblcdf(xx,1,5);
        plot(xx,yy)
        xlim([0,2])
    hold off

    plot(output.As)
    plot(output.Omegas)
    plot(output.Etas)
    plot(output.Kappas)
    plot(output.As,output.Omegas,'o')    

    [GAM,grid,thetas] = plot_surv(Y,X,thetree,10,50,[],[],[],1);

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


addpath(genpath('../bayes-treesurv-gp'));
dat = readtable('pbc.csv','Delimiter',',','ReadVariableNames',true);
Y = [dat.time, dat.status];
%X = dat(:, [4:20]);
X = dat(:, [4:6, 10, 11, 13, 19]); % based on abstract from https://aasldpubs.onlinelibrary.wiley.com/doi/abs/10.1002/hep.1840100102
% X = dat(:, [4, 13]);

Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc01/', 'saveall', 1, 'seed', 7320, 'n_parallel_temp', 8);
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc02/', 'saveall', 1, 'seed', 7321, 'n_parallel_temp', 8, 'resume', '../output/pbc01/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc03/', 'saveall', 1, 'seed', 7322, 'n_parallel_temp', 8, 'resume', '../output/pbc02/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc04/', 'saveall', 1, 'seed', 7323, 'n_parallel_temp', 8, 'resume', '../output/pbc03/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc05/', 'saveall', 1, 'seed', 7324, 'n_parallel_temp', 8, 'resume', '../output/pbc04/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc06/', 'saveall', 1, 'seed', 7325, 'n_parallel_temp', 8, 'resume', '../output/pbc05/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc07/', 'saveall', 1, 'seed', 7326, 'n_parallel_temp', 8, 'resume', '../output/pbc06/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc08/', 'saveall', 1, 'seed', 7327, 'n_parallel_temp', 8, 'resume', '../output/pbc07/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc09/', 'saveall', 1, 'seed', 7328, 'n_parallel_temp', 8, 'resume', '../output/pbc08/');
Tree_Surv(Y, X, 'burn', 0, 'nmcmc', 1000, 'filepath', '../output/pbc10/', 'saveall', 1, 'seed', 7329, 'n_parallel_temp', 8, 'resume', '../output/pbc09/');

ind = 0;
themax0 = -Inf;
I0 = 0;
for ii=1:10
    load(['../output/pbc', num2str(ii, '%02i'), '/mcmc_id1.mat'])
    [themax, I] = max(output.llike + output.lprior);
    if themax > themax0
        themax0 = themax;
        I0 = I;
        ind = ii;
    end
end
load(['../output/pbc', num2str(ind, '%02i'), '/mcmc_id1.mat'])

thetree = output.Trees{I0};
Treeplot(thetree, Y, X, 0)
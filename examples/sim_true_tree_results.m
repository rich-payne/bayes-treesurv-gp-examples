addpath(genpath('../bayes-treesurv-gp/'));

% Posterior Analysis (and figures)
bname = '../output/sim_true_tree';
n_sims = 50;
n_rep = 10;
alloutput = cell(n_rep,n_sims);
mlpost = zeros(n_sims, 1);
maxrep = zeros(n_sims, 1);
maxtree = cell(n_sims, 1);
% pdraws = cell(n_sims);
seed = 1:50; % data generating seed from simulations: MUST MATCH!
for sim=1:n_sims
    sim
    for rep=1:n_rep
        rep
        fname = strcat([bname, '_simnum_', num2str(sim), '_rep_', num2str(rep), '/mcmc_id1.mat']);
        load(fname)
        alloutput{rep, sim} = output;
        if rep == 1
            [mlpost(sim), I] = max(output.llike + output.lprior);
            maxrep(sim) = 1;
            maxtree{sim} = output.Trees{I};
        else
            [tmpmax, I] = max(output.llike + output.lprior);
            if tmpmax > mlpost(sim)
                mlpost(sim) = tmpmax;
                maxrep(sim) = rep;
                maxtree{sim} = output.Trees{I};
            end
        end
    end
    %[Y, X] = gen_data_true_tree(seed(sim));
    %pdraw = get_surv_tree(maxtree{sim}, Y, X, 10000, 0, X(sim, :), [], .05, []);
    %pdraws{sim} = pdraw;
    % need to save pdraws for each sim
end

% Generate data to evalute quantiles for evaluation (comment when complete)
% [Y, X] = gen_data_true_tree(seed(sim));
[~, X_new, Y_new] = gen_data_true_tree(12345);
n_data = size(X_new, 1);
% quantiles for evaluation of brier scores
% ts = quantile(Y_new, [.1, .25, .5, .75, .9]);
% ts = round(ts, 2);
ts = [.09, .62, 2.04, 4.6, 6.46];

%n_new_data = 100;
%X_new = X_all(1:n_new_data, :);
% y_max = max(Y(:, 1));
% PS COULD BE WRONG (should time be consistent rather than percentiles?)
%ps = [0.1, 0.25, 0.5, 0.75, 0.9]; % percentiles to calculate quantiles
% info for arbitrary survival curve
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
ht = -log(surv);
theta1 = 1;
surv1 = exp(-ht * theta1);
qsurv = zeros(n_data, n_sims, length(ts));
qsurv_lb = qsurv;
qsurv_ub = qsurv;
qtrue = qsurv;
brier = zeros(n_sims, length(ts));

% Bias, MSE, Brier, Coverage
for sim=1:n_sims
    [Y, X] = gen_data_true_tree(seed(sim));
    [~, X_new_brier, Y_new_brier] = gen_data_true_tree(seed(sim) + 100);
    thetree = maxtree{sim};
    [pdraws_all, term_node_ind] = get_surv_tree(thetree, Y, X, 10000, 0, [], ts, .05, []);
    brier(sim, :) = get_brier_score(thetree, Y, X, Y_new_brier, X_new_brier, ts);
    % Bias and MSE
    for ind = 1:n_data
        % loop the following over all sims
        [~, t_ind] = get_termnode(thetree, X_new(ind, :));
        pdraws = pdraws_all{t_ind == term_node_ind};
        qsurv(ind, sim, :) = pdraws.pmean;
        qsurv_lb(ind, sim, :) = pdraws.CI(1, :);
        qsurv_ub(ind, sim, :) = pdraws.CI(2, :);        
        % get true value
        x1 = X_new{ind, 1};
        x2 = X_new{ind, 2};
        if ismember(x2,{'A','B'}) && x1 > 5
            qtrue(ind, sim, :) = 1 - wblcdf(ts, 5,2);
        elseif ismember(x2,{'A'}) && x1 <= 5
            qtrue(ind, sim, :) = 1 - wblcdf(ts, 1,5);
        elseif ismember(x2,{'B'}) && x1 <= 5
            qtrue(ind, sim, :) = 1 - wblcdf(ts, .5,.9);
        elseif ismember(x2,{'C','D'}) && x1 <= 3
            qtrue(ind, sim, :) = 1 - wblcdf(ts, 5,5);
        elseif ismember(x2,{'C','D'}) && x1 > 3 && x1 <= 7
            qtrue(ind, sim, :) = 1 - wblcdf(ts, .5,.5);
        elseif ismember(x2,{'C','D'}) && x1 >  7
            qtrue(ind, sim, :) = interp1(times, surv1, ts);
        end
    end
end

bias = squeeze(mean(qsurv - qtrue, 2, 'omitnan'));
% rmse = squeeze(sqrt(mean(squeeze(mean((qsurv - qtrue), 1, 'omitnan') .^ 2), 1, 'omitnan')))';
rmse = squeeze(mean((qsurv - qtrue) .^ 2, 2, 'omitnan'));
coverage = squeeze(mean(qsurv_lb < qtrue & qtrue < qsurv_ub, 2, 'omitnan'));
brier_mean = mean(brier, 1);

tab = table(ts', brier_mean', 'VariableNames', ["time", "brier"]);

writetable(tab, 'true_tree_results_brier_mean.csv');
writematrix(brier, 'true_tree_results_brier.csv');
writematrix(bias, 'true_tree_results_bias.csv');
writematrix(rmse, 'true_tree_results_rmse.csv');
writematrix(coverage, 'true_tree_results_coverage.csv');


addpath(genpath('../bayes-treesurv-gp/'));

% Posterior Analysis (and figures)
bname = '../output/sim_cox';
n_sims = 50;
n_rep = 10;
% alloutput = cell(n_rep,n_sims);
mlpost = zeros(n_sims, 1);
maxrep = zeros(n_sims, 1);
maxtree = cell(n_sims, 1);
seed = 101:150; % data generating seed from simulations: MUST MATCH!
for sim=1:n_sims
    sim
    for rep=1:n_rep
        rep
        fname = strcat([bname, '_simnum_', num2str(sim), '_rep_', num2str(rep), '/mcmc_id1.mat']);
        load(fname)
        %alloutput{rep, sim} = output;
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
[~, X_new, Y_new] = gen_data_cox(12345);
n_data = size(X_new, 1);
% quantiles for evaluation of brier scores
% ts = quantile(Y_new, [.1, .25, .5, .75, .9]);
% ts = round(ts, 2);
ts = [.11, .36, 1.24, 4.08, 10.36];

% true data generating
betas = [-1 1 2 0 0 -2, -1, 1, 1.5 -1.5]';


qsurv = zeros(n_data, n_sims, length(ts));
qsurv_lb = qsurv;
qsurv_ub = qsurv;
qtrue = qsurv;
brier = zeros(n_sims, length(ts));

% Bias, MSE, Brier, Coverage
for sim=1:n_sims
    [Y, X] = gen_data_cox(seed(sim));
    [~, X_new_brier, Y_new_brier] = gen_data_cox(seed(sim) + 200);
    thetree = maxtree{sim};
    [pdraws_all, term_node_ind] = get_surv_tree(thetree, Y, X, 10000, 0, [], ts, .05, []);
    brier(sim, :) = get_brier_score(thetree, Y, X, Y_new_brier, X_new_brier, ts);
    X_new_array = table2array(X_new);
    % Bias and MSE
    for ind = 1:n_data
        % loop the following over all sims
        [~, t_ind] = get_termnode(thetree, X_new(ind, :));
        pdraws = pdraws_all{t_ind == term_node_ind};
        qsurv(ind, sim, :) = pdraws.pmean;
        qsurv_lb(ind, sim, :) = pdraws.CI(1, :);
        qsurv_ub(ind, sim, :) = pdraws.CI(2, :);        
        
        % survival true value
        nu = X_new_array(ind, :) * betas;
        thehaz = 0.5 * exp(nu); % Hazard function
        qtrue(ind, sim, :) = exp(-thehaz * ts);
    end
end

bias = squeeze(mean(qsurv - qtrue, 2, 'omitnan'));
% rmse = squeeze(sqrt(mean(squeeze(mean((qsurv - qtrue), 1, 'omitnan') .^ 2), 1, 'omitnan')))';
rmse = squeeze(mean((qsurv - qtrue) .^ 2, 2, 'omitnan'));
coverage = squeeze(mean(qsurv_lb < qtrue & qtrue < qsurv_ub, 2, 'omitnan'));
brier_mean = mean(brier, 1);

tab = table(ts', brier_mean', 'VariableNames', ["time", "brier"]);

writetable(tab, 'cox_results_brier_mean.csv');
writematrix(brier, 'cox_results_brier.csv');
writematrix(bias, 'cox_results_bias.csv');
writematrix(rmse, 'cox_results_rmse.csv');
writematrix(coverage, 'cox_results_coverage.csv');


addpath(genpath('../bayes-treesurv-gp/'));

% Posterior Analysis (and figures)
bname = '../output/sim_pred_prog';
n_sims = 50;
n_rep = 10;
mlpost = zeros(n_sims, 1);
maxrep = zeros(n_sims, 1);
maxtree = cell(n_sims, 1);
seed = 201:250; % data generating seed from simulations: MUST MATCH!
for sim=1:n_sims
    sim
    for rep=1:n_rep
        rep
        fname = strcat([bname, '_simnum_', num2str(sim), '_rep_', num2str(rep), '/mcmc_id1.mat']);
        load(fname)
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
[~, X_new, Y_new] = gen_data_pred_prog(12345);
n_data = size(X_new, 1);
% quantiles for evaluation of brier scores
% ts = quantile(Y_new, [.1, .25, .5, .75, .9]);
% ts = round(ts, 2);
% ts = [.55, .81, 1.13, 1.84, 2.50]; % ensure it is up to date.

% true data generating
a_shape = 1;
b_shape = 5;
a_scale = 1;
b_scale = 2;

qsurv = zeros(n_data, n_sims, length(ts));
qsurv_lb = qsurv;
qsurv_ub = qsurv;
qtrue = qsurv;
brier = zeros(n_sims, length(ts));
miss = brier;
brier_cens = brier;
miss_cens = brier;
miss_detail = [];

% Bias, MSE, Brier, Coverage
for sim=1:n_sims
    [Y, X] = gen_data_pred_prog(seed(sim));
    [Y_new_brier_cens, X_new_brier, Y_new_brier] = gen_data_pred_prog(sim + 100); % ensure it is different
    thetree = maxtree{sim};
    [pdraws_all, term_node_ind] = get_surv_tree(thetree, Y, X, 10000, 0, [], ts, .05, []);
    [brier(sim, :), miss(sim, :)] = get_brier_score(thetree, Y, X, Y_new_brier, X_new_brier, ts);
    [brier_cens(sim, :), miss_cens(sim, :), miss_detail0] = get_brier_score_cens(thetree, Y, X, Y_new_brier_cens, X_new_brier, ts);
    miss_detail = [miss_detail; [ones(size(miss_detail0, 1), 1) * sim, miss_detail0]];
    % Bias and MSE
    for ind = 1:n_data
        % loop the following over all sims
        [~, t_ind] = get_termnode(thetree, X_new(ind, :));
        pdraws = pdraws_all{t_ind == term_node_ind};
        qsurv(ind, sim, :) = pdraws.pmean;
        qsurv_lb(ind, sim, :) = pdraws.CI(1, :);
        qsurv_ub(ind, sim, :) = pdraws.CI(2, :);        
        
        % survival true value
        x_trt = X_new{ind, 'trt'};
        x_biomarker = X_new{ind, 'biomarker'};
        w_shape = a_shape + b_shape * x_biomarker;
        w_scale = a_scale + b_scale * x_trt .* x_biomarker;
        qtrue(ind, sim, :) = 1 - wblcdf(ts, w_scale, w_shape);
    end
end

bias = squeeze(mean(qsurv - qtrue, 2, 'omitnan'));
% rmse = squeeze(sqrt(mean(squeeze(mean((qsurv - qtrue), 1, 'omitnan') .^ 2), 1, 'omitnan')))';
rmse = sqrt(squeeze(mean((qsurv - qtrue) .^ 2, 2, 'omitnan')));
coverage = squeeze(mean(qsurv_lb < qtrue & qtrue < qsurv_ub, 2, 'omitnan'));
brier_mean = mean(brier, 1);
miss_detail_tab = array2table(miss_detail, 'VariableNames', {'sim', 'time1', 'time2', 'time3', 'time4', 'time5'});
tab = table(ts', brier_mean', 'VariableNames', ["time", "brier"]);

writetable(tab, 'pred_prog_results_brier_mean.csv');
writematrix(brier, 'pred_prog_results_brier.csv');
writematrix(brier_cens, 'pred_prog_results_brier_cens.csv');
writematrix(miss, 'pred_prog_results_miss.csv');
writematrix(miss_cens, 'pred_prog_results_miss_cens.csv');
writematrix(bias, 'pred_prog_results_bias.csv');
writematrix(rmse, 'pred_prog_results_rmse.csv');
writematrix(coverage, 'pred_prog_results_coverage.csv');
writetable(miss_detail_tab, 'pred_prog_results_missing_detail.csv');

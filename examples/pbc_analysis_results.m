addpath(genpath('../bayes-treesurv-gp/'));

n_kfold = 10;
n_reps = 10;
mlpost = zeros(n_kfold, 1);
maxrep = zeros(n_kfold, 1);
maxtree = cell(n_kfold, 1);
data = readtable("pbc_kfold.csv");
%times = quantile(data.time, [.1, .25, .5, .75, .9]);
times = [617, 1181, 1788, 2691, 3608];
briers = zeros(n_kfold, length(times));
miss = briers;
miss_detail = [];
for fold=1:n_kfold
    fold
    for rep=1:n_reps
        fname = strcat(['../output/pbc', num2str(rep, "%02d"), '_kfold_', num2str(fold), '/mcmc_id1.mat']);
        load(fname)
        if rep == 1
            [mlpost(fold), I] = max(output.llike + output.lprior);
            maxrep(fold) = 1;
            maxtree{fold} = output.Trees{I};
        else
            [tmpmax, I] = max(output.llike + output.lprior);
            if tmpmax > mlpost(fold)
                mlpost(fold) = tmpmax;
                maxrep(fold) = rep;
                maxtree{fold} = output.Trees{I};
            end
        end
    end
    % Calculate brier scores for hold-out datasets
    data_orig = data(data.k_fold ~= fold, :);
    data_holdout = data(data.k_fold == fold, :);
    Y = [data_orig.time, data_orig.status];
    X = data_orig(:, [4:6, 10, 11, 13, 19]);
    Y_new = [data_holdout.time, data_holdout.status];
    X_new = data_holdout(:, [4:6, 10, 11, 13, 19]);
    thetree = fatten_tree(maxtree{fold}, X);
    [briers(fold, :), miss(fold, :), miss_detail0] = get_brier_score_cens(thetree, Y, X, Y_new, X_new, times);
    miss_detail = [miss_detail; [ones(size(miss_detail0, 1), 1) * fold, miss_detail0]];
end
miss_detail_tab = array2table(miss_detail, 'VariableNames', {'fold', 'time1', 'time2', 'time3', 'time4', 'time5'});

writematrix(briers, 'pbc_analysis_results_brier_cens.csv');
writematrix(miss, 'pbc_analysis_results_missing.csv');
writetable(miss_detail_tab, 'pbc_analysis_results_missing_detail.csv');


ts = [617, 1181, 1788, 2691, 3608]; % survival quantiles
n_times = length(ts);
data_all = readtable("pbc_kfold.csv");
data_miss = readtable('pbc_analysis_results_missing_detail.csv');
% n_subj = size(data_ref, 1);
n_folds = 10;
n_params = 1e4;
brier = zeros(n_folds, n_times);
brier_cens = brier;
brier_cens_miss = brier;
surv_hat_indep = [];
% surv_hat_indep = zeros(n_folds, n_subj, n_times);
% surv_hat_ref = surv_hat_indep;
% surv_true_indep = surv_hat_indep;
% surv_true_ref = surv_hat_indep;
% surv_hat_ref_lb = surv_hat_indep;
% surv_hat_ref_ub = surv_hat_indep;
for fold=1:n_folds
    fold
    data = data_all(data_all.k_fold ~= fold, :);
    data_indep = data_all(data_all.k_fold == fold, :);
    X = get_x_pbc(data);
    X_indep = get_x_pbc(data_indep);
    y = data.time;
    censored = 1 - data.status;
    % Get the estimate from the model
    [b, logl, H, stats] = coxphfit(X, y, 'Censoring', censored);
    params = mvnrnd(b, stats.covb, n_params)';
    H2 = [[0, 0]; H];
    % interpolate hazard function (and handle duplicates)
    [c, ~, ic] = unique(H2(:, 1), 'stable');
    val = accumarray(ic, H2(:, 2), [], @mean);
    cumhaz = interp1(c, val, ts);
    % calculate censoring disribution for brier scores
    [f, x] = ecdf(data_all.time, 'censoring', data_all.status);
    % get interpolation, take care of duplicate x values
    [c, ~, ic] = unique(x, 'stable');
    val = accumarray(ic, 1 - f, [], @mean);
    G = griddedInterpolant(c, val, 'previous');
    n_subj = size(data_indep, 1);
    brier_sim = zeros(n_subj, n_times);
    for subj_index=1:n_subj
        samps_indep = zeros(n_params, n_times);
        samps_ref = samps_indep;
        for ii=1:n_params
            bb = params(:, ii);
            cum_hazard_indep = cumhaz * exp(X_indep(subj_index, :) * bb);
            surv_func_indep = exp(-cum_hazard_indep);
            samps_indep(ii, :) = surv_func_indep;
        end
        surv_hat_indep0 = mean(samps_indep);
        surv_hat_indep = [surv_hat_indep; [fold, surv_hat_indep0]];
        % brier score
        for t_ind=1:length(ts)
            if (data_indep.time(subj_index) > ts(t_ind))
                brier_sim(subj_index, t_ind) = (1 - surv_hat_indep0(t_ind)) ^ 2 / G(ts(t_ind));
            elseif (data_indep.time(subj_index) <= ts(t_ind)) && (data_indep.status(subj_index) == 1)
                brier_sim(subj_index, t_ind) = (surv_hat_indep0(t_ind)) ^ 2 / G(data_indep.time(subj_index));
            end % otherwise zero
        end
    end
    brier_cens(fold, :) = mean(brier_sim);
    miss_index = table2array(data_miss(data_miss.fold == fold, 2:6));
    brier_cens_miss(fold, :) = [mean(brier_sim(miss_index(:, 1) == 0, 1)),...
        mean(brier_sim(miss_index(:, 2) == 0, 2)),...
        mean(brier_sim(miss_index(:, 3) == 0, 3)),...
        mean(brier_sim(miss_index(:, 4) == 0, 4)),...
        mean(brier_sim(miss_index(:, 5) == 0, 5))];
end

bias = squeeze(mean(surv_hat_ref - surv_true_ref, 1));
rmse = sqrt(squeeze(mean((surv_hat_ref - surv_true_ref) .^ 2, 1)));
coverage = squeeze(mean(surv_hat_ref_lb < surv_true_ref & surv_true_ref < surv_hat_ref_ub, 1));

% writematrix(brier, 'pred_prog_results_brier_cox_model.csv');
% writematrix(bias, 'pred_prog_results_bias_cox_model.csv');
% writematrix(rmse, 'pred_prog_results_rmse_cox_model.csv');
% writematrix(coverage, 'pred_prog_results_coverage_cox_model.csv');
% writematrix(brier_cens, 'pred_prog_results_brier_cens_cox_model.csv');
% writematrix(brier_cens_miss, 'pred_prog_results_brier_cens_miss_cox_model.csv');

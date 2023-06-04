% Truth vs cox PH
a_shape = 1;
b_shape = 5;
a_scale = 1;
b_scale = 2;
ts = [.55, .81, 1.13, 1.84, 2.50]; % survival quantiles
n_times = length(ts);
data_ref = readtable("../competingmethods/data/pred_prog_12345.csv");
data_miss = readtable('pred_prog_results_missing_detail.csv');
n_subj = size(data_ref, 1);
seed_data = 201:250;
seed_indep = 101:150;
n_sims = 50;
n_params = 1e4;
brier = zeros(n_sims, n_times);
brier_cens = brier;
brier_cens_miss = brier;
surv_hat_indep = zeros(n_sims, n_subj, n_times);
surv_hat_ref = surv_hat_indep;
surv_true_indep = surv_hat_indep;
surv_true_ref = surv_hat_indep;
surv_hat_ref_lb = surv_hat_indep;
surv_hat_ref_ub = surv_hat_indep;
for sim=1:n_sims
    sim
    data = readtable(strcat(['../competingmethods/data/pred_prog_', num2str(seed_data(sim)), '.csv']));
    data_indep = readtable(strcat(['../competingmethods/data/pred_prog_', num2str(seed_indep(sim)), '.csv']));
    X = [data.biomarker, data.trt, data.biomarker .* data.trt];
    X_indep = [data_indep.biomarker, data_indep.trt, data_indep.biomarker .* data_indep.trt];
    X_ref = [data_ref.biomarker, data_ref.trt, data_ref.biomarker .* data_ref.trt];
    y = data.Y_1;
    censored = abs(1 - data.Y_2);
    % Get the estimate from the model
    [b, logl, H, stats] = coxphfit(X, y, 'Censoring', censored);
    params = mvnrnd(b, stats.covb, n_params)';
    H2 = [[0, 0]; H];
    % interpolate hazard function (and handle duplicates)
    [c, ~, ic] = unique(H2(:, 1), 'stable');
    val = accumarray(ic, H2(:, 2), [], @mean);
    cumhaz = interp1(c, val, ts);
    % calculate censoring disribution for brier scores
    [f, x] = ecdf(data.Y_1, 'censoring', data.Y_2);
    % get interpolation, take care of duplicate x values
    [c, ~, ic] = unique(x, 'stable');
    val = accumarray(ic, 1 - f, [], @mean);
    G = griddedInterpolant(c, val, 'previous');
    brier_sim = zeros(n_subj, n_times);
    for subj_index=1:n_subj
        samps_indep = zeros(n_params, n_times);
        samps_ref = samps_indep;
        for ii=1:n_params
            bb = params(:, ii);

            cum_hazard_indep = cumhaz * exp(X_indep(subj_index, :) * bb);
            surv_func_indep = exp(-cum_hazard_indep);
            samps_indep(ii, :) = surv_func_indep;

            cum_hazard_ref = cumhaz * exp(X_ref(subj_index, :) * bb);
            surv_func_ref = exp(-cum_hazard_ref);
            samps_ref(ii, :) = surv_func_ref;
        end
        surv_hat_indep(sim, subj_index, :) = mean(samps_indep);
        surv_hat_ref(sim, subj_index, :) = mean(samps_ref);
        surv_hat_ref_lb(sim, subj_index, :) = quantile(samps_ref, .025);
        surv_hat_ref_ub(sim, subj_index, :) = quantile(samps_ref, .975);
        w_shape_indep = a_shape + b_shape * data_indep.biomarker(subj_index);
        w_scale_indep = a_scale + b_scale * data_indep.trt(subj_index) .* data_indep.biomarker(subj_index);
        surv_true_indep(sim, subj_index, :) = 1 - wblcdf(ts, w_scale_indep, w_shape_indep);
        w_shape_ref = a_shape + b_shape * data_ref.biomarker(subj_index);
        w_scale_ref = a_scale + b_scale * data_ref.trt(subj_index) .* data_ref.biomarker(subj_index);
        surv_true_ref(sim, subj_index, :) = 1 - wblcdf(ts, w_scale_ref, w_shape_ref);
        % brier score
        for t_ind=1:length(ts)
            if (data_indep.Y_1(subj_index, 1) > ts(t_ind))
                brier_sim(subj_index, t_ind) = (1 - surv_hat_indep(sim, subj_index, t_ind)) ^ 2 / G(ts(t_ind));
            elseif (data_indep.Y_1(subj_index, 1) <= ts(t_ind)) && (data_indep.Y_2(subj_index) == 1)
                brier_sim(subj_index, t_ind) = (surv_hat_indep(sim, subj_index, t_ind)) ^ 2 / G(data_indep.Y_1(subj_index));
            end % otherwise zero
        end
    end
    % brier scores for each simulated dataset
    Y_indep = data_indep.Y_true;
    Y_hat = zeros(size(Y_indep, 1), n_times);
    for ind = 1:n_times
        Y_hat(:, ind) = Y_indep > ts(ind);
    end
    brier(sim, :) = mean((Y_hat - squeeze(surv_hat_indep(sim, :, :))) .^ 2);
    brier_cens(sim, :) = mean(brier_sim);
    miss_index = table2array(data_miss(data_miss.sim == sim, 2:6));
    brier_cens_miss(sim, :) = [mean(brier_sim(miss_index(:, 1) == 0, 1)),...
        mean(brier_sim(miss_index(:, 2) == 0, 2)),...
        mean(brier_sim(miss_index(:, 3) == 0, 3)),...
        mean(brier_sim(miss_index(:, 4) == 0, 4)),...
        mean(brier_sim(miss_index(:, 5) == 0, 5))];
end

bias = squeeze(mean(surv_hat_ref - surv_true_ref, 1));
rmse = sqrt(squeeze(mean((surv_hat_ref - surv_true_ref) .^ 2, 1)));
coverage = squeeze(mean(surv_hat_ref_lb < surv_true_ref & surv_true_ref < surv_hat_ref_ub, 1));

writematrix(brier, 'pred_prog_results_brier_cox_model.csv');
writematrix(bias, 'pred_prog_results_bias_cox_model.csv');
writematrix(rmse, 'pred_prog_results_rmse_cox_model.csv');
writematrix(coverage, 'pred_prog_results_coverage_cox_model.csv');
writematrix(brier_cens, 'pred_prog_results_brier_cens_cox_model.csv');
writematrix(brier_cens_miss, 'pred_prog_results_brier_cens_miss_cox_model.csv');

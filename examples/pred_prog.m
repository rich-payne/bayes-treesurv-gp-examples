addpath(genpath('../bayes-treesurv-gp'));

rng(97891, 'twister');
n = 1000;
x_biomarker = rand(n, 1);
x_trt = binornd(1, .5, n, 1);
X = [x_biomarker, x_trt, x_biomarker .* x_trt];
a_shape = 1;
b_shape = 5;
a_scale = 1;
b_scale = 2;
w_shape = a_shape + b_shape * x_biomarker;
w_scale = a_scale + b_scale * x_trt .* x_biomarker;

y_true = wblrnd(w_scale, w_shape, n, 1);
y_cens = exprnd(20, n, 1);
censored = y_cens < y_true;
y = y_true;
y(censored) = y_cens(censored);

X_table = array2table(X(:, [1, 2]), 'VariableNames', {'biomarker', 'trt'});
Y = [y, abs(censored - 1)];
writematrix([Y, X(:, [1, 2])], 'pred_prog.csv');
Tree_Surv(Y, X_table, 'burn', 0, 'nmcmc', 10000, 'filepath', '../output/pred_prog/', 'saveall', 1, 'seed', 1990, 'n_parallel_temp', 8);

% plots
load('../output/pred_prog/mcmc_id1.mat');
[~, I] = max(output.llike + output.lprior);
thetree = output.Trees{I};
Treeplot(thetree, Y, X_table, 0)

% Truth vs cox PH
% Get the estimate from the model vs true hazard
[b, logl, H, stats] = coxphfit(X, y, 'Censoring', censored);
stats.p; % p-values
b; % coefficients
params = mvnrnd(b, stats.covb, 1e5)';
subj_index = 436;
H2 = [[0, 0]; H];
samps = zeros(size(params, 2), size(H2, 1));
for ii=1:size(params, 2)
    bb = params(:, ii);
    hazard = H2(:, 2) * exp(X(subj_index, :) * bb);
    cum_hazard = cumtrapz(H2(:, 1), hazard);
    surv_func = exp(-cum_hazard);
    samps(ii, :) = surv_func;
end
qtls = quantile(samps, [.025, .975]);
muhat = mean(samps);
w_shape_pat = a_shape + b_shape * x_biomarker(subj_index);
w_scale_pat = a_scale + b_scale * x_trt(subj_index) .* x_biomarker(subj_index);
subplot(1, 4, 1);
plot(H2(:, 1), muhat, 'k:', H2(:, 1), qtls(1, :), 'k--', H2(:, 1), qtls(2, :), 'k--')
hold on
plot(H2(:, 1), 1 - wblcdf(H2(:, 1), w_scale_pat, w_shape_pat), 'k');
hold off
title('Cox Proportional Hazard');
ylabel('survival');
xlabel('time');

% BART
subplot(1, 4, 2);
db = readmatrix('../competingmethods/BART/bart_pred.csv');
plot(db(:, 1), db(:, 2), 'k:', db(:, 1), db(:, 3), 'k--', db(:, 1), db(:, 4), 'k--', db(:, 1), db(:, 5), 'k')
title('BART');
ylabel('survival');
xlabel('time');

% Random forest
subplot(1, 4, 3)
rf = readtable('../competingmethods/random_forest/random_forest_pred.csv');
plot(rf.time, rf.mean, 'k:', rf.time, rf.lb, 'k--', rf.time, rf.ub, 'k--');
hold on;
plot(H2(:, 1), 1 - wblcdf(H2(:, 1), w_scale_pat, w_shape_pat), 'k')
hold off;
title('Random Forest');
ylabel('survival');
xlabel('time');

% tree model
subplot(1, 4, 4);
x0 = X_table(subj_index, :);
get_surv_tree(thetree,Y,X_table,10000,1,x0,[],.05,'THM')
hold on;
plot(H2(:, 1), 1 - wblcdf(H2(:, 1), w_scale_pat, w_shape_pat), 'k')
hold off;
xlim([0, 4]);
ylim([0, 1]);
ylabel('survival');
xlabel('time');







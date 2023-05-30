function [brier, miss] = get_brier_score(thetree, Y, X, Y_new, X_new, times)
    % Y_new is uncensored
    [pdraws, term_node_ind] = get_surv_tree(thetree, Y, X, 10000, 0, [], times, .05, []);
    n_x = size(X_new, 1);
    n_times = length(times);
    surv_hat = zeros(n_x, n_times);
    miss = surv_hat;
    for ind = 1:n_x
        [~, t_ind] = get_termnode(thetree, X_new(ind, :));
        pd = pdraws{t_ind == term_node_ind};
        surv_hat(ind, :) = pd.pmean;
        miss(ind, ~isfinite(surv_hat(ind, :))) = 1;
    end
    Y_hat = zeros(n_x, n_times);
    for ind = 1:n_times
        Y_hat(:, ind) = Y_new > times(ind);
    end
    brier = mean((Y_hat - surv_hat) .^ 2, 'omitnan');
    miss = sum(miss);
end
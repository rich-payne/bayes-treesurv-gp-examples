function [brier, miss, miss_complete] = get_brier_score_cens(thetree, Y, X, Y_new, X_new, times)
    % Y_new may be censored
    % from Graf 1999 paper in Stats in Medicine
    [pdraws, term_node_ind] = get_surv_tree(thetree, Y, X, 10000, 0, [], times, .05, []);
    [f, x] = ecdf(Y(:, 1), 'censoring', Y(:, 2));
    % get interpolation, take care of duplicate x values
    [c, ~, ic] = unique(x, 'stable');
    val = accumarray(ic, 1 - f, [], @mean);
    G = griddedInterpolant(c, val, 'previous');
    n_x = size(X_new, 1);
    n_times = length(times);
    brier_all = zeros(n_x, n_times);
    miss = brier_all;
    for ind = 1:n_x
        [~, t_ind] = get_termnode(thetree, X_new(ind, :));
        pd = pdraws{t_ind == term_node_ind};
        surv_hat = pd.pmean;
        miss(ind, ~isfinite(surv_hat)) = 1;
        for t_ind=1:length(times)
            if (Y_new(ind, 1) > times(t_ind))
                brier_all(ind, t_ind) = (1 - surv_hat(t_ind)) ^ 2 / G(times(t_ind));
            elseif (Y_new(ind, 1) <= times(t_ind)) && (Y_new(ind, 2) == 1)
                brier_all(ind, t_ind) = (surv_hat(t_ind)) ^ 2 / G(Y_new(ind));
            end % otherwise zero
        end
    end
    brier = mean(brier_all, 1, 'omitnan');
    miss_complete = miss;
    miss = sum(miss);
end

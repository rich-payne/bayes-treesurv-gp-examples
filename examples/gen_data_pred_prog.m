function [Y, X, Y_uncensored] = gen_data_pred_prog(seed)
    rng(seed)
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

    Y_uncensored = wblrnd(w_scale, w_shape, n, 1);
    y_cens = unifrnd(0, 5, n, 1);
    censored = y_cens < Y_uncensored;
    y = Y_uncensored;
    y(censored) = y_cens(censored);

    X = array2table(X(:, [1, 2]), 'VariableNames', {'biomarker', 'trt'});
    Y = [y, abs(censored - 1)];
end

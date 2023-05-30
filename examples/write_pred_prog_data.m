for seed = [201:250, 101:150, 12345]
    [Y, X, Y_true] = gen_data_pred_prog(seed);
    dat = table(Y, X, Y_true);
    dat = splitvars(dat);
    writetable(dat, strcat('../competingmethods/data/pred_prog_', num2str(seed), '.csv'));
end
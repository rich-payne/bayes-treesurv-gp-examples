for seed = [101:150, 201:250, 12345]
    [Y, X, Y_true] = gen_data_cox(seed);
    dat = table(Y, X, Y_true);
    dat = splitvars(dat);
    writetable(dat, strcat('../competingmethods/data/cox_', num2str(seed), '.csv'));
end
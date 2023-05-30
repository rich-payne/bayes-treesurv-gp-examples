for seed = [1:50, 101:150, 12345]
    [Y, X, Y_true] = gen_data_true_tree(seed);
    dat = table(Y, X, Y_true);
    dat = splitvars(dat);
    writetable(dat, strcat('../competingmethods/data/true_tree_', num2str(seed), '.csv'));
end
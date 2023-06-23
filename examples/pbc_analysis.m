function [] = pbc_analysis(data_file, k_fold, seed, nmcmc)
    if ischar(k_fold)
        k_fold = str2double(k_fold);
    end
    if ischar(seed)
        seed = str2double(seed);
    end
    if ischar(nmcmc)
        nmcmc = str2double(nmcmc);
    end

    dat_full = readtable(data_file,'Delimiter',',','ReadVariableNames',true);
    dat = dat_full(dat_full.k_fold ~= k_fold, :);
    Y = [dat.time, dat.status];
    X = dat(:, [4:6, 10, 11, 13, 19]);
    n_k = 100;
    n_reps = 10;
    n_temp = 8;
    label = strcat(['pbc_kfold_', num2str(k_fold)]);
    run_mcmc(Y, X, nmcmc, n_k, n_reps, label, k_fold, seed, n_temp)
end

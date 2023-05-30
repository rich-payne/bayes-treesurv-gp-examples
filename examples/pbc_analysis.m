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
    %X = dat(:, [4:20]);
    X = dat(:, [4:6, 10, 11, 13, 19]); % based on abstract from https://aasldpubs.onlinelibrary.wiley.com/doi/abs/10.1002/hep.1840100102
    % X = dat(:, [4, 13]);
    n_k = 100;
    n_reps = 10;
    n_temp = 8;
    label = strcat(['pbc_kfold_', num2str(k_fold)]);
    run_mcmc(Y, X, nmcmc, n_k, n_reps, label, k_fold, seed, n_temp)

%     fname1 = strcat(['../output/pbc01_kfold_', num2str(k_fold), '/']);
%     fname2 = strcat(['../output/pbc02_kfold_', num2str(k_fold), '/']);
%     fname3 = strcat(['../output/pbc03_kfold_', num2str(k_fold), '/']);
%     fname4 = strcat(['../output/pbc04_kfold_', num2str(k_fold), '/']);
%     fname5 = strcat(['../output/pbc05_kfold_', num2str(k_fold), '/']);
%     fname6 = strcat(['../output/pbc06_kfold_', num2str(k_fold), '/']);
%     fname7 = strcat(['../output/pbc07_kfold_', num2str(k_fold), '/']);
%     fname8 = strcat(['../output/pbc08_kfold_', num2str(k_fold), '/']);
%     fname9 = strcat(['../output/pbc09_kfold_', num2str(k_fold), '/']);
%     fname10 = strcat(['../output/pbc10_kfold_', num2str(k_fold), '/']);
%     
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname1, 'saveall', 1, 'seed', seed + 1, 'n_parallel_temp', 8);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname2, 'saveall', 1, 'seed', seed + 2, 'n_parallel_temp', 8, 'resume', fname1);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname3, 'saveall', 1, 'seed', seed + 3, 'n_parallel_temp', 8, 'resume', fname2);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname4, 'saveall', 1, 'seed', seed + 4, 'n_parallel_temp', 8, 'resume', fname3);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname5, 'saveall', 1, 'seed', seed + 5, 'n_parallel_temp', 8, 'resume', fname4);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname6, 'saveall', 1, 'seed', seed + 6, 'n_parallel_temp', 8, 'resume', fname5);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname7, 'saveall', 1, 'seed', seed + 7, 'n_parallel_temp', 8, 'resume', fname6);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname8, 'saveall', 1, 'seed', seed + 8, 'n_parallel_temp', 8, 'resume', fname7);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname9, 'saveall', 1, 'seed', seed + 9, 'n_parallel_temp', 8, 'resume', fname8);
%     Tree_Surv(Y, X, 'burn', 0, 'nmcmc', nmcmc, 'filepath', fname10, 'saveall', 1, 'seed', seed + 10, 'n_parallel_temp', 8, 'resume', fname9);
end
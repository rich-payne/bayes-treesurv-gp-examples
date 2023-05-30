function [] = sim_true_tree(seed, sim_number, MCMC, n_reps)
    if ischar(seed)
        seed = str2double(seed);
    end
    if ischar(sim_number)
        sim_number = str2double(sim_number);
    end
    if ischar(MCMC)
        MCMC = str2double(MCMC);
    end
    if ischar(n_reps)
        n_reps = str2double(n_reps);
    end
    
    [Y, X] = gen_data_true_tree(seed);
    n_k = 100;
    label = 'sim_true_tree';
    n_temp = 8;
    run_mcmc(Y, X, MCMC, n_k, n_reps, label, sim_number, seed, n_temp);
end

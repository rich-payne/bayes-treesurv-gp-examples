function [] = run_mcmc(Y, X, MCMC, n_k, n_reps, label, sim_number, seed, n_temp)
    SEEDS = seed .* (1:n_reps);
    cntr = 1;
    for rep = 1:n_reps
        basename = ['../output/', label, '_simnum_', num2str(sim_number), '_rep_'];
        fname = strcat([basename, num2str(rep), '/']);
        % Run the code
        if rep > 1
            fname_resume = strcat([basename,num2str(rep-1),'/']);
            Tree_Surv(Y,X,'nmcmc',MCMC,'burn',0,'filepath',fname,...
                'seed',SEEDS(cntr),'bigK',n_k,'nprint',1000,'saveall',1,...
                'swapfreq',10,'resume',fname_resume,...
                'n_parallel_temp', n_temp);
        else
            Tree_Surv(Y,X,'nmcmc',MCMC,'burn',0,'filepath',fname,...
                'seed',SEEDS(cntr),'bigK',n_k,'nprint',1000,'saveall',1,...
                'swapfreq',10,'n_parallel_temp', n_temp);
        end
        cntr = cntr + 1;
    end
    disp(strcat(['Simulation ', label, ' #', num2str(sim_number), ' is complete!']));
end
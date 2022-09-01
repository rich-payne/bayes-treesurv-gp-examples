function [] = sim_true_tree(seed, sim_number)
    addpath(genpath('../bayes-treesurv-gp/'));
    [Y, X] = gen_data_true_tree(seed);

    % SETTINGS FOR THE USER TO CHANGE
    lowrep = 1;
    highrep = 5;

    MCMC = 10000; % MCMC iterations
    Kvals = 20;
    SEEDS = seed .* (1:highrep); % Seed for K = 20
    cntr = 1;
    for rep = lowrep:highrep
        basename = ['../output/sim_true_tree_simnum_', num2str(sim_number), '_rep_'];
        fname = strcat([basename, num2str(rep), '/']);
        % Run the code
        if rep > 1
            fname_resume = strcat([basename,num2str(rep-1),'/']);
            Tree_Surv(Y,X,'nmcmc',MCMC,'burn',0,'filepath',fname,...
                'seed',SEEDS(cntr),'bigK',Kvals,'nprint',1000,'saveall',1,...
                'swapfreq',10,'resume',fname_resume);
        else
            Tree_Surv(Y,X,'nmcmc',MCMC,'burn',0,'filepath',fname,...
                'seed',SEEDS(cntr),'bigK',Kvals,'nprint',1000,'saveall',1,...
                'swapfreq',10);
        end
        cntr = cntr + 1;
    end
end

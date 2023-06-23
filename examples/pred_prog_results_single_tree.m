% Posterior Analysis (and figures)
bname = '../output/sim_pred_prog';
sim = 1;
n_rep = 10;
seed = 201:250; % data generating seed from simulations
for rep=1:n_rep
    rep
    fname = strcat([bname, '_simnum_', num2str(sim), '_rep_', num2str(rep), '/mcmc_id1.mat']);
    load(fname)
    if rep == 1
        [mlpost, I] = max(output.llike + output.lprior);
        maxrep = 1;
        maxtree = output.Trees{I};
    else
        [tmpmax, I] = max(output.llike + output.lprior);
        if tmpmax > mlpost(sim)
            mlpost = tmpmax;
            maxrep = rep;
            maxtree = output.Trees{I};
        end
    end
end

[Y, X] = gen_data_pred_prog(seed(sim));
Treeplot(maxtree, Y, X, [])

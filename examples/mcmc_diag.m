function mcmc_diag(output)
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I}; 
    Treeplot(thetree)
end
function compare_cens(alloutput,ns,k_ind,nrep) 
    % ns: vector of sample sizes
    % k_ind: index on grid size (K)
    % Compare censoring (assumes <= 3 censoring cases) for trees
    for ii=1:length(ns)
        for jj=1:size(alloutput,2)
            % Find the best tree over all reps
            bestoutput = [];
            maxlpost = -Inf;
            for ll=1:nrep
                output = alloutput{ii,jj,k_ind,ll};
                [maxlpostmp,I] = max(output.llike + output.lprior);
                if maxlpostmp > maxlpost
                    bestoutput = output;
                    maxlpost = maxlpostmp;
                    thetree = output.Trees{I};
                end
            end
            Treeplot(thetree)
            title(strcat(['N = ',num2str(ns(ii)),', Cens = ',num2str(jj)]));
        end
    end
end
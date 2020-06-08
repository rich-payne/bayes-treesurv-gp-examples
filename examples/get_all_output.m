function out = get_all_output(simnum,ns,ncens,Ks,nrep,id)
    out = cell(length(ns),ncens,length(Ks),nrep);
    for ii=1:length(ns)
        for jj=1:ncens
            for kk=1:length(Ks)
                for ll=1:nrep
                    fname = strcat(['../output/sim',num2str(simnum),...
                        'final_N',num2str(ns(ii)),'_cens',num2str(jj),...
                        '_K',num2str(Ks(kk)),'_rep',num2str(ll),'/mcmc_id',...
                        num2str(id),'.mat']);
                    load(fname);
                    out{ii,jj,kk,ll} = output;
                end
            end
        end
    end
end
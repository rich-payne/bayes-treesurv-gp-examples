addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

% Generate some data from exponential proportional hazards model
ns = [1000,500,100]; % sample sizes
censparms = [50,25,2.5]; % censoring parameters
DATA = []; % struct holding all the data

rng(351199)
cntr = 1; 
for ii=1:length(ns)
    for jj=1:length(censparms)
        n = ns(ii);
        p = 5;
        X = rand(n,p);
        % Add in a categorical variable
        X = [X, binornd(1,.5,n,1)];
        beta = 2; % Hazard rate over time (exponential hazard, scale parameter)
        betas = [-1 1 2 0 0 -2]';
        nu = X * betas;
        thehaz = (1/beta) * exp(nu); % Hazard function
        % Therefore each individual comes from a different exponential failure
        % distribution.
        % Failure times
        y = exprnd(1./thehaz,n,1); % the new beta is 1/hazard
        % Add some censoring
        cens = exprnd(censparms(jj),n,1);
        ind = y <= cens;
        disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (n))));
        y(~ind) = cens(~ind);
        ds = double(ind);
        Y = [y,ds];
        X = array2table(X);
        DATA(ii,jj).Y = Y;
        DATA(ii,jj).X = X;
        cntr = cntr + 1;
    end
end

parpool(8);
% SETTINGS FOR THE USER TO CHANGE
lowrep = 1;
highrep = 10;
simnum = 7;

MCMC = [10000]; % MCMC iterations for the small and large K respectively
Kvals = [20];
nn = length(Kvals) * numel(DATA) * highrep;
SEEDS = simnum .* (1:nn); % Seed for K = 20
cntr = 1;
for rep = lowrep:highrep
    for kk=1:length(Kvals) % two different K values (potentially)
        for ii=1:length(ns)
            for jj=1:length(censparms)
                if (rep == 9 && ns(ii) == 500 && jj == 3) || ...
                   (rep == 9 && ns(ii) < 500) || ...
                   rep == 10
                    % Get Data
                    Y = DATA(ii,jj).Y;
                    X = DATA(ii,jj).X;
                    % get filnames
                    basename = strcat(['../output/sim',num2str(simnum),'final_N',num2str(ns(ii)),'_cens',num2str(jj),'_K',num2str(Kvals(kk)),'_rep']);
                    fname = strcat([basename,num2str(rep),'/']);
                    % Run the code
                    if rep > 1
                        fname_resume = strcat([basename,num2str(rep-1),'/']);
                        Tree_Surv(Y,X,'nmcmc',MCMC(kk),'burn',0,'filepath',fname,...
                            'seed',SEEDS(cntr),'bigK',Kvals(kk),'nprint',1000,'saveall',1,...
                            'swapfreq',10,'resume',fname_resume);
                    else
                        Tree_Surv(Y,X,'nmcmc',MCMC(kk),'burn',0,'filepath',fname,...
                            'seed',SEEDS(cntr),'bigK',Kvals(kk),'nprint',1000,'saveall',1,...
                            'swapfreq',10);
                    end
                end
                cntr = cntr + 1;
            end
        end
    end
end

if 0
    load('../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/mcmc_id1.mat')
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    Treeplot(thetree)
    get_surv_tree(thetree,Y,X,1000,1,[],[])
    
    x0 = [.5 .1 .8 .5 .5];
    res = get_surv_tree(thetree,Y,X,1000,1,x0,[]);
    x0hazard = 1/beta * exp(x0 * betas);
    orig_grid = res.ystar*max(Y(:,1));
    hold on;
    plot(orig_grid,1-expcdf(orig_grid,1/x0hazard))
    hold off;
    
    
    
    % Subset of the data wth 500 points
    load('../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    
    % Subset of the data wth 100 points
    load('../output/sim7_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub2);
    Treeplot(thetree)
    
    
    
    
    
    
    
    % Quick check on something...
    Xtmp = X{(X{:,2} <= .22947) & (X{:,3} > .76316),:};
    allinvhazards = 1./(1/beta * exp(Xtmp *  betas));
    allsurvs = zeros(size(Xtmp,1),length(orig_grid));
    for ii=1:size(Xtmp,1)
        allsurvs(ii,:) = 1 - expcdf(orig_grid,allinvhazards(ii));
    end
    plot(allsurvs')
    avgsurv = mean(allsurvs);
    hold on
    plot(res.ystar,avgsurv)
    hold off
    
    

    % Ensemble method...
    llikes = output.llike;
    [ul,ui] = unique(llikes);
    %[ul_sorted,uis] = sort(ul);
    a = max(ul);
    tree_probs = exp(ul - (a + log(sum(exp(ul - a))))); % in smallest to largest order
    nmc = 10000; % number of monte carlo runs;
    x0 = array2table([.1 .1 .9 .9 .5]); % Value of the covariate
    % sample the trees
    treeind = randsample(length(tree_probs),nmc,true,tree_probs);
    %treeind_u = unique(treeind);
    tograph = 0;
    ystar = linspace(.01,.99,102);
    pdraws = zeros(nmc,length(ystar));
    theind = 1;
    for ii=1:length(treeind)
        ndraw = sum(treeind == ii);
        if ndraw > 0
            tmp = get_surv_tree(output.Trees{ui(ii)},Y,X,ndraw,tograph,x0,ystar);
            pdraws(theind:(theind + ndraw -1),:) = tmp.surv;
            theind = theind + ndraw;
        end
    end
    % Now plot posterior draws...
    pmean = mean(pdraws);
    qtiles = quantile(pdraws,[.025,.975],1);
    orig_grid = ystar*max(Y(:,1));
    plot(orig_grid,pmean,'-k',orig_grid,qtiles(1,:),'--k',orig_grid,qtiles(2,:),'--k')
    % True survival function
    x0hazard = 1/beta * exp(x0{1,:} * betas);
    hold on;
    plot(orig_grid,1-expcdf(orig_grid,1/x0hazard))
    hold off;
end









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
        DATA(ii,jj).X = X; DATA(ii,jj).mean = 1 ./ thehaz;
        cntr = cntr + 1;
    end
end

parpool(8);
% SETTINGS FOR THE USER TO CHANGE
% SETTINGS FOR THE USER TO CHANGE
% Original settings
% lowrep = 1;
% highrep = 10;
% simnum = 7;
% Updated settings when error occurred
lowrep = 2;
highrep = 10;
simnum = 7;


MCMC = [10000]; % MCMC iterations for the small and large K respectively
Kvals = [100];
nn = length(Kvals) * numel(DATA) * highrep;
SEEDS = Kvals .* simnum .* (1:nn); % Seed for K = 100
cntr = 1;
for rep = lowrep:highrep
    for kk=1:length(Kvals) % two different K values (potentially)
        for ii=1:length(ns)
            for jj=1:length(censparms)
                if (rep == 8 && ns(ii) == 100 && jj == 3) || rep >= 9
                    % Get Data
                    Y = DATA(ii,jj).Y;
                    X = DATA(ii,jj).X;
                    % get filnames
                    basename = strcat(['../output/sim',num2str(simnum),'final_bigk_N',num2str(ns(ii)),'_cens',num2str(jj),'_K',num2str(Kvals(kk)),'_rep']);
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
    % Posterior Analysis (and figures) of the first result
    bname = '../output/sim7final_bigk_N1000_cens1_K100_rep';
    nrep = 10; % change this later to 10!!!!
    alloutput = cell(nrep,1);
    mlpost = [];
    maxrep = [];
    for ii=1:nrep
        ii
        fname = strcat([bname, num2str(ii), '/mcmc_id1.mat']);
        load(fname)
        alloutput{ii} = output;
        if ii == 1
            mlpost = max(output.llike + output.lprior);
            maxrep = 1;
        else
            tmpmax = max(output.llike + output.lprior);
            if tmpmax > mlpost
                mlpost = tmpmax;
                maxrep = ii;
            end
        end
    end
    output = alloutput{maxrep};
    [~,I] = max(output.llike + output.lprior);
    Y = DATA(1,1).Y;
    X = DATA(1,1).X;
    themean = DATA(1,1).mean;
        % Mallows tree selection
        mtrees = mallow(output);
        mind = 10;
        thetree = mtrees{mind};
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4];
    saveas(fig,'../figs/sim7tree.jpg')
    print('../figs/sim7tree','-depsc')
    % Treeplot(thetree,Y,X)
    
    
    get_surv_tree(thetree,Y,X,1000,1,[],[])
    
    % Obtain which terminal node index each datapoint belongs to
    nodeindex = zeros(size(X,1),1);
    for ii=1:size(X,1)
        [~,tmpind] = get_termnode(thetree,X(ii,:));
        nodeindex(ii) = tmpind;
    end
    uinds = unique(nodeindex); % Unique indices
    ntermnodes = length(uinds); % number of terminal nodes
    % Obtain data so we can plot each terminal node
    Xplot_ind = [];
    for ii=1:ntermnodes
        for jj=1:size(X,1)
            if nodeindex(jj) == uinds(ii)
                Xplot_ind(ii) = jj;
            	continue;
            end
        end
    end
    Xplot = X(Xplot_ind,:);
    % subplot dimensions
    nrow = 2;
    ncol = 4; 
    xlimmax = [15, 4, 20, 7, 50, 12, 10,50,50,50,50];
    figure()
    for ii=1:ntermnodes
        subplot(nrow,ncol,ii)
        get_surv_tree(thetree,Y,X,10000,1,Xplot(ii,:),[],.05,[])
        % Plot the median survival function of the observed data
        node = thetree.Allnodes{uinds(ii)};
        med_mean = median(themean(node.Xind));
        xx = linspace(.01,max(Y(:,1)),10000);
        hold on;
            plot(xx,1-expcdf(xx,med_mean))
            xlim([0,xlimmax(ii)]);
            title('');
            xlabel('t');
            ylabel('S(t)');
        hold off
        % Plot Kaplan Meier Curve
        hold on
            ypart = Y(node.Xind,1);
            cens = Y(node.Xind,2) == 0;
            [f,x,flo,fup] = ecdf(ypart,'censoring',cens);
            plot(x,1-f,':r')
            plot(x,1-flo,'--r');
            plot(x,1-fup,'--r');         
        hold off
    end
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4];
    saveas(fig,'../figs/sim7survival.jpg')
    print('../figs/sim7survival','-depsc')
                
    % Posterior analysis for all 9 simulations
    n_num = 2;     % index on sample sizes
    cens_num = 3;  % index on censoring rates
    K = 100;
    nreps = 10;
    basename = '../output/sim7final_bigk';
    
    Y = DATA(n_num,cens_num).Y;
    X = DATA(n_num,cens_num).X;
    ns = [1000,500,100];
    fname0 = strcat([basename,'_N',num2str(ns(n_num)),...
        '_cens',num2str(cens_num),...
        '_K',num2str(K),...
        '_rep']);
    % Read in all of the reps
    alloutput = cell(nreps,1);
    for ii=1:nreps
        fname = strcat([fname0,num2str(ii),'/mcmc_id1.mat']);
        load(fname)
        alloutput{ii} = output;
    end
    [mtrees,llike,whichout,mtrees_lpost,lpost,whichout_lpost] = mallow(alloutput{1},alloutput{2},alloutput{3},alloutput{4},alloutput{5},...
        alloutput{6},alloutput{7},alloutput{8},alloutput{9},alloutput{10});
    [~,I] = max(lpost); % choose largest posterior tree
    thetree = mtrees_lpost{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree)   
    [~,I2] = max(llike);
    I2 % llike selection of number of termnodes
    
    
    
    
    
    
    
    load('../output/sim7final_bigk_N1000_cens1_K100_rep1/mcmc_id1.mat')
    Y = DATA(1,1).Y;
    X = DATA(1,1).X;
    
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









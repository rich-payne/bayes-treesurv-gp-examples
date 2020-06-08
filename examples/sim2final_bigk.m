addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

% Create my own survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
% plot(times,surv)
ht = -log(surv);

% Data from first function
theta1 = 1;
surv1 = exp(-ht*theta1);
theta2 = 2;
surv2 = exp(-ht*theta2);
% plot(times,surv1);
% hold on;
% plot(times,surv2);
% hold off;

ns = [1000,500,100]; % sample sizes
censparms = [50,25,5]; % censoring parameters
DATA = []; % struct holding all the data

rng(559)
cntr = 1; 
for ii=1:length(ns)
    for jj=1:length(censparms)
        n1 = ns(ii)/2;
        n2 = ns(ii)/2;
        X = table([rand(n1,1)/2; rand(n2,1)/2 + .5]);
        y = [interp1(surv1,times,rand(n1,1)); interp1(surv2,times,rand(n2,1))];
        % Random censoring function.
        %cens = rand(n1 + n2,1) + 5; % Original Censoring...
        cens = exprnd(censparms(jj),n1+n2,1);
        ind = y < cens;
        disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (n1 + n2))));
        ds = double(ind);
        y(~ind) = cens(~ind);
        Y = [y,ds];
        DATA(ii,jj).Y = Y;
        DATA(ii,jj).X = X;
        cntr = cntr + 1;
    end
end

% % GENERATE TEST DATASETS
% rng(18351)
% cntr = 1; 
% DATA_TEST = [];
% for ii=1:length(ns)
%     for jj=1:length(censparms)
%         n1 = ns(ii)/2;
%         n2 = ns(ii)/2;
%         X = table([rand(n1,1)/2; rand(n2,1)/2 + .5]);
%         y = [interp1(surv1,times,rand(n1,1)); interp1(surv2,times,rand(n2,1))];
%         % Random censoring function.
%         %cens = rand(n1 + n2,1) + 5; % Original Censoring...
%         cens = exprnd(censparms(jj),n1+n2,1);
%         ind = y < cens;
%         disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (n1 + n2))));
%         ds = double(ind);
%         y(~ind) = cens(~ind);
%         Y = [y,ds];
%         DATA_TEST(ii,jj).Y = Y;
%         DATA_TEST(ii,jj).X = X;
%         cntr = cntr + 1;
%     end
% end


parpool(8);
% SETTINGS FOR THE USER TO CHANGE
lowrep = 1;
highrep = 10;
simnum = 2;

MCMC = [10000]; % MCMC iterations for the small and large K respectively
Kvals = [100];
nn = length(Kvals) * numel(DATA) * highrep;
SEEDS = Kvals .* simnum .* (1:nn); % Seed for K = 100
cntr = 1;
for rep = lowrep:highrep
    for kk=1:length(Kvals) % two different K values (potentially)
        for ii=1:length(ns)
            for jj=1:length(censparms)
                if ns(ii) == 100 && jj == 3 && rep == 10  % DELETE! TEMPORARY IF STATEMENT
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
    bname = '../output/sim2final_bigk_N1000_cens1_K100_rep';
    nrep = 10;
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
    figure()
    [mtrees,mlikes] = mallow(output);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree);
    thenode = thetree.Allnodes{1};
    ruleval  = thenode.Rule{2};
    size(X{X{:,1} < .5 & X{:,1} > ruleval,1},1) % number of misclassified
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4];
    saveas(fig,'../figs/sim2tree.jpg')
    print('../figs/sim2tree','-depsc')
    
    figure()
    subplot(1,2,1)
    test = get_surv_tree(thetree,Y,X,10000,1,.25,[],.05,[])
    %y_max = max(Y(:,1));
    hold on
        plot(times,surv1)
        xlabel('t');
        ylabel('S(t)');
        title('');
    hold off
    
    subplot(1,2,2)
    get_surv_tree(thetree,Y,X,10000,1,.75,[],.05,[])
    hold on
        plot(times,surv2);
        xlabel('t');
        ylabel('S(t)');
        title('');
    hold off
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4];
    saveas(fig,'../figs/sim2survival.jpg')
    print('../figs/sim2survival','-depsc')
    
    
    
    % Posterior analysis for all 9 simulations
    n_num = 3;     % index on sample sizes
    cens_num = 3;  % index on censoring rates
    K = 100;
    nreps = 10;
    basename = '../output/sim2final_bigk';
    
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
   
    figure()
    subplot(1,2,1)
    get_surv_tree(thetree,Y,X,1000,1,.25,[])
    hold on
    plot(times,surv1,'r')
    hold off
    subplot(1,2,2)
    get_surv_tree(thetree,Y,X,1000,1,.75,[])
    hold on
    plot(times,surv2,'r');
    hold off
    
    % NOTES
    % N=1000, all three censoring find the partition and models it well
    % N=500, cens=1 doesn't find the right partition 
    % N=500, cens=2 finds the partition and models it well
    % N=500, cens=3 log-posterior favors no partition, llike favors
    %                a partition with 3 elements, which models things well.
    %               See follow-up code.
    % N=100, cens=1, no splitting favored at all
    % N=100, cens=2, log-partition favors no partition.  llike-favors  
    %                two partitions.  Models well, but with wider CI.
    %                See code below.
    % N=100, cens=3, log-partition favors no partition.  llike favors
    %                three partitions, models well but widely.  See below.
   
    % FOLLOW UP CODE
    if n_num == 2 & cens_num == 3
        % Plot the highest likelihood tree (has a tree size of 3)
        [~,I] = max(llike);
        thetree = mtrees{I};
        thetree = fatten_tree(thetree,X);
        Treeplot(thetree)
        figure()
        subplot(2,2,1)
        get_surv_tree(thetree,Y,X,1000,1,.25,[])
        hold on
        plot(times,surv1,'r')
        hold off
        subplot(2,2,2)
        get_surv_tree(thetree,Y,X,1000,1,.48,[])
        hold on
        plot(times,surv1,'r')
        hold off
        subplot(2,2,3)
        get_surv_tree(thetree,Y,X,1000,1,.75,[])
        hold on
        plot(times,surv2,'r');
        hold off
    elseif n_num == 3 & cens_num == 2
        % Plot the highest likelihood tree (has a tree size of 3)
        [~,I] = max(llike);
        thetree = mtrees{I};
        thetree = fatten_tree(thetree,X);
        Treeplot(thetree)
        figure()
        subplot(1,2,1)
        get_surv_tree(thetree,Y,X,1000,1,.25,[])
        hold on
        plot(times,surv1,'r')
        hold off
        subplot(1,2,2)
        get_surv_tree(thetree,Y,X,1000,1,.75,[])
        hold on
        plot(times,surv2,'r');
        hold off
    elseif n_ind == 3 && cens_ind == 3
        % Plot the highest likelihood tree (has a tree size of 3)
        [~,I] = max(llike);
        thetree = mtrees{I};
        thetree = fatten_tree(thetree,X);
        Treeplot(thetree)
        figure()
        subplot(2,2,1)
        get_surv_tree(thetree,Y,X,1000,1,.25,[])
        hold on
        plot(times,surv1,'r')
        hold off
        subplot(2,2,2)
        get_surv_tree(thetree,Y,X,1000,1,.75,[])
        hold on
        plot(times,surv2,'r')
        hold off
        subplot(2,2,3)
        get_surv_tree(thetree,Y,X,1000,1,.9,[])
        hold on
        plot(times,surv2,'r');
        hold off
    end
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    load('../output/sim2final_N1000_cens1_K20_rep1/mcmc_id8.mat');
    [~,I] = max(output.llike + output.lprior);
    Y = DATA(1,1).Y;
    X = DATA(1,1).X;
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    figure()
    get_surv_tree(thetree,Y,X,1000,1,.25,[])
    %y_max = max(Y(:,1));
    hold on
    plot(times,surv1,'r')
    hold off
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,.75,[])
    hold on
    plot(times,surv2,'r');
    hold off
    
    
    
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1) ./ max(Ystd(:,1));
    nugget = 0;
    l = .01;
    tau =.1;
    EB = 1;
    K = 20;
    [marg_y,res] = get_marginal(Ystd,K,[],1e-9,tau,l,nugget,EB);
    [marg_y1,res1] = get_marginal(Ystd(X{:,1} < .5,:),K,[],1e-9,tau,l,1e-9,1);
    [marg_y2,res2] = get_marginal(Ystd(X{:,1} >= .5,:),K,[],1e-9,tau,l,1e-9,1);
    marg_y
    marg_y1 + marg_y2
    figure()
    out = get_surv(Y,res1,1000,1,[],.05);
    newgrid = out.ystar * max(Y(:,1));
    hold on
        plot(newgrid,1 - expcdf(out.ystar*max(Y(:,1)),1/2),'-r')
    hold off

    figure()
    out = get_surv(Y,res2,1000,1,[],.05);
    newgrid = out.ystar * max(Y(:,1));
    hold on
        plot(newgrid,1 - expcdf(newgrid,100),'-r')
    hold off
    
    
    
    
    
    
    % Read in all the data
    nrep = 1;
    simnum = 2;
    
    ncens = length(censparms);
    K = [20];
    id = 1;
    alloutput = get_all_output(simnum,ns,ncens,K,nrep,id);
    
    % Compare censoring
    k_ind = 1; % index for K = 20;
    compare_cens(alloutput,ns,k_ind,nrep)
            
                
                
            output = alloutput(ii,jj);
         
            
            %subplot(1,3,jj);
            %suptitle(strcat(['N = ',num2str(ns(ii))]));
            
    
    
    
    simnum = 2;
    K = 20
    out = cell(length(ns),ncens,nrep);
    for ii=1:length(ns)
        for jj=1:ncens
            for kk=1:nrep
                fname = strcat(['../output/sim',num2str(simnum),...
                    'final_N',num2str(ns(ii)),'_cens',num2str(jj),...
                    '_K',num2str(K),'_rep',num2str(rep),'mcmc_id',...
                    num2str(kk),'.mat']);
                load(fname);
                out(ii,jj,kk) = output;
            end
        end
    end
    
    load('../output/sim2final_N1000_cens1_K20_rep1/mcmc_id1.mat')
    output11 = output;
    Y11 = DATA(1,1).Y;
    X11 = DATA(1,1).X;
    
    load('../output/sim2final_N1000_cens2_K20_rep1/mcmc_id1.mat')
    output12 = output;
    Y12 = DATA(1,2).Y;
    X12 = DATA(1,2).X;
    
    load('../output/sim2final_N1000_cens3_K20_rep1/mcmc_id1.mat')
    output13 = output;
    Y13 = DATA(1,3).Y;
    X13 = DATA(1,3).X;
    
    % Compare final trees for large data
    mcmc_diag(output11)
    mcmc_diag(output12)
    mcmc_diag(output13)
    
    
    
    
    
    %tabulate(output.treesize)

    figure()
    get_surv_tree(thetree,Y,X,1000,1,.25,[])
    %y_max = max(Y(:,1));
    hold on
    plot(times,surv1,'r')
    hold off
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,.75,[])
    hold on
    plot(times,surv2,'r');
    hold off
    
    
    
    % Subset of the data wth 250 points
    load('../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    %tabulate(output.treesize)
    figure()
    get_surv_tree(thetree,Ysub1,Xsub1,10000,1,.25,[])
    hold on
    plot(times,surv1,'r')
    hold off
    
    % Subset of the data wth 100 points
    load('../output/sim2_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub2);
    Treeplot(thetree)
    figure()
    get_surv_tree(thetree,Ysub2,Xsub2,10000,1,.25,[])
    %y_max = max(Y(:,1));
    hold on
    plot(times,surv1,'r')
    hold off
    
    
    
    
    
    
    % Test out a fixed grid!
    
    
    truetree = thetree;
    truetree.Allnodes{2}.Rule(2) = {.5};
    truetree = fatten_tree(truetree,X);
    truetree = llike_termnodes(truetree,Y);
    
    
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1) ./ max(Ystd(:,1));
    ysub1 = Ystd(X{:,1} < .5,:);
    ysub2 = Ystd(X{:,1} > .5,:);
    ysub3 = Ystd(X{:,1} < .31786,:);
    ysub4 = Ystd(X{:,1} > .31786,:);
    
    % ysub(:,1) = ysub(:,1) ./ max(Y(:,1));
    K = 100;
    s = linspace(0,1.01,K); % Fixed grid
    eps = 1e-10;
    tau = 1;
    l = .01;
    nugget = 1e-10;
    EB = 1;
    [marg_y1,res1] = get_marginal(ysub1,K,s,eps,tau,l,nugget,EB);
    [marg_y2,res2] = get_marginal(ysub2,K,s,eps,tau,l,nugget,EB);
    [marg_y3,res3] = get_marginal(ysub3,K,s,eps,tau,l,nugget,EB);
    [marg_y4,res4] = get_marginal(ysub4,K,s,eps,tau,l,nugget,EB);
    marg_y1 + marg_y2
    marg_y3 + marg_y4
    
    
    get_surv(Y,res,1000,1,[],.05)
    hold on
      plot(times,surv1,'r')
    hold off
    
    
    ysub = Y(X{:,1} > .5,:);
    ysub(:,1) = ysub(:,1) ./ max(Y(:,1));
    K = 20;
    s = [];
    eps = 1e-10;
    tau = 1;
    l = 1;
    mu = 0;
    nugget = 1e-10;
    EB = 1;
    [marg_y,res2] = get_marginal(ysub,K,s,eps,tau,l,mu,nugget,EB);
    get_surv(max(Y(:,1)),res2,1000,1,[])
    hold on
      plot(times,surv2,'r')
    hold off
    


    plot(output.As)
    plot(output.Omegas)

    [GAM,grid,thetas] = plot_surv(Y,X,thetree,10,50,[],[],[],1000);
    [GAM2,grid2,thetas2] = plot_surv(Y,X,thetree,10,50,[],10,[],1000);
    [GAM3,grid3,thetas3] = plot_surv(Y,X,thetree,10,50,[],1,[],1);

    ii = 2;
    cumhaz = squeeze(GAM(ii,:,:));
    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,min(exp(-cumhaz)))
    plot(grid,max(exp(-cumhaz)))
    plot(times,surv1)
    plot(times,surv2)
    hold off


    ii = 2;
    cumhaz2 = squeeze(GAM2(ii,:,:));
    plot(grid,mean(exp(-cumhaz2)))
    hold on
    plot(grid,min(exp(-cumhaz2)))
    plot(grid,max(exp(-cumhaz2)))
    plot(times,surv1)
    plot(times,surv2)
    hold off



    plot(grid,exp(-cumhaz(2,:)))

    options = optimoptions('fminunc','Display','off','Algorithm','quasi-newton');
    [~,lmax] = fminunc(@(theta) log_ptheta_t(theta,Y,thetree.a,thetree.omega,thetree.theta_shape,thetree.theta_rate),1,options);
    afunc = @(a) ptheta_t(30,Y,a,1,1,.001,lmax+1000);
    fplot(afunc,[0 25])
    fplot(@(theta) ptheta_t(theta,Y,1000,1,1,.001,lmax + 1000),[.01 5])

    afunc(10)

end            
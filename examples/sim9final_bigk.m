% Same as sim8 but we have changed one of the Weibull distributions

addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

% Create my own survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
%plot(times,surv)
ht = -log(surv);
% Data from first function
theta1 = 1;
surv1 = exp(-ht*theta1);

% CART example (Figure 2)
ns = [1000,500,100]; % sample sizes
censparms = [50,25,5]; % censoring parameters
DATA = []; % struct holding all the data

rng(3199771)
cntr = 1; 
for ii=1:length(ns)
    for jj=1:length(censparms)
        N = ns(ii);
        x1 = unifrnd(0,10,N,1);
        grps = {'A','B','C','D'};
        x2 = randsample(grps,N,true)';
        dat = table(x1,x2);
        y = zeros(N,1);
        I1 = ismember(x2,{'A','B'}) & x1 > 5;
        y(I1) = wblrnd(5,2,sum(I1),1);
        I2 = ismember(x2,{'A'}) & x1 <= 5;
        y(I2) = wblrnd(1,5,sum(I2),1);
        I3 = ismember(x2,{'B'}) & x1 <= 5;
        y(I3) = wblrnd(.5,.9,sum(I3),1);
        I4 = ismember(x2,{'C','D'}) & x1 <= 3;
        y(I4) = wblrnd(5,5,sum(I4),1);
        I5 = ismember(x2,{'C','D'}) & x1 > 3 & x1 <= 7;
        y(I5) = wblrnd(.5,.5,sum(I5),1);
        I6 = ismember(x2,{'C','D'}) & x1 >  7;
        y(I6) = [interp1(surv1,times,rand(sum(I6),1))];
        % Do some censoring
        cens = exprnd(censparms(jj),N,1);
        ind = y <= cens;
        y(~ind) = cens(~ind);
        Y = [y,double(ind)];
        X = dat;
        disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (N))));
        DATA(ii,jj).Y = Y;
        DATA(ii,jj).X = X;
        cntr = cntr + 1;
    end
end

parpool(8);
% SETTINGS FOR THE USER TO CHANGE
% SETTINGS FOR THE USER TO CHANGE
lowrep = 1;
highrep = 10;
simnum = 9;

MCMC = [10000]; % MCMC iterations for the small and large K respectively
Kvals = [100];
nn = length(Kvals) * numel(DATA) * highrep;
SEEDS = Kvals .* simnum .* (1:nn); % Seed for K = 100
cntr = 1;
for rep = lowrep:highrep
    for kk=1:length(Kvals) % two different K values (potentially)
        for ii=1:length(ns)
            for jj=1:length(censparms)
                if (rep == 10 && ns(ii) == 1000 && jj == 3) || ...
                   (rep == 10 && ns(ii) < 1000)
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
    % You ought to check to see if there are any gaps or clusters of data 
    %   which may be fixed by adjusting the number of bins, etc.
    %load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/mcmc_id1.mat');
    %load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep9/mcmc_id1.mat');
    
    % Posterior Analysis (and figures) of the first result
    bname = '../output/sim9final_bigk_N1000_cens1_K100_rep';
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
    mtrees = mallow(output);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4];
    saveas(fig,'../figs/sim9tree.jpg')
    print('../figs/sim9tree','-depsc')

    % get_surv_tree(thetree,Y,X,1000,1,[],[])
    % Plot Survival Functions
    figure()
    subplot(2,3,1)
    get_surv_tree(thetree,Y,X,10000,1,table(6,{'A'}),[],.05,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),10000);
        yy = 1 - wblcdf(xx,5,2);
        plot(xx,yy)
        xlim([0,15])
        xlabel('t')
        ylabel('S(t)')
        title('W(5,2)')
    hold off
    
    subplot(2,3,2)
    get_surv_tree(thetree,Y,X,10000,1,table(4,{'A'}),[],.05,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),10000);
        yy = 1 - wblcdf(xx,1,5);
        plot(xx,yy)
        xlim([0,3])
        xlabel('t')
        ylabel('S(t)')
        title('W(1,5)')
    hold off
    
    subplot(2,3,3)
    get_surv_tree(thetree,Y,X,10000,1,table(4,{'B'}),[],.05,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),10000);
        yy = 1 - wblcdf(xx,.5,.9);
        plot(xx,yy)
        xlim([0,5])
        xlabel('t')
        ylabel('S(t)')
        title('W(.5,.9)')
    hold off
    
    subplot(2,3,4)
    get_surv_tree(thetree,Y,X,10000,1,table(1,{'C'}),[],.05,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),10000);
        yy = 1 - wblcdf(xx,5,5);
        plot(xx,yy)
        xlim([0,10])
        xlabel('t')
        ylabel('S(t)')
        title('W(5,5)')
    hold off
    
    subplot(2,3,5)
    get_surv_tree(thetree,Y,X,10000,1,table(5,{'C'}),[],.05,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),10000);
        yy = 1 - wblcdf(xx,.5,.5);
        plot(xx,yy)
        xlim([0,15])
        xlabel('t')
        ylabel('S(t)')
        title('W(.5,.5)')
    hold off
    
    subplot(2,3,6)
    get_surv_tree(thetree,Y,X,10000,1,table(8,{'C'}),[],.05,[])
    hold on
        plot(times,surv1)
        xlim([0,15])
        xlabel('t')
        ylabel('S(t)')
        title('G')
    hold off
    
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4];
    saveas(fig,'../figs/sim9survival.jpg')
    print('../figs/sim9survival','-depsc')
    
    % Posterior analysis for all 9 simulations
    n_num = 3;     % index on sample sizes
    cens_num = 3;  % index on censoring rates
    K = 100;
    nreps = 10;
    basename = '../output/sim9final_bigk';
    
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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % Subset of the data wth 500 points
    load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    
    % Subset of the data wth 100 points
    load('../output/sim9_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub2_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub2);
    Treeplot(thetree)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    % 9 i s 7.7085
    % 10 is "c"
    newnode = thetree.Allnodes{6};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 100;
    %s = linspace(0,1.01,K+1);
    s = [];
    eps = 1e-10;
    tau = 1;
    l = .01;
    %mu = 4.5;
    nugget= 1e-10;
    EB = 1;
    lchild = thetree.Allnodes{nodeind(thetree,newnode.Lchild)};
    rchild = thetree.Allnodes{nodeind(thetree,newnode.Rchild)};
    [marg_y,res] = get_marginal(Ystd(newnode.Xind,:),K,s,eps,tau,l,nugget,EB);
    [marg_y_lc,res_lc] = get_marginal(Ystd(lchild.Xind,:),K,s,eps,tau,l,nugget,EB);
    [marg_y_rc,res_rc] = get_marginal(Ystd(rchild.Xind,:),K,s,eps,tau,l,nugget,EB);
    marg_y
    marg_y_lc + marg_y_rc
    
    
    % Find the true marginal
    I1 = ismember(x2,{'A','B'}) & x1 > 5;
    [marg_y1,res1] = get_marginal(Ystd(I1,:),K,s,eps,tau,l,nugget,EB);
    I2 = ismember(x2,{'A'}) & x1 <= 5;
    [marg_y2,res2] = get_marginal(Ystd(I2,:),K,s,eps,tau,l,nugget,EB);
    I3 = ismember(x2,{'B'}) & x1 <= 5;
    [marg_y3,res3] = get_marginal(Ystd(I3,:),K,s,eps,tau,l,nugget,EB);
    I4 = ismember(x2,{'C','D'}) & x1 <= 3;
    [marg_y4,res4] = get_marginal(Ystd(I4,:),K,s,eps,tau,l,nugget,EB);
    I5 = ismember(x2,{'C','D'}) & x1 > 3 & x1 <= 7;
    [marg_y5,res5] = get_marginal(Ystd(I5,:),K,s,eps,tau,l,nugget,EB);
    I6 = ismember(x2,{'C','D'}) & x1 >  7;
    [marg_y6,res6] = get_marginal(Ystd(I6,:),K,s,eps,tau,l,nugget,EB);
    
    thetree.Lliketree
    marg_y1 + marg_y2 + marg_y3 + marg_y4 + marg_y5 + marg_y6
   
    
    
    
   
    
    newnode
    lchild
    rchild
    
    montecarlo_int(thetree,res,100000)
    marg_y - get_lprior(newnode.tau,newnode.l,newnode.mu)
    montecarlo_int(thetree,res_lc,100000)
    marg_y_lc - get_lprior(lchild.tau,lchild.l,lchild.mu)
    montecarlo_int(thetree,res_rc,100000)
    marg_y_rc - get_lprior(rchild.tau,rchild.l,rchild.mu)
    
    %get_surv(Y_orig,res,ndraw,graph,ystar)
    get_surv(Y,res,1000,1,[],.05)
    hold on 
        plot(times,surv)
    hold off 
    figure()
    get_surv(Y,res_lc,1000,1,[],.05)
    hold on 
        plot(times,surv)
    hold off 
    figure()
    get_surv(Y,res_rc,1000,1,[],.05)
    hold on 
        plot(times,surv)
    hold off 
    
    plot(res.f .* res.ns)
    hold on
        plot(res_lc.f .* res_lc.ns)
        plot(res_rc.f .* res_rc.ns)
    hold off
    
    
    figure()
    plot(-exp(res.f) .* (res.a + res.b))
    hold on
        plot(-exp(res_lc.f) .* (res_lc.a + res_lc.b))
        plot(-exp(res_rc.f) .* (res_rc.a + res_rc.b))
    hold off
    
    figure()
    plot(res.f .* res.ns -exp(res.f) .* (res.a + res.b))
    hold on
        plot(res_lc.f.*res_lc.ns -exp(res_lc.f) .* (res_lc.a + res_lc.b))
        plot(res_rc.f.*res_rc.ns-exp(res_rc.f) .* (res_rc.a + res_rc.b))
    hold off
    
    llike = sum(res.f .* res.ns -exp(res.f) .* (res.a + res.b));
    llike_lc = sum(res_lc.f.*res_lc.ns -exp(res_lc.f) .* (res_lc.a + res_lc.b));
    llike_rc = sum(res_rc.f.*res_rc.ns-exp(res_rc.f) .* (res_rc.a + res_rc.b));
    
    Sigma = get_sigma(res.Z,res.tau,res.l,sqrt(nugget));
    Sigma_lc = get_sigma(res_lc.Z,res_lc.tau,res_lc.l,sqrt(nugget));
    Sigma_rc = get_sigma(res_rc.Z,res_rc.tau,res_rc.l,sqrt(nugget));
    eprior = -.5 * (res.f - res.mu)' * (Sigma \ (res.f - res.mu));
    eprior_lc = -.5 * (res_lc.f - res_lc.mu)' * (Sigma_lc \ (res_lc.f - res_lc.mu));
    eprior_rc = -.5 * (res_rc.f - res_rc.mu)' * (Sigma_rc \ (res_rc.f - res_rc.mu));
    
    thedet = -.5*ldet(res.Omegainv) - .5*ldet(Sigma);
    thedet_lc = -.5*ldet(res_lc.Omegainv) - .5*ldet(Sigma_lc);
    thedet_rc = -.5*ldet(res_rc.Omegainv) - .5*ldet(Sigma_rc);
    
    plot([llike,eprior,thedet,marg_y],'b')
    hold on
        plot([llike_lc,eprior_lc,thedet_lc,marg_y_lc],'r')
        plot([llike_rc,eprior_rc,thedet_rc,marg_y_rc],'g')
        plot([llike_lc + llike_rc,...
              eprior_lc + eprior_rc,...
              thedet_lc + thedet_rc,marg_y_lc + marg_y_rc],'m')
    hold off
    
    marg_y - get_lprior(res.tau,res.l,res.mu)
    llike + eprior + thedet
    
end


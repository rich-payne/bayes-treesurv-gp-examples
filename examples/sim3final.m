addpath('/home/grad/richard/Documents/mallick/prophaz2/src');

ns = [1000,500,100]; % sample sizes
censparms = [50,25,1]; % censoring parameters
DATA = []; % struct holding all the data

rng(3511)
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
        y(I1) = exprnd(1,sum(I1),1);
        I2 = ismember(x2,{'A','B'}) & x1 <= 5;
        y(I2) = exprnd(10,sum(I2),1);
        I3 = ismember(x2,{'C','D'}) & x1 <= 3;
        y(I3) = exprnd(.5,sum(I3),1);
        I4 = ismember(x2,{'C','D'}) & x1 > 3 & x1 <= 7;
        y(I4) = exprnd(3,sum(I4),1);
        I5 = ismember(x2,{'C','D'}) & x1 >  7;
        y(I5) = exprnd(10,sum(I5),1);
        % Do some censoring
        cens = exprnd(censparms(jj),N,1);
        ind = y <= cens;
        disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (N))));
        y(~ind) = cens(~ind);
        Y = [y,double(ind)];
        [~,I] = sort(Y(:,1));
        Y = Y(I,:);
        dat = dat(I,:);
        DATA(ii,jj).Y = Y;
        DATA(ii,jj).X = dat;
        cntr = cntr + 1;
    end
end

parpool(8);
% SETTINGS FOR THE USER TO CHANGE
lowrep = 1;
highrep = 10;
simnum = 3;

MCMC = [10000]; % MCMC iterations for the small and large K respectively
Kvals = [20];
nn = length(Kvals) * numel(DATA) * highrep;
SEEDS = simnum .* (1:nn); % Seed for K = 20
cntr = 1;
for rep = lowrep:highrep
    for kk=1:length(Kvals) % two different K values (potentially)
        for ii=1:length(ns)
            for jj=1:length(censparms)
                if (rep == 8 && ns(ii) == 1000 && jj == 3) || ...
                   (rep == 8 && ns(ii) < 1000) || ...
                   rep >= 9
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
                    cntr = cntr + 1;
                end
            end
        end
    end
end



if 0
    nind = 1;
    censind = 1;
    k_ind = 1;
    rep = 1;
    worker = 1;
    disp(strcat(['Sample size: ',num2str(ns(nind)),', Censoring level: ',...
        num2str(censind),', K = ',num2str(Kvals(k_ind))]))
    fname = ['../output/sim3final_N' num2str(ns(nind)) '_cens',...
        num2str(censind) '_K' num2str(Kvals(k_ind)) '_rep' num2str(rep),...
        '/mcmc_id' num2str(worker) '.mat'];
    Y = DATA(nind,censind).Y;
    X = DATA(nind,censind).X;
    
    load(fname)
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree)
    %plot(output.treesize)
    %tabulate(output.treesize)
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,table(6,{'A'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - expcdf(xx,1);
        plot(xx,yy)
        xlim([0,10])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,table(4,{'A'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - expcdf(xx,10);
        plot(xx,yy)
        %xlim([0,10])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,table(1,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,.5);
        plot(xx,yy)
        xlim([0,5])
    hold off
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,table(5,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,3);
        plot(xx,yy)
        xlim([0,15])
        % Kaplan Meijer Curve
%         thenode = thetree.Allnodes{9};
%         [f,x,flo,fup] = ecdf(Y(thenode.Xind,1),'censoring',Y(thenode.Xind,2) == 0,'bounds','on');
%         plot(x,1-f,x,1-flo,'--',x,1-fup,'--')
    hold off
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,table(8,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,10);
        plot(xx,yy)
    hold off
    
    
    % Subset of the data wth 600 points
    load('../output/sim3_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_sub1_K20_rep1/mcmc_id1.mat');
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,Xsub1);
    Treeplot(thetree)
    
    
    
    %ind = ismember(dat{:,2},{'C','D'}) & dat{:,1} > 3.0274 & dat{:,1} <= 6.9197;
    % sum(Y(ind,:))
    thenode = thetree.Allnodes{9};
    Ystd = Y;
    Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
    K = 20;
    eps=1e-10;
    tau=1e-10;
    l = thenode.l;
    tau = thenode.tau;
    nugget = 1e-10;
    EB = 0;
    %(Y,K,s,eps,tau,l,nugget,EB)
    [marg_y,res] = get_marginal(Ystd(thenode.Xind,:),K,[],eps,tau,l,nugget,EB);
    
    
    figure()
    get_surv_tree(thetree,Y,X,1000,1,table(8,{'C'}),[])
    hold on
        xx = linspace(.01,max(Y(:,1)),1000);
        yy = 1 - expcdf(xx,10);
        plot(xx,yy)
        % xlim([0,15])
    hold off
    
    

    get_surv_tree(thetree,Y,X,1000,1,[],[])

    tind = output.treesize == 2;
    [~,I] = max(output.llike(tind));
    thetree = output.Trees(tind);
    thetree = thetree{I};
    Treeplot(thetree);
    
    
    plot(output.As)
    plot(output.Omegas)

    [GAM,grid,thetas] = plot_surv(Y,X,thetree,10,50,[],[],[],100);

    ii = 5;
    cumhaz = squeeze(GAM(ii,:,:));
    qtiles = quantile(exp(-cumhaz),[.025,.975]);
    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,qtiles(1,:))
    plot(grid,qtiles(2,:));
    fplot(@(x) exp(-x/10),[0,60])
    hold off
end

addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
dat = readtable('../data/LUAD_RPPA_CLEAN.csv','Delimiter',',','ReadVariableNames',true);

% Important proteins: PCADHERIN, FOXO3A_pS318S321, DIRAS3, RAD51, GAB2, SF2, CRAF_pS338, MAPK_pT202Y204
y = [dat.new_death, dat.death_event];
X = dat(:,3:end);

parpool(8);
% Quick test case
% Tree_Surv(y,X,'nmcmc',10^2,'burn',10^2,...
%     'saveall',1,'filepath','../output/LUAD_RPPA_rep1/','seed',3003,...
%     'nprint',10);
bname = '../output/LUAD_RPPA_rep';
for ii=1:10
    fname = strcat([bname, num2str(ii),'/']);
    if ii > 1
        fname_resume = strcat([bname, num2str(ii-1),'/']);
    else
        fname_resume = [];
    end
    Tree_Surv(y,X,'nmcmc',10^4,'burn',0,'filepath',fname,'seed',90610 + ii,...
        'saveall',1,'nprint',1000,'bigK',100,'swapfreq',10,...
        'resume',fname_resume);
end

if 0    
    load('../output/LUAD_RPPA_rep1/mcmc_id1.mat')
    output1 = output;
    load('../output/LUAD_RPPA_rep2/mcmc_id1.mat')
    output2 = output;
    load('../output/LUAD_RPPA_rep3/mcmc_id1.mat')
    output3 = output;
    load('../output/LUAD_RPPA_rep4/mcmc_id1.mat')
    output4 = output;
    load('../output/LUAD_RPPA_rep5/mcmc_id1.mat')
    output5 = output;
    load('../output/LUAD_RPPA_rep6/mcmc_id1.mat')
    output6 = output;
    load('../output/LUAD_RPPA_rep7/mcmc_id1.mat')
    output7 = output;
    load('../output/LUAD_RPPA_rep8/mcmc_id1.mat')
    output8 = output;
    load('../output/LUAD_RPPA_rep9/mcmc_id1.mat')
    output9 = output;
    load('../output/LUAD_RPPA_rep10/mcmc_id1.mat')
    output10 = output;
    
    
    [mtree,llike,outind] = mallow(output1,output2,output3,output4,output5,output6,output7,output8,output9,output10);
    mind = 5;
    thetree = mtree{mind};
    thetree = fatten_tree(thetree,X);
    %X_new = X; % for changing names
    %X_new.Properties.VariableNames = {'X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9'};
    thetree.Varnames = {'X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8','X_9'};
    %Treeplot(thetree)
    %figure()
    %get_surv_tree(thetree,y,X,1000,1,[],[],.05,[])
    Treeplot(thetree,y,X)
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 15 8];
    thetexts = findall(fig,'-property','FontSize');
    axind = [];
    texind = [];
    for ii=1:length(thetexts)
        if strcmp(class(thetexts(ii)),'matlab.graphics.axis.Axes')
            axind = [axind,ii];
        else
            texind = [texind,ii];
        end
    end       
    set(thetexts(texind),'FontSize',25)
    saveas(fig,'../figs/LUADtree.jpg')
    print('../figs/LUADtree','-depsc')
    
    
    % Plot with Kaplan Meier curves
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
    ncol = 3; 
    xlimmax = [15, 4, 20, 7, 50, 12, 10,50,50,50,50];
    figure()
    for ii=1:ntermnodes
        subplot(nrow,ncol,ii)
        get_surv_tree(thetree,y,X,10000,1,Xplot(ii,:),[],.05,[])
        title('');
        % Plot the median survival function of the observed data
        node = thetree.Allnodes{uinds(ii)};
        % Plot Kaplan Meier Curve
        hold on
            ypart = y(node.Xind,1);
            cens = y(node.Xind,2) == 0;
            [f,x,flo,fup] = ecdf(ypart,'censoring',cens);
            plot(x,1-f,':r')
            plot(x,1-flo,'--r');
            plot(x,1-fup,'--r');         
        hold off
    end
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 6 4];
    saveas(fig,'../figs/LUADsurvivaldetail.jpg')
    print('../figs/LUADsurvivaldetail','-depsc')
    
    
    
    % Analysis of weird curve (both 8 and 9 are a problem)
    node = thetree.Allnodes{9};
    y(node.Xind,:)
    
    
    % Further analysis
    termind = termnodes(thetree);
    for ii=1:length(termind)
        thenode = thetree.Allnodes{termind(ii)};
        ytmp = y(thenode.Xind,:);
        mean(ytmp(:,2))
    end
    % For the one with 100% censoring
    ii = 1;
    thenode = thetree.Allnodes{termind(ii)};
    nind = nodeind(thetree,thenode.Parent);
    pnode = thetree.Allnodes{nind}
    
    
    plot(output.llike)
%     [~,I] = max(output.llike);
%     thetree = output.Trees{I};
    Treeplot(thetree)
    %plot(output.treesize)
    tabulate(output.treesize)

    figure()
    get_surv_tree(thetree,y,X,1000,1,[],[])
    
    
    thenode = thetree.Allnodes{10};
    thenode.Updatellike = 0;
    y(node.Xind,:)
    Ystd = y;
    Ystd(:,1) = Ystd(:,1) / max(Ystd(:,1));
    delta_z = 1;
    delta_pi = 2.5;
    [fullmarg,ZZ] = get_fullmarginal(Ystd,thenode,thetree,delta_z,delta_pi);
    
    scatter3(ZZ(:,1),ZZ(:,2),ZZ(:,3))
    
    
    % Diagnose a node...
    node = thetree.Allnodes{9};
    Ystd = y;
    Ystd(:,1) = Ystd(:,1) / max(Ystd(:,1));
    K = 100;
    tau = 5;
    l = .000001;
    EB = 1;
    % s,eps,tau,l,nugget,EB)
    [marg_y,res] = get_marginal(Ystd(node.Xind,:),K,[],thetree.eps,tau,l,0,EB);
    llike = -sum(exp(res.f) .* (res.a + res.b)); % likelihood
    lprior = res.marg_y - llike;
    
    
    % Figure out prior and likelihood contributions on a grid
    taus = linspace(.00001,.01,10);
    ls = linspace(.00001,.01,10);
    [taugrid,lgrid] = meshgrid(taus,ls);
    node = thetree.Allnodes{8};
    Ystd = y;
    Ystd(:,1) = Ystd(:,1) / max(Ystd(:,1));
    K = 100;
    EB = 0;
    LLIKE = zeros(length(ls),length(taus));
    LPRIOR = LLIKE;
    for ii=1:length(taus)
        for jj=1:length(ls)
            [~,res] = get_marginal(Ystd(node.Xind,:),K,[],thetree.eps,taugrid(ii,jj),lgrid(ii,jj),0,EB);
            llike = -sum(exp(res.f) .* (res.a + res.b)); % likelihood
            lprior = res.marg_y - llike;
            LLIKE(ii,jj) = llike;
            LPRIOR(ii,jj) = lprior;
        end
    end

    subplot(1,2,1)
    contour(taugrid,lgrid,LLIKE)
    xlabel('tau');
    ylabel('l');
    subplot(1,2,2)
    contour(taugrid,lgrid,LPRIOR)
    xlabel('tau');
    ylabel('l');
    
	% Check out determinant of the product to see how much it changes.
    -.5*ldet(Omegainv,'chol')
    .5*ldet(Sigmainv,'chol')
    Sigmainv = 
    -.5 * ldet(res.Omegainv
    
    
        
    plot(res.f)
    subplot(2,1,1)
    plot(exp(res.f))
    subplot(2,1,2)
    plot(res.a)
    hold on 
        plot(res.b)
        plot(res.a + res.b,'--')
    hold off
    llike = sum(exp(res.f) .* (res.a + res.b)); % likelihood
    lprior = res.marg_y - llike;
    
    % Look at marginal
    tau = .00001;
    l = .000001;
    dtau = .0001;
    dl = .000001;
    ntau = 50;
    nl = 50;
    [truemarginal,marges,XX,YY] = get_truemarginal(Ystd(node.Xind,:),tau,l,dtau,dl,ntau,nl,K);
    figure()
    s = surf(XX,YY,marges)
    s.EdgeColor = 'none';
    xlabel('tau')
    ylabel('l')
    zlabel('marginal')
    
    
end
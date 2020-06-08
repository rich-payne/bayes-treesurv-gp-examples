addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
dat = readtable('../data/KIRC.csv','Delimiter',',','ReadVariableNames',true);

% Important proteins: PCADHERIN, FOXO3A_pS318S321, DIRAS3, RAD51, GAB2, SF2, CRAF_pS338, MAPK_pT202Y204
y = [dat.days, dat.event];
X = array2table([dat.PCADHERIN, dat.FOXO3A_pS318S321, dat.DIRAS3,...
    dat.RAD51, dat.GAB2, dat.SF2, dat.CRAF_pS338, dat.MAPK_pT202Y204]);

parpool(8);
% Tree_Surv(y,X,'nmcmc',10^3,'burn',0,'filepath','../output/KIRC_rep1/','seed',632,...
%     'saveall',1,'nprint',100,'bigK',100,'swapfreq',10);
% Tree_Surv(y,X,'nmcmc',10^4,'burn',0,'filepath','../output/KIRC_rep2/','seed',236,...
%     'saveall',1,'nprint',1000,'bigK',100,'swapfreq',10,...
%     'resume','../output/KIRC_rep1/');
bname = '../output/KIRC_rep';
for ii=3:11
    fname = strcat([bname, num2str(ii),'/']);
    fname_resume = strcat([bname, num2str(ii-1),'/']);
    Tree_Surv(y,X,'nmcmc',10^4,'burn',0,'filepath',fname,'seed',111 + ii,...
        'saveall',1,'nprint',1000,'bigK',100,'swapfreq',10,...
        'resume',fname_resume);
end


if 0 
    load('../output/KIRC_rep2/mcmc_id1.mat')
    trees = mallow(output);
    tind = 9; % chosen by user based on plot
    thetree = trees{tind};
    thetree = fatten_tree(thetree,X);
    
    plot(output.llike + output.llike)
    %[~,I] = max(output.llike + output.lprior);
    %thetree = output.Trees{I};
    Treeplot(thetree)
    plot(output.treesize)
    tabulate(output.treesize)
    
    plot(output.treesize,output.lprior,'o')

    get_surv_tree(thetree,y,X,1000,1,[],[])
    
    
    % Diagnose a node...
    node = thetree.Allnodes{21};
    Ystd = y;
    Ystd(:,1) = Ystd(:,1) / max(Ystd(:,1));
    K = 1000;
    tau = 1;
    l = .1;
    EB = 1;
    % s,eps,tau,l,nugget,EB)
    [marg_y,res] = get_marginal(Ystd(node.Xind,:),K,[],thetree.eps,tau,l,0,EB);
    
    % get Kaplan Meier curve
    t = y(node.Xind,1);
    censored = y(node.Xind,2) == 0;
    [f,x,flo,fup] = ecdf(t,'censoring',censored,'bounds','on');
    get_surv(y,res,1000,1,[],.05)
    hold on
        plot(x,1-f);
        plot(x,1-flo);
        plot(x,1-fup);
    hold off
    
    
end
addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
dat = readtable('../data/cancer.csv','Delimiter',',','ReadVariableNames',true);

y = [dat.time, dat.status - 1];
X = dat(:,4:size(dat,2)); % Ignoring institution

parpool(8);
% Tree_Surv(y,X,'nmcmc',10^2,'burn',10^2,...
%     'saveall',1,'filepath','../output/cancer_rep1/','seed',3003,...
%     'nprint',10000);
bname = '../output/cancer_rep';
for ii=1:10
    fname = strcat([bname, num2str(ii),'/']);
    if ii > 1
        fname_resume = strcat([bname, num2str(ii-1),'/']);
    else
        fname_resume = [];
    end
    Tree_Surv(y,X,'nmcmc',10^4,'burn',0,'filepath',fname,'seed',5200 + ii,...
        'saveall',1,'nprint',1000,'bigK',100,'swapfreq',10,...
        'resume',fname_resume);
end

if 0 
    load('../output/cancer_rep1/mcmc_id1.mat')
    output1 = output;
    load('../output/cancer_rep2/mcmc_id1.mat')
    output2 = output;
    load('../output/cancer_rep3/mcmc_id1.mat')
    output3 = output;
    load('../output/cancer_rep4/mcmc_id1.mat')
    output4 = output;
    load('../output/cancer_rep5/mcmc_id1.mat')
    output5 = output;
    load('../output/cancer_rep6/mcmc_id1.mat')
    output6 = output;
    load('../output/cancer_rep7/mcmc_id1.mat')
    output7 = output;
    load('../output/cancer_rep8/mcmc_id1.mat')
    output8 = output;
    load('../output/cancer_rep9/mcmc_id1.mat')
    output9 = output;
    load('../output/cancer_rep10/mcmc_id1.mat')
    output10 = output;
    [mtree,llike,outind] = mallow(output1,output2,output3,output4,output5,output6,output7,output8,output9,output10);
    mind = 6;
    thetree = mtree{mind};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree,y,X)
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
    node = thetree.Allnodes{8};
    Ystd = y;
    Ystd(:,1) = Ystd(:,1) / max(Ystd(:,1));
    K = 100;
    tau = 5;
    l = .000001;
    EB = 0;
    % s,eps,tau,l,nugget,EB)
    [marg_y,res] = get_marginal(Ystd(node.Xind,:),K,[],thetree.eps,tau,l,0,EB);
    marg_y
    node.Llike
    res.l
    res.tau
    node.l
    node.tau
    
    plot(res.f)
    subplot(2,1,1)
    plot(exp(res.f))
    subplot(2,1,2)
    plot(res.a)
    hold on 
        plot(res.b)
        plot(res.a + res.b,'--')
    hold off
    sum(exp(res.f) .* (res.a + res.b)) % likelihood
    
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
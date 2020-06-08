% There is too much censoring in this dataset

addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
dat = readtable('../data/KICH.csv','Delimiter',',','ReadVariableNames',true);

% Important proteins: PCADHERIN, FOXO3A_pS318S321, DIRAS3, RAD51, GAB2, SF2, CRAF_pS338, MAPK_pT202Y204
y = [dat.days, dat.event];
X = array2table([dat.PCADHERIN, dat.FOXO3A_pS318S321, dat.DIRAS3,...
    dat.RAD51, dat.GAB2, dat.SF2, dat.CRAF_pS338, dat.MAPK_pT202Y204]);

parpool(8);
Tree_Surv(y,X,'nmcmc',10^3,'burn',0,'filepath','../output/KICH_rep1/','seed',103325,...
    'saveall',1,'leafmin',10,'bigK',100,'nprint',100,'swapfreq',10);

if 0 
    load('../output/KICH_rep1/mcmc_id1.mat')
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike);
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
    Treeplot(thetree)
    tabulate(output.treesize)
    
    plot(output.treesize,output.llike + output.lprior,'o')

    get_surv_tree(thetree,y,X,1000,1,[],[])
    
    % Find all unique trees based on log likelihood
    [~,IU] = unique(output.llike); % The uniqueness of trees 
    model_probs = exp(output.llike(IU))/sum(exp(output.llike(IU)));
    nmodels = 10;
    themods = randsample(IU,nmodels,true,model_probs);
    for ii=1:nmodels
        get_surv
    
end
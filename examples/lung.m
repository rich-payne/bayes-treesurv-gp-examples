addpath('/home/grad/richard/Documents/mallick/prophaz2');
dat = readtable('../data/lung.txt','Delimiter',' ','ReadVariableNames',true);
y = [dat.time, dat.status - 1];
X = array2table([dat.age, dat.sex, dat.ph_ecog]);

parpool(8);
TreePH_MCMCparalleltemp(y,X,'nmcmc',1000,'burn',1000,'filepath','../output/lung/','seed',3003);

load('../output/lung/mcmc_id1.mat')
plot(output.llike)
[~,I] = max(output.llike);
thetree = output.Trees{I};
Treeplot(thetree)
plot(output.treesize)
tabulate(output.treesize)

get_surv_tree(thetree,y,X,10,1)
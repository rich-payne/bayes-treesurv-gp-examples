% Generate a single covariate on which to partition
rng(424)
n = 1000;
x1 = rand(n,1);
X = table(x1);

% Create a non-parametric (arbitrary) survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
% Plot survival function
plot(times,surv);
xlabel('t');
ylabel('S(t)')
title('Survival Function')

% Generate data from two partitions
ind1 = X{:,1} < .5;
ind2 = ~ind1;
y = zeros(n,1);
y(ind1) = exprnd(1,sum(ind1),1);
y(ind2) = interp1(surv,times,rand(sum(ind2),1));

% Induce some censoring
cens = exprnd(50,n,1);
ind_cens = y <= cens;
delta = double(ind_cens);
Y = [y, delta];

MCMC = 100; % MCMC iterations
burn = 0; % burn in
Kval = 100; % Number of bins for the hazard function
fpath = './output/'; % directory to save results to
saveall = 1; % save all the output from each tempered chain
nprint = 10; % print MCMC info every 10 iterations
swapfreq = 10; % propose a swap of parallel chains every 10 iterations

% Start a parallel pool with 4 cores
parpool(4);

% Run MCMC
Tree_Surv(Y,X,'nmcmc',MCMC,'burn',burn,'filepath',fpath,...
    'seed',1990,'bigK',Kval,'saveall',saveall,'swapfreq',swapfreq,...
    'nprint',nprint);

% Load data from untempered chain
load([fpath,'mcmc_id1.mat']);
% Find tree with highest likelihood
[~,I] = max(output.llike);
% plot tree with highest likelihood (along with posterior survival curves)
Treeplot(output.Trees{I}) % Simple plot of tree
Treeplot(output.Trees{I},Y,X,0); % Tree with survival curves at terminal nodes
                                 % 0 indicates no Kaplan Meier curves overlayed
mallow(output); % Mallow plot indicates 2 partition elements is best


% Continue the MCMC chain using the 'resume' feature
fpath2 = './output2/';
Tree_Surv(Y,X,'nmcmc',MCMC,'burn',burn,'filepath',fpath2,...
    'seed',1991,'bigK',Kval,'saveall',saveall,'swapfreq',swapfreq,...
    'nprint',nprint,'resume',fpath);
output1 = output; % store output from first run;
% Load data from untempered chain (from continued chain)
load([fpath2,'mcmc_id1.mat']);
output2 = output; % store output from second (continued) run;

% Look at second half of output;
[~,I2] = max(output2.llike); % get highest likelihood 
Treeplot(output2.Trees{I2},Y,X,1); % 1 indicates Kaplan Meier curves overlayed
mallow(output1,output2); % Mallow plot indicates 2 partition elements is best from both chains
  
addpath(genpath('../bayes-treesurv-gp/'));

% Generate some data from exponential proportional hazards model
n = 1000; % sample sizes
censparms = 50; % censoring parameters
rng(351199)
cntr = 1; 
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
cens = exprnd(censparms,n,1);
ind = y <= cens;
disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (n))));
y(~ind) = cens(~ind);
ds = double(ind);
Y = [y,ds];
X = array2table(X);

% SETTINGS FOR THE USER TO CHANGE
lowrep = 1;
highrep = 1;
simnum = 7;

MCMC = 10000; % MCMC iterations
Kvals = 20;
SEEDS = simnum .* (1:highrep);
cntr = 1;
for rep = lowrep:highrep
    basename = strcat(['../output/sim',num2str(simnum),'final_N',num2str(n),'_rep']);
    fname = strcat([basename,num2str(rep),'/']);
    % Run the code
    if rep > 1
        fname_resume = strcat([basename,num2str(rep-1),'/']);
        Tree_Surv(Y,X,'nmcmc',MCMC,'burn',0,'filepath',fname,...
            'seed',SEEDS(cntr),'bigK',Kvals,'nprint',1000,'saveall',1,...
            'swapfreq',10,'resume',fname_resume);
    else
        Tree_Surv(Y,X,'nmcmc',MCMC,'burn',0,'filepath',fname,...
            'seed',SEEDS(cntr),'bigK',Kvals,'nprint',1000,'saveall',1,...
            'swapfreq',10);
    end
    cntr = cntr + 1;
end

load('../output/sim7final_N1000_rep1/mcmc_id1.mat')
% plot(output.llike + output.lprior)
[~,I] = max(output.llike + output.lprior);
thetree = output.Trees{I};
Treeplot(thetree)
ylim([-3, 0])
thetree = fatten_tree(thetree, X);
tnodes = termnodes(thetree);

xmax = [5, 3, 12, 21, 9, 60, 15, 21];
for ii = 1:length(tnodes)
    subplot(2, 4, ii);
    ind = thetree.Allnodes{tnodes(ii)}.Xind;
    Xsub = X(ind, :);
    Ysub = Y(ind, :);
    x0 = Xsub(1, :);
    true_hazards = 1 / beta * exp(table2array(Xsub) * betas);
    med_hazard = median(true_hazards);
    %xgrid = 0:.01:max(Ysub(:, 1));
    xgrid = 0:.01:60;
    the_surv = 1 - expcdf(xgrid, 1 / med_hazard);
    get_surv_tree(thetree,Y,X,10000,1,x0,[],.05,' ');
    xlim([0, xmax(ii)]);
    xlabel('time');
    if ii == 1 || ii == 5
        ylabel('S(t)');
    end
    hold on;
        plot(xgrid, the_surv, 'k')
    hold off;
end

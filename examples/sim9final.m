% Same as sim8 but we have changed one of the Weibull distributions

addpath(genpath('../bayes-treesurv-gp/'));

% Create my own survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
%plot(times,surv)
ht = -log(surv);
% Data from first function
theta1 = 1;
surv1 = exp(-ht*theta1);

% CART example (Figure 2)
N = 1000; % sample sizes
censparm = 50; % censoring parameters

rng(3199771)
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
cens = exprnd(censparm,N,1);
ind = y <= cens;
y(~ind) = cens(~ind);
Y = [y,double(ind)];
X = dat;
disp(strcat('Censoring rate: ',num2str(1 - sum(ind) / (N))));

% SETTINGS FOR THE USER TO CHANGE
lowrep = 1;
highrep = 1;
simnum = 9;

MCMC = 10000; % MCMC iterations for the small and large K respectively
Kvals = 20;
SEEDS = simnum .* (1:highrep); % Seed for K = 20
cntr = 1;
for rep = lowrep:highrep
    basename = strcat(['../output/sim',num2str(simnum),'final_N',num2str(N),'_rep']);
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


% Posterior Analysis (and figures)
bname = '../output/sim9final_N1000_rep';
nrep = 1;
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
thetree = output.Trees{I};
thetree = fatten_tree(thetree,X);
Treeplot(thetree);

thetree.K = 100; % smooth things out for plots

figure()
subplot(2,3,1)
get_surv_tree(thetree,Y,X,1000,1,table(6,{'A'}),[],.05,'W(5, 2)')
hold on
    xx = linspace(.01,max(Y(:,1)),100);
    yy = 1 - wblcdf(xx,5,2);
    plot(xx,yy,'k')
    xlim([0,15])
    ylabel('S(t)');
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
hold off

subplot(2,3,2)
get_surv_tree(thetree,Y,X,1000,1,table(4,{'A'}),[],.05,'W(1, 5)')
hold on
    xx = linspace(.01,max(Y(:,1)),100);
    yy = 1 - wblcdf(xx,1,5);
    plot(xx,yy,'k')
    xlim([0,3])
hold off

subplot(2,3,3)
get_surv_tree(thetree,Y,X,1000,1,table(4,{'B'}),[],.05,'W(.5, .9)')
hold on
    xx = linspace(.01,max(Y(:,1)),100);
    yy = 1 - wblcdf(xx,.5,.9);
    plot(xx,yy,'k')
    xlim([0,5])
hold off

subplot(2,3,4)
get_surv_tree(thetree,Y,X,1000,1,table(1,{'C'}),[],.05,'W(5, 5)')
hold on
    xx = linspace(.01,max(Y(:,1)),1000);
    yy = 1 - wblcdf(xx,5,5);
    plot(xx,yy,'k')
    xlim([0,10])
    ylabel('S(t)');
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right')
hold off

subplot(2,3,5)
get_surv_tree(thetree,Y,X,1000,1,table(5,{'C'}),[],.05,'W(.5, .5)')
hold on
    xx = linspace(.01,max(Y(:,1)),1000);
    yy = 1 - wblcdf(xx,.5,.5);
    plot(xx,yy,'k')
    xlim([0,15])
    xlabel('t');
hold off

subplot(2,3,6)
get_surv_tree(thetree,Y,X,1000,1,table(8,{'C'}),[],.05,'G')
hold on
    plot(times,surv1,'k')
    xlim([0,15])
hold off


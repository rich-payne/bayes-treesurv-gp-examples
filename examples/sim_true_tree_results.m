addpath(genpath('../bayes-treesurv-gp/'));

% Posterior Analysis (and figures)
bname = '../output/sim_true_tree_rep';
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
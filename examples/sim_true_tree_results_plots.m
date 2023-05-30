% NOT FINISHED...


x1 = unifrnd(0,10,N,1);
grps = {'A','B','C','D'};
for 
pdraws = get_surv_tree(maxtree{sim}, Y, X, 10000, 0, X(sim, :), [], .05, []);


% use get_surv_tree to get survival function for a specific covariate to
% calculate mean, MSE, bias, etc.
output = alloutput{maxrep(sim)};
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
% Compares two exponential survival functions

addpath('/home/grad/richard/Documents/mallick/prophaz2/src');
% Simulate data
rng(33330424);
n = 100;
%n = 1000;
x = rand(n,1);
ind = x < .5;
y = zeros(n,1);
y(ind) = gamrnd(1,1/2,sum(ind),1);
y(~ind) = gamrnd(1,1/10,sum(~ind),1);
%y(~ind) = gamrnd(1,100,sum(~ind),1);
% Do some censoring
% cens = gamrnd(1,1,n,1);
cens = 1000 * ones(n,1);
ind = y < cens;
ds = double(ind); % 1 if not censored, 0 if censored;
y(~ind) = cens(~ind);

X = table(x);
[~,I] = sort(y);
y = y(I);
ds = ds(I);
X = X(I,:);
Y = [y, ds];
% 
Ystd = Y;
Ystd(:,1) = Ystd(:,1) ./ max(Ystd(:,1));
nugget = 0;
l = .01;
tau =.1;
EB = 1;
K = 20;
tic 
%for ii=1:100
%    ii
    [marg_y,res] = get_marginal(Ystd,K,[],1e-9,tau,l,nugget,EB);
%end
toc
% get_f(res.ns,res.a,res.b,res.Z,res.tau,res.l,nugget,1e-10)
% 
% 
% 
[marg_y1,res1] = get_marginal(Ystd(X{:,1} < .5,:),K,[],1e-9,tau,l,1e-9,1);
[marg_y2,res2] = get_marginal(Ystd(X{:,1} >= .5,:),K,[],1e-9,tau,l,1e-9,1);
marg_y
marg_y1 + marg_y2
figure()
out = get_surv(Y,res1,1000,1,[],.05);
newgrid = out.ystar * max(Y(:,1));
hold on
    plot(newgrid,1 - expcdf(out.ystar*max(Y(:,1)),1/2),'-r')
hold off

figure()
out = get_surv(Y,res2,1000,1,[],.05);
newgrid = out.ystar * max(Y(:,1));
hold on
    plot(newgrid,1 - expcdf(newgrid,100),'-r')
hold off

% [tau,l]= get_thetas_EM(Y,1,.1,1e-6,20,10,1e-3);





parpool(8);
Tree_Surv(Y,X,'nmcmc',10,'burn',0,...
    'filepath','../output/sim1_fastprofile_rep1/','seed',1,'saveall',1,'bigK',20,...
    'swapfreq',10,'nprint',10,'parallelprofile',1);


% Testing sequential save
Tree_Surv(Y,X,'nmcmc',10,'burn',0,...
    'filepath','../output/sim1_seqsave_rep1/','seed',1,'saveall',1,'bigK',20,...
    'swapfreq',10,'nprint',10);

% Verify saving worked just fine
for ii=1:8
    load(['../output/sim1_seqsave_rep1/mcmc_id' num2str(ii) '.mat'])
    figure()
    plot(output.llike)
end

Tree_Surv(Y,X,'nmcmc',10,'burn',0,...
    'filepath','../output/sim1_seqsave_rep2/','seed',2,'saveall',1,'bigK',20,...
    'swapfreq',10,'nprint',10,'resume','../output/sim1_seqsave_rep1/');

Tree_Surv(Y,X,'nmcmc',10,'burn',0,...
    'filepath','../output/sim1_seqsave_rep3/','seed',2,'saveall',0,'bigK',20,...
    'swapfreq',10,'nprint',10,'resume','../output/sim1_seqsave_rep2/');



% Tree_Surv(Y,X,'nmcmc',100,'burn',100,...
%     'filepath','../output/sim1_K20_nugget/','seed',1,'saveall',1,'bigK',20);
% Tree_Surv(Y,X,'nmcmc',100,'burn',100,...
%     'filepath','../output/sim1_K100_nugget/','seed',512,'saveall',1,'bigK',100,.../
%     'nprint',5);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',10000,...
%     'filepath','../output/sim1_rep1/','seed',5123,'saveall',1,'bigK',20,.../
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_rep2/','seed',1325,'saveall',1,'bigK',20,...
%     'nprint',1000,...
%     'resume','../output/sim1_rep1/');
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_1overlprior_rep3ish/','seed',9876,'saveall',1,'bigK',20,...
%     'nprint',1000,...
%     'resume','../output/sim1_rep2/');

% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_intercept_rep1/','seed',1,'saveall',1);

% Full prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_fullprior_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);

% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim1_fullprior_rep2/','seed',6332,'saveall',1,...
%     'nprint',1000,'resume','../output/sim1_fullprior_rep1/');

% Full prior, mu marginalized
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',80000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_rep2/','seed',801,'saveall',1,...
%     'nprint',1000,'resume','../output/sim1_fullprior_mumarg_rep1/');

% Reversed prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_revprior_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_revprior_bigvartau_rep1/','seed',1,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_revprior_bigvartau_tl_rep1/','seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative
% Tree_Surv(Y,X,'nmcmc',20000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_revprior_bigvartau_tl_fixed_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);

% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid, K=20
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_K20_rep1/',...
%     'seed',12,'saveall',1,'bigK',20,...
%     'nprint',1000);


% Reversed prior, larger var for tau, l has t-prior, fixed derivative,
% fixed grid
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/',...
%     'seed',12,'saveall',1,...
%     'nprint',1000);
% Tree_Surv(Y,X,'nmcmc',10000,'burn',0,...
%     'filepath','../output/sim1_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep2/',...
%     'seed',12113,'saveall',1,...
%     'nprint',1000,'resume','../output/sim1_fullprior_mumarg_revprior_bigvartau_tl_fixed_grid_rep1/');


if 0

    load('../output/sim1_fast_rep1/mcmc_id1.mat')
    plot(output.llike + output.lprior)
    [~,I] = max(output.llike + output.lprior);
    thetree = output.Trees{I};
    thetree = fatten_tree(thetree,X);
%     Ystd = Y;
%     Ystd(:,1) = Y(:,1) ./ max(Y(:,1));
%     thetree = llike_termnodes(thetree,Ystd);
    Treeplot(thetree)

    figure()
    get_surv_tree(thetree,Y,X,1000,1,.3,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - gamcdf(xx,1,.5);
        plot(xx,yy)
    hold off

    figure()
    get_surv_tree(thetree,Y,X,1000,1,.7,[])
    hold on
        xx = linspace(.01,max(Y(:,1)),100);
        yy = 1 - gamcdf(xx,1,1/10);
        plot(xx,yy)
    hold off



    plot(output.As)
    plot(output.Omegas)
    plot(output.As,output.Omegas,'o')    


    % Get posterior of the survival function...
    % Suppose we know omega, a,theta, and the true partition...
    omega = thetree.omega;
    a = thetree.a;
    thetas = [2.85,.2802];
    prt = double(X{:,1} > .5) + 1; % The final partition
    nsamp = 10000; % number of Monte Carlo samples at each tt

    A = zeros(size(Y,1),1);
    for jj=1:size(Y,1)
        tmp = 0;
        for kk=1:length(thetas)
            tmp = tmp + sum(Y(prt==kk,1) >= Y(jj,1)) * thetas(kk);
        end
        A(jj) = tmp;
    end




    % Estimate the survival function at the vector tt
    tt = linspace(.01,2,25);
    GAM = zeros(nsamp,length(tt));
    for ii=1:length(tt)
        ii
        thegam = zeros(nsamp,1);
        U = zeros(nsamp,1);
        jj = 1;
        while tt(ii) >= Y(jj,1) && jj < size(Y,1);
    %         A = 0;
    %         A2 = 0;
    %         for kk=1:length(thetas)
    %             A = A + sum(Y(prt==kk,1) >= Y(jj,1)) * thetas(kk);
    %             A2 = A2 + sum(Y(prt==kk,1) >= Y(jj+1,1)) * thetas(kk);
    %         end
            if jj == 1
                thediff = Y(jj,1);
            else
                thediff = Y(jj,1) - Y(jj-1,1);
            end
            thegam = thegam + gamrnd(a*omega*thediff,1/(a + A(jj)),nsamp,1);
            if Y(jj,2) % only do for non-censored observations
                %U = U + slicesampler(nsamp,a + A(jj),a + A(jj+1),.5,10);
                U = U + rejectsamp(nsamp,a + A(jj),a + A(jj+1));
                %U = U + slicesample(.5,nsamp,'pdf',@(x) fu(x,a + A(jj),a + A(jj+1),0));
            end
            jj = jj + 1;
        end
    %     for kk=1:length(thetas)
    %         A = A + sum(Y(prt==kk,1) >= Y(jj,1)) * thetas(kk);
    %     end
        if jj == 1
            thediff = tt(ii);
        else
            thediff = tt(ii) - Y(jj-1,1);
        end
        delta = gamrnd(a*omega*thediff,1/(a + A(jj)),nsamp,1);
        GAM(:,ii) = thegam + U + delta;
    end


    [GAM,grid,thetas] = plot_surv(Y,X,thetree,10,50,[],[],[],100);

    ii = 2;
    cumhaz = squeeze(GAM(ii,:,:));
    plot(grid,mean(exp(-cumhaz)))
    hold on
    plot(grid,min(exp(-cumhaz)))
    plot(grid,max(exp(-cumhaz)))
    fplot(@(x) exp(-x*10),[0,2])
    hold off


    ii = 2;
    cumhaz = GAM;
    plot(tt,mean(exp(-cumhaz*thetas(ii))))
    hold on
    plot(tt,min(exp(-cumhaz*thetas(ii))))
    plot(tt,max(exp(-cumhaz*thetas(ii))))
    fplot(@(x) exp(-x*1),[0,2])
    hold off








    ypart = Y(X{:,1} < .5,:);
    fplot(@(x) ptheta_t(x,ypart,thetree.a,thetree.omega,thetree.theta_shape,thetree.theta_rate,0),[0,15])
    ypart = Y(X{:,1} >= .5,:);
    fplot(@(x) ptheta_t(x,ypart,thetree.a,thetree.omega,thetree.theta_shape,thetree.theta_rate,0),[0,5])

    thetapart = @(x) ptheta_t(x,ypart,thetree.a,thetree.a0,thetree.b0,thetree.theta_shape,thetree.theta_rate,0);
    thing = slicesample(1,10000,'pdf',thetapart);
    histogram(thing)
end






            
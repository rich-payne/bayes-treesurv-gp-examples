% seeing the effect of the K on various survival functions and the
% likelihood values

% Exponential data
% Simulate data
rng(333);
n = 100;
y = exprnd(1/2,n,1);
cens = exprnd(5,n,1);
ind = y < cens;
ds = double(ind); % 1 if not censored, 0 if censored;
y(~ind) = cens(~ind);
Y = [y, ds];

% Obtain the marginal of Y for various number of K
K_max = 20; % Maximum K
Ystd = Y;
Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
eps =1e-10;
l = .01;
tau = .01;
nugget = 1e-10;
EB = 1;
taus = [.01,.5,1,5];
ls = [.01,.5,.1,5];
%(Y,K,s,eps,tau,l,nugget,EB)
MARG = zeros(K_max,1);
for ii = 1:K_max
    K = ii
    for jj=1:length(taus)
        for kk=1:length(ls)
            [marg_y_new,~] = get_marginal(Ystd,K,[],eps,taus(jj),ls(kk),nugget,EB);
            marg_y = max([marg_y, marg_y_new]);
        end
    end
    %MARG(ii) = marg_y;
    %[marg_y,res] = get_marginal(Ystd,K,[],eps,tau,l,nugget,EB);
    MARG(ii) = marg_y;
end
plot(MARG)

MARGK20 = MARG;

K_max = 100; % Maximum K
Ystd = Y;
Ystd(:,1) = Ystd(:,1)/max(Y(:,1));
eps =1e-10;
l = .01;
tau = .01;
nugget = 1e-10;
EB = 1;
%(Y,K,s,eps,tau,l,nugget,EB)
MARG = zeros(K_max,1);
for ii = 1:K_max
    K = ii
    for jj=1:length(taus)
        for kk=1:length(ls)
            [marg_y_new,~] = get_marginal(Ystd,K,[],eps,taus(jj),ls(kk),nugget,EB);
            marg_y = max([marg_y, marg_y_new]);
        end
    end
    %MARG(ii) = marg_y;
    %[marg_y,res] = get_marginal(Ystd,K,[],eps,tau,l,nugget,EB);
    MARG(ii) = marg_y;
end
figure()
plot(MARG)




% Now with a nonparametric version
% Create my own survival function
times = [0,1,2,3,4,5,6,7,8,9,10,11];
surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];
%plot(times,surv)
ht = -log(surv);

% Data from first function
theta1 = 1;
surv1 = exp(-ht*theta1);
theta2 = 2;
surv2 = exp(-ht*theta2);
% plot(times,surv1);
% hold on;
% plot(times,surv2);
% hold off;

% Simulate survival times from each
rng(33366)
n1 = 250;
n2 = 250;
X = table([rand(n1,1)/2; rand(n2,1)/2 + .5]);
y = [interp1(surv1,times,rand(n1,1)); interp1(surv2,times,rand(n2,1))];
% Random censoring function.
%cens = rand(n1 + n2,1) + 5; % Original Censoring...
cens = exprnd(50,n1+n2,1);
ind = y < cens;
ds = double(ind);
y(~ind) = cens(~ind);
Y = [y,ds];
[~,yind] = sort(y);
Y = Y(yind,:);
X = X(yind,:);

Ystd = Y;
Ystd(:,1) = Ystd(:,1)/max(Ystd(:,1));

Ystd1 = Ystd(X{:,1} < .5,:);
Ystd2 = Ystd(X{:,1} > .5,:);

% Try several starting points
taus = [.01,.5,1,5];
ls = [.01,.5,.1,5];

MARG = zeros(K_max,1);
for ii = 1:K_max
    K = ii;
    marg_y = [];
    for jj=1:length(taus)
        for kk=1:length(ls)
            [marg_y_new,res] = get_marginal(Ystd,K,[],eps,taus(jj),ls(kk),nugget,EB);
            marg_y = max([marg_y, marg_y_new]);
        end
    end
    MARG(ii) = marg_y;
end
plot(MARG)


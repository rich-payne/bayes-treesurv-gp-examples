function out = get_marginal_vararg(tau,l,varargin)
    varargin = varargin{1};
    Ystd = varargin{1};
    K = varargin{2};
    s = varargin{3};
    eps = varargin{4};
    nugget = varargin{5};
    EB = varargin{6};
    out = get_marginal(Ystd,K,s,eps,tau,l,nugget,EB);
end
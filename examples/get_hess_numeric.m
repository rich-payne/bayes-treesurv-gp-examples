% This function approximates the Hessian of a function
% in two dimensions by supplying the function handle and the
% point at which the Hessian should be approximated.
% Additional arguments (without argnames) to fun can also be specified 
% but MUST BE IN ORDER since fun will use them in order;
% By convention, it is best to put the 'delta_tau' and
% 'delta_l' arguments first if they are to be used.
function [theHess,MAP] = get_hess_numeric(tau,l,fun,varargin) 
    % Custom "parsing" of varargin
    ind_tau = find(strcmp(varargin,'delta_tau'));
    ind_l = find(strcmp(varargin,'delta_l'));

    delta_l = min(.01,l/2);
    delta_tau = min(.01,tau/2);
    if ~isempty(ind_tau)
        delta_tau = varargin{ind_tau + 1};
    end
    if ~isempty(ind_l)
        delta_l = varargin{ind_l + 1};
    end
    ind = [ind_tau , ind_l];
    delind = [ind,ind + 1];
    parameters_passed = varargin;
    parameters_passed(delind) = [];
    
    % Choose points to evaluate the function to estimate Hessian
    thegrid = [tau,l;
        tau + delta_tau, l;
        tau - delta_tau, l;
        tau, l + delta_l;
        tau, l - delta_l;
        tau - delta_tau/2, l - delta_l/2;
        tau + delta_tau/2, l - delta_l/2;
        tau - delta_tau/2, l + delta_l/2;
        tau + delta_tau/2, l + delta_l/2];
    
    % Evaluate the posterior at these points
    ngrid = size(thegrid,1);
    post_vals = zeros(ngrid,1);
    %post_vals(1) = node.Llike;
    for ii=1:ngrid
        post_vals(ii) = fun(thegrid(ii,1),thegrid(ii,2),parameters_passed);
    end
    MAP = post_vals(1);
      
    % Calculate derivatives
    taup1 = (post_vals(1) - post_vals(3)) / delta_tau;
    taup2 = (post_vals(2) - post_vals(1)) / delta_tau;
    taupl = (post_vals(7) - post_vals(6)) / delta_tau;
    taupu = (post_vals(9) - post_vals(8)) / delta_tau;
    % taupcent = (post_vals(1) - post_vals(2))/(2*delta_tau);

    lp1 = (post_vals(1) - post_vals(5)) / delta_l;
    lp2 = (post_vals(4) - post_vals(1)) / delta_l;
    lpl = (post_vals(8) - post_vals(6)) / delta_l;
    lpu = (post_vals(9) - post_vals(7)) / delta_l;
    % lpcent = (post_vals(4) - post_vals(5))/(2*delta_l);
        
    % Calculate second derivatives
    dtaudtau = (taup2 - taup1) / delta_tau;
    dldl = (lp2 - lp1) / delta_l;
    dtaudl = (lpu - lpl) / delta_l;
    dldtau = (taupu - taupl) / delta_tau;
       
    % This sigma corresponds to the one in Rue paper 2009
    theHess = [dtaudtau, dtaudl;
             dldtau, dldl];
    
end




function out = theta_z(z1,z2,tau,l,VsqrtLambda)
    z = [z1; z2];
    out = [tau; l] + VsqrtLambda * z;
end



function get_post(x,y)
    [marg_y,out] = get_marginal(Y,K,s,eps,tau,l,nugget,EB);
end
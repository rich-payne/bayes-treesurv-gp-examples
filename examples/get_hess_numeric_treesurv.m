function [approx_Hess,MAP] = get_hess_numeric_treesurv(Ystd,node,thetree)
    if node.Updatellike
        error('Node must have updated llike before calling this function.');
    end
    % Posterior mode
    tau0 = node.tau;
    l0 = node.l;
    
    Ynode = Ystd(node.Xind,:);
    K = thetree.K;
    s = [];
    eps = thetree.eps;
    nugget = thetree.nugget;
    EB = 0;  

    [approx_Hess, MAP] = get_hess_numeric(tau0,l0,@get_marginal_vararg,Ynode,K,s,eps,nugget,EB);
end
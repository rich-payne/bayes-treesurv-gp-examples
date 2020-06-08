function [fullmarg,ZZ] = get_fullmarginal(Ystd,thenode,thetree,delta_z,delta_pi)
    tau0 = thenode.tau;
    l0 = thenode.l;
    
    % This Sigma is not the same as the covariance function;
    %   This sigma is from the paper by Rue (2009)
    [theHess,MAP] = get_hess_numeric_treesurv(Ystd,thenode,thetree);
    Sigmainv = -theHess;
    % Enforce Symmetry
    dxdyavg = (Sigmainv(1,2) + Sigmainv(2,1))/2;
    Sigmainv(1,2) = dxdyavg;
    Sigmainv(2,1) = dxdyavg;
    %Sigmainv
    %Sigma = inv(Sigmainv)
    %Sigma * Sigmainv
    %Sigmainv * Sigma
    inv(Sigmainv)
    [V,Lambda] = eig(inv(Sigmainv));
    %[V0,Lambdainv] = eig(Sigmainv);
    % Lambda = diag(1 ./ diag(Lambdainv));
    % Lambda = Lambdainv;
    % Lambda = diag(1./diag(Lambdainv));
    %V_sqrtLambda = V0' \ Lambda
    %inv(Sigmainv)
    %V0' \ Lambda * inv(V0)
    
    %V
    %Lambda
    %V
    %Lambda
    V_sqrtLambda = V * sqrt(Lambda)
    
    %theHess
    
    Ynode = Ystd(thenode.Xind,:);
    K = thetree.K;
    s = [];
    eps = thetree.eps;
    nugget = thetree.nugget;
    EB = 0; 
    
    done = 0;
    ii = 1;
    cntr = 2;
    nz = 100;
    ZZ = zeros(nz,3); % first two columns are tau and l, last column is the posterior value
    % Put in the MAP estimate
    ZZ(1,:) = [0,0,thenode.Llike];
    
    
    % Do this in the positive direction of tau;
    while ~done
        z = ii*delta_z;
        theta_z1 = get_thetaz(tau0,l0,z,0,V_sqrtLambda);
        %theta_z1 = [tau0; l0] + V_sqrtLambda * [ii*delta_z; 0];
        if any(theta_z1 <= 0)
            break;
        end
        m1 = get_marginal(Ynode,K,s,eps,theta_z1(1),theta_z1(2),nugget,EB);
        ZZ(cntr,:) = [z, 0, m1];
        ii = ii + 1;
        cntr = cntr + 1;
        if (MAP - m1) > delta_pi
            done = 1;
        end
    end
    
    % Do this in the negative direction of tau;
    done = 0;
    ii = 1;
    while ~done
        z = -ii*delta_z;
        theta_z1 = get_thetaz(tau0,l0,z,0,V_sqrtLambda);
        if any(theta_z1 <= 0)
            break;
        end
        %theta_z1 = [tau0; l0] + V_sqrtLambda * [-ii*delta_z; 0];
        m1 = get_marginal(Ynode,K,s,eps,theta_z1(1),theta_z1(2),nugget,EB);
        ZZ(cntr,:) = [z, 0, m1];
        ii = ii + 1;
        cntr = cntr + 1;
        if (MAP - m1) > delta_pi
            done = 1;
        end
    end
    
    % Do this in the positive direction of l;
    done = 0;
    ii = 1;
    while ~done
        z = ii*delta_z;
        theta_z1 = get_thetaz(tau0,l0,0,z,V_sqrtLambda);
        if any(theta_z1 <= 0)
            break;
        end
        %theta_z1 = [tau0; l0] + V_sqrtLambda * [0; ii*delta_z];
        m1 = get_marginal(Ynode,K,s,eps,theta_z1(1),theta_z1(2),nugget,EB);
        ZZ(cntr,:) = [0, z, m1];
        ii = ii + 1;
        cntr = cntr + 1;
        if (MAP - m1) > delta_pi
            done = 1;
        end
    end
    
    % Do this in the negative direction of l;
    done = 0;
    ii = 1;
    while ~done
        z = -ii*delta_z;
        %theta_z1 = [tau0; l0] + V_sqrtLambda * [0; -ii*delta_z];
        theta_z1 = get_thetaz(tau0,l0,0,z,V_sqrtLambda);
        if any(theta_z1 <= 0)
            break;
        end
        m1 = get_marginal(Ynode,K,s,eps,theta_z1(1),theta_z1(2),nugget,EB);
        ZZ(cntr,:) = [0,z,m1];
        ii = ii + 1;
        cntr = cntr + 1;
        if (MAP - m1) > delta_pi
            done = 1;
        end
    end
      
    minlim = min(ZZ(:,1:2));
    maxlim = max(ZZ(:,1:2));
    
    tauind = minlim(1):maxlim(1);
    tauind(tauind == 0) = [];
    for ii=tauind
        % for a given tau, do positive l values
        for jj=1:maxlim(2)
            z = jj*delta_z;
            theta_z1 = get_thetaz(tau0,l0,ii,z,V_sqrtLambda);
            if any(theta_z1 <= 0)
                break;
            end
            m1 = get_marginal(Ynode,K,s,eps,theta_z1(1),theta_z1(2),nugget,EB);
            ZZ(cntr,:) = [ii,z,m1];
            cntr = cntr + 1;
%             if (MAP - m1) > delta_pi
%                 break;
%             end
        end
        % for a given tau, do negative l values
        for jj=-(1:abs(minlim(2)))
            z = jj*delta_z;
            theta_z1 = get_thetaz(tau0,l0,ii,z,V_sqrtLambda);
            if any(theta_z1 <= 0)
                break;
            end
            m1 = get_marginal(Ynode,K,s,eps,theta_z1(1),theta_z1(2),nugget,EB);
            ZZ(cntr,:) = [ii,z,m1];
            cntr = cntr + 1;
%             if (MAP - m1) > delta_pi
%                 break;
%             end
        end
    end
    
    % shrink ZZ if necessary (if overallocated)
    if size(ZZ,1) > (cntr-1)
        ZZ = ZZ(1:(cntr-1),:);
    end
    
    %ZZ = sortrows(ZZ,[1,2]);
    %ZZ
     
    fullmarg=[];
    
    %interp2(
    %interp_func = scatteredInterpolant(ZZ(:,1),ZZ(:,2),exp(ZZ(:,3)));
    % Now create an interpolation function and integrate.
    %myinterp = @(z1,z2) interp_func(z1,z2); %z1 is z for tau, z2 is tau for l;
    %fullmarg = integral2(myinterp,minlim(1),maxlim(1),minlim(2),maxlim(2));
end


function out = get_thetaz(tau0,l0,z_tau,z_l,V_sqrtLambda)
    out = [tau0; l0] + V_sqrtLambda * [z_tau; z_l];
    % outputs a vector of [tau; l]
end


    
    
    
    
     
    
    
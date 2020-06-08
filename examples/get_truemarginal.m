function [truemarginal,marges,X,Y] = get_truemarginal(Ystd,tau,l,dtau,dl,ntau,nl,K)
    gridtau = linspace(tau - dtau .* ntau/2,tau + dtau .* ntau/2,ntau);
    gridl = linspace(l - dl .* nl/2,l + dl .* nl/2,nl);

    % Restrict to only positive values of tau and l
    gridtau(gridtau <= 0) = [];
    gridl(gridl <= 0) = [];

    %thegrid = combvec(gridtau,gridl);
    [X,Y] = meshgrid(gridtau,gridl);

    marges = zeros(size(X,1),size(X,2));
    for ii=1:size(X,1)
        for jj=1:size(X,2)
            tmpmarg = get_marginal(Ystd,K,[],1e-10,X(ii,jj),Y(ii,jj),0,0);
            marges(ii,jj) = tmpmarg;
        end
    end
    truemarginal = dtau .* dl .* sum(sum(exp(marges)));
end
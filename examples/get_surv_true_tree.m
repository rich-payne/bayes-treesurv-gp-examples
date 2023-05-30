function survs = get_surv_true_tree(x, t)
    times = [0,1,2,3,4,5,6,7,8,9,10,11];
    surv = [1,.9,.85,.5,.45,.44,.43,.1,.09,.05,.01,0];

    I1 = ismember(x{:, 2},{'A','B'}) & x{:, 1} > 5;
    I2 = ismember(x{:, 2},{'A'}) & x{:, 1} <= 5;
    I3 = ismember(x{:, 2},{'B'}) & x{:, 1} <= 5;
    I4 = ismember(x{:, 2},{'C','D'}) & x{:, 1} <= 3;
    I5 = ismember(x{:, 2},{'C','D'}) & x{:, 1} > 3 & x{:, 1} <= 7;
    I6 = ismember(x{:, 2},{'C','D'}) & x{:, 1} >  7;
    I = I1 + 2 * I2 + 3 * I3 + 4 * I4 + 5 * I5 + 6 * I6;
    
    surv1 = 1 - wblcdf(t, 5, 2);
    surv2 = 1 - wblcdf(t, 1, 5);
    surv3 = 1 - wblcdf(t, .5, .9);
    surv4 = 1 - wblcdf(t, 5, 5);
    surv5 = 1 - wblcdf(t, .5, .5);
    surv6 = [interp1(times, surv, t)]; 

    survs = zeros(length(t), size(x, 1));    
    for ii=1:size(survs, 1)
       if I(ii) == 1
          tmp = surv1;
       elseif I(ii) == 2
          tmp = surv2;
       elseif I(ii) == 3
           tmp = surv3;
       elseif I(ii) == 4
           tmp = surv4;
       elseif I(ii) == 5
           tmp = surv5;
       elseif I(ii) == 6
           tmp = surv6;
       end
       survs(:, ii) = tmp;
    end
end
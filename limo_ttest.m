function [m,dfe,ci,sd,n,t,p] = limo_ttest(type,data1,data2,alphav)

% [m,dfe,ci,sd,n,t,p] = limo_ttest(type,data1,data2,alpha)
% implement the one-sample t-test, paired t-test (type = 1) 
% and the two samples t-test (type = 2). The ttest is performed on the last
% non-singleton dimension, which can be subjects or trials. If vectors are
% compared they should have dimensions (1,N).
% 
% INPUTS:
%
% type  = 1 for one sample and paired t-test
%       = 2 for two samples t-test
% data1 = a matrix of data - the test is performed on the last non-singleton dimension
% data2 = a matrix of data or 0 for the one-sample t-test
% alphav = test performed at the (100*alpha)% significance level
%
% OUTPUTS:
%
% [m,dfe,ci,sd,n,t,p] = means, dfe, confidence interavals, std, 
%                   number of observations, t values, p values 
%
% Cyril Pernet 28-08-2009
% GAR 02-09-2009: made dimension flexible and updated documentation
% Cyril 07-09-2009 changed the code to return the dfe
% -----------------------------
%  Copyright (C) LIMO Team 2010

nd = ndims(data1);

switch (type)
 
    case(1)
        if data2 ~= 0
            if ~any(size(data1)==size(data2))
                error('data1 and data2 are of different dimensions')
            else
                data1 = data1 - data2; % paired t-test
            end
        end

        n = size(data1,nd);
        m = mean(data1,nd);
        sd = std(data1,0,nd);
        t = m ./ (sd ./ sqrt(n));
        c = tinv((1 - alphav / 2), n - 1) .* (sd ./ sqrt(n));
        ci = [(m - c) (m + c)];
        dfe = n-1;
        p = 2 * tcdf(-abs(t), n - 1); % two tailed     
        
    case(2)
       
        n = NaN(2,1);
        n(1) = size(data1,nd);
        n(2) = size(data2,nd);

        m = (mean(data1,nd) - mean(data2,nd));
        s1 = var(data1,0,nd) ./ n(1);
        s2 = var(data2,0,nd) ./ n(2);
        se = sqrt(s1 + s2);
        sd = sqrt([s1.*n(1); s2.*n(2)]);
        t = m ./ se;
        
        dfe = ((s1 + s2).^2) ./ ((s1.^2 ./ (n(1)-1) + s2.^2 ./ (n(2)-1)));
        p = 2 * tcdf(-abs(t),dfe); % two tailed
        c = tinv(1-alphav./2, dfe) .* se;
        ci = [(m - c) (m + c)];
end



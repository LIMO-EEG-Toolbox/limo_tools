function [m,dfe,ci,sd,n,t,p] = limo_bootttest1(varargin)

% implement the one-sample t-test using bootstrap.
% The ttest is performed on the last non-singleton dimension,
% which can be subjects or trials. If vectors are compared they
% should have dimensions (1,N).
%
% [m,dfe,ci,sd,n,t,p] = limo_bootttest1(data,alpha,boot)
%
% INPUTS:
%
% data = a matrix of data to be used in the one sample t-test
% alpha = test performed at the (100*alpha) significance level
% boot  = matrix of bootstrap index or N the number of bootstrap to do.
%
% OUTPUTS:
%
% [m,ci,sd,n,t,p] = means, dfe, confidence interavals, std,
%                   number of observations, t values, p values
%
% Based on limo_ttest
% Cyril 01-03-2011
% removed electrode loop in 3D case: GAR 06-12-2011
% -----------------------------
%  Copyright (C) LIMO Team 2010
 
%% check inputs
if nargin == 1
    data = varargin{1};
    nd = numel(size(data));
    alpha = 5/100;
    Nboot = 1000;
    if isvector(data)
        boot = ceil(rand(Nboot,length(data)).*length(data))';
    else
        boot = ceil(rand(Nboot,size(data,nd)).*size(data,nd))';
    end
elseif nargin == 2
    data = varargin{1};
    nd = numel(size(data));
    alpha = varargin{2};
    Nboot = 1000;
    boot = ceil(rand(Nboot,size(data,nd)).*size(data,nd))';
elseif nargin == 3
    data = varargin{1};
    nd = numel(size(data));
    alpha = varargin{2};
    tmp = varargin{3};
    if size(tmp,1) == 1 && size(tmp,2) == 1
        Nboot = tmp;
        boot = ceil(rand(Nboot,size(data,nd)).*size(data,nd))';
    elseif size(tmp,1) == size(data,nd)
        Nboot = size(tmp,2);
        boot = tmp;
    elseif size(tmp,2) == size(data,nd)
        disp('boot matrix transposed')
        boot = tmp';
        Nboot = size(boot,2);
    else
        error('error in boot input')
    end
else
    error('wrong number of arguments')
end
clear tmp

if nd > 3
    error('max data dimension = 3')
end


%% compute

n = size(data,nd);
dfe = n-1;

switch (nd)
    case(1) % 1 dim
        m = NaN(Nboot);
        sd = NaN(Nboot);
        t = NaN(Nboot);
        ci = NaN(Nboot,2);
        p = NaN(Nboot);
        
        for B=1:Nboot
            boot_data = data(boot(:,B));
            m(B) = mean(boot_data);
            sd(B) = std(boot_data,0);
            t(B) = m(B) ./ (sd(B) ./ sqrt(n));
            c = tinv((1 - alpha / 2), n - 1) .* (sd(B) ./ sqrt(n));
            ci(B,:) = [(m(B) - c) (m(B) + c)];
            dfe(B) = n-1;
            p(B) = 2 * tcdf(-abs(t(B)), n - 1); % two tailed
        end
        
    case(2) % 2 dim
        m = NaN(size(data,1),Nboot);
        sd = NaN(size(data,1),Nboot);
        t = NaN(size(data,1),Nboot);
        ci = NaN(size(data,1),Nboot,2);
        p = NaN(size(data,1),Nboot);
        
        for B=1:Nboot
            boot_data = squeeze(data(:,boot(:,B)));
            m(:,B) = mean(boot_data,nd);
            sd(:,B) = std(boot_data,0,nd);
            t(:,B) = m(:,B) ./ (sd(:,B) ./ sqrt(n));
            c = tinv((1 - alpha / 2), n - 1) .* (sd(:,B) ./ sqrt(n));
            ci(:,B,:) = [(m(:,B) - c) (m(:,B) + c)];
            p(:,B) = 2 * tcdf(-abs(t(:,B)), n - 1); % two tailed
        end
        
    case(3) % 3 dim
        m = NaN(size(data,1),size(data,2),Nboot);
        sd = NaN(size(data,1),size(data,2),Nboot);
        t = NaN(size(data,1),size(data,2),Nboot);
        ci = NaN(size(data,1),size(data,2),Nboot,2);
        p = NaN(size(data,1),size(data,2),Nboot);
        
        for B=1:Nboot
            boot_data = squeeze(data(:,:,boot(:,B)));
            m(:,:,B) = mean(boot_data,nd);
            sd(:,:,B) = std(boot_data,0,nd);
            t(:,:,B) = m(:,:,B) ./ (sd(:,:,B) ./ sqrt(n));
            c = tinv((1 - alpha / 2), n - 1) .* (sd(:,:,B) ./ sqrt(n));
            ci(:,:,B,1) = (m(:,:,B) - c);
            ci(:,:,B,2) = (m(:,:,B) + c);
            p(:,:,B) = 2 * tcdf(-abs(t(:,:,B)), n - 1); % two tailed
        end

end



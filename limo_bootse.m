function bse = limo_bootse(varargin)

% Compute a bootstrapped standard error
% We estimate the standard error of the est statistic by the standard
% deviation of the bootstrap replications Efron & Tibshinari 1993, chapter 6
%
% INPUT
% limo_bse = bootse(x,nboot,est,q)
%            x is a vector
%            nboot the number of bootstrap to perform
%            est is the estimator to use, either Harrell-Davis (='limo_hd')
%            or median (='median')
%            q is the decile if est = limo_hd
%
% Wilcox 2005, p.44-45 - See Wilcox p.44-46
% GA Rousselet, University of Glasgow, Dec 2007
% C Pernet make it work with LIMO_EEG
% Version 1 June 2010
% -----------------------------
%  Copyright (C) LIMO Team 2010

x = varargin{1};
nboot = varargin{2};
est = varargin{3};
if length(varargin) == 4
    q = varargin{4};
end

rand('state',sum(100*clock));
n=length(x); % number of subjects/tirals
boot = zeros(1,nboot);
switch lower(est)
    case {'hd'} % Harrell-Davis estimator
        for kk=1:nboot % do bootstrap
            eval(['boot(kk)=',est,'(randsample(x,n,true),q);'])
        end
    otherwise
        for kk=1:nboot % do bootstrap
            eval(['boot(kk)=',est,'(randsample(x,n,true));'])
        end
end
bse = std(boot,0); % normalize by (n-1)


end


%% Harell-Davis subfunction

function thetaq = hd(x,q)

% Compute the Harrell-Davis estimate of the qth quantile
% The vector x contains the data, and the desired quantile is q
% The default value for q is .5.
% % Original R code by Rand Wilcox
% See Wilcox p.71, 139
% GAR, University of Glasgow, Dec 2007

n=length(x);
m1=(n+1).*q;
m2=(n+1).*(1-q);
vec=1:length(x);
w=betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
y=sort(x);
thetaq=sum(w(:).*y(:));

end
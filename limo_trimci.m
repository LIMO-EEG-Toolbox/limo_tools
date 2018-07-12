function [t,tmdata,trimci,se,p,tcrit,df]=limo_trimci(data, percent, alphav, nullvalue)

% [t,tmdata,trimci,se,p,tcrit,df]=limo_trimci(data, percent, alpha, nullvalue)
%
% This function is based on Rand Wilcox's trimci R function.
% Computes a one-sample ttest at all electrodes and time points using 
% trimmed means and windsorised standard error.
%
% This code assumes no missing values.
% 
% INPUT: 
%        data (electrodes*frames*trials)
%        percent (percentage of trimming [0-100], default=20)
%        alphav (default=.05)
%        nullvalue (mean of nullhypothesis, default=0)
%
% OUTPUT:
%        t      (2D matrix of t statistics at every electrode and time frame)
%        tmdata (2D matrix of trimmed means)
%        trimci (confidence intervals around the trimmed means)
%        se      standard error
%        p      (p values)
%        tcrit  (1-alpha/2 quantile of the Student's t distribution with
%               adjusted degrees of freedom)
%        df     (degrees of freedom)
%
%   See also LIMO_TTEST, LIMO_TRIMMED_MEAN, LIMO_WINVAR, LIMO_YUEN_TTEST, LIMO_YUEND_TTEST.
%
% Luisa Frei, 23/11/2011: wrote first version
% GAR, 23/11/2011: edited calculation of confidence intervals & help 
% -----------------------------
%  Copyright (C) LIMO Team 2011

%% data checking
if nargin<4;nullvalue=0;end
if nargin<3;alphav=.05;end
if nargin<2;percent=20;end

if (percent >= 100) || (percent < 0)
    error('limo_trimci:InvalidPercent', 'PERCENT must be between 0 and 100.');
end

if percent == 50
    error('limo_trimci:InvalidPercent', 'PERCENT cannot be 50, use a method for medians instead.');
end

%% ttest
% assumes there are no NaNs in the data

% number of trials
na=size(data,3);

% number of items to winsorize and trim
ga=floor((percent/100)*na);

% windsorize variance of data
asort=sort(data,3);
wa=asort;
wa(:,:,1:ga+1)=repmat(asort(:,:,ga+1),[1 1 ga+1]);
wa(:,:,na-ga:end)=repmat(asort(:,:,na-ga),[1 1 ga+1]);
wva=var(wa,0,3);

% get standard error and df
se=sqrt(wva)/((1-2*(percent/100))*sqrt(na));

df=na-2*floor((percent/100)*na)-1;

% trimmed mean and CI
tmdata=mean(asort(:,:,(ga+1):(na-ga)),3);
trimci(:,:,1)=tmdata+tinv(alphav./2,df).*se; %
trimci(:,:,2)=tmdata-tinv(alphav./2,df).*se; %

t=(tmdata-nullvalue)./se;
p=2.*(1-tcdf(abs(t),df)); % 2-tailed probability

tcrit=tinv(1-alphav./2,df); % 1-alpha/2 quantile of Student's distribution with df degrees of freedom




function TM = limo_trimmed_mean(varargin)

% Compute a trimmed mean and CI on the 3rd dimension of a 3D matrix
%
% INPUT
% TM = limo_trimmed_mean(data,percent,alpha)
%      data is a 3D matrix
%      percent [0-100] is the percentage of trimming to do
%        e.g. 20 will trim 20% on each side, thus a total of 40% of the data.
%      alpha is the error rate of the confidence interval
%      if alpha is not used, only the trimmed mean is computed
%
% OUTPUT
% TM is a 2D or 3D matrix with the lower CI, the trimmed mean and the high CI

% -----------------------------
%  Copyright (C) LIMO Team 2010

% Original R code by Rand Wilcox - See also Wilcox p.71, 139
% GA Rousselet, University of Glasgow, Dec 2007
% C. Pernet rewrote for 3D matrices 
% Version 1 June 2010



%% checkings

data = varargin{1};
percent = 20/100;
if nargin > 1
    percent = varargin{2};
end

if nargin == 3
    option = 1;
    alpha = varargin{3};
else
    option = [];
end

reduced_dim = 0;
if length(size(data)) ==2
    tmp = zeros(2,size(data,1),size(data,2));
    tmp(1,:,:) = data; tmp(2,:,:) = data;
    clear data; data = tmp; clear tmp
    reduced_dim = 1;
end

if length(size(data)) ~=3
    error('data in must be a 2D or 3D matrix')
end


%% do the trimming
n = size(data,3); 
g=floor((percent/100)*n); % g trimmed elements
datasort=sort(data,3);
if option == 1
    TM = NaN(size(data,1),size(data,2),3);
    TM(:,:,2) = nanmean(datasort(:,:,(g+1):(n-g)),3);
else
    TM = NaN(size(data,1),size(data,2));
    TM = nanmean(datasort(:,:,(g+1):(n-g)),3);
end
    
%% compute confidence intervals
if option == 1
    for i=1:size(data,1)
            tmp_data = squeeze(datasort(i,:,:)); tmp = tmp_data(~isnan(tmp_data));
            [tv,g] = tvar(reshape(tmp,size(tmp_data,1),length(tmp)/size(tmp_data,1)),percent); % trimmed squared standard error + g trimmed elements
            se = sqrt(tv); % trimmed standard error
            df = n - 2.*g - 1; % n-2g = number of observations left after trimming
            TM(i,:,1) = TM(i,:,2)+tinv(alpha./2,df).*se; %
            TM(i,:,3) = TM(i,:,2)-tinv(alpha./2,df).*se; %
    end
end

if reduced_dim == 1
   TM = TM(1,:,:);
end

end


function [tv,g]=tvar(x,percent)

% function [tv]=tvar(x,percent)
% The trimmed variance is calculated from the winsorized variance, vw,
% according to the formula vw/(k*length(x)), where k=(1-2*percent/100)^2.
% Original code provided by Prof. Patrick J. Bennett, McMaster University
% See Rand R. Wilcox (2001), Fundamentals of Modern Statisical Methods, page 164.
% See also Rand R. Wilcox (2005), p.61-63
% Edit input checks: GA Rousselet - University of Glasgow - Nov 2008
% Merge function tvar, winvar, winsample - make it work 3D - C Pernet June 2010

g=floor((percent/100)*size(x,2));
xsort=sort(x,2);
loval=xsort(:,g+1);
hival=xsort(:,size(x,2)-g);
for i=1:size(x,1)
    x(i,find(x(i,:)<=loval(i))) = loval(i);
    x(i,find(x(i,:)>=hival(i))) = hival(i);
    wv(i)=var(x(i,:));
end
k=(1-2*percent/100)^2;
tv=wv/(k*length(x));

end




function HD = limo_harrell_davis(varargin)

% Compute the Harrell-Davis estimate and 95% CI of the qth quantile
%
% HD = limo_harrell_davis(data,q,nboot)
%
% INPUT
%      data is a 3D matrix
%      q is a quantile [0-1], e.g. enter 0.1 for 10%
%      nboot if used determines the nb of bootstraps to compute the
%      standard error in limo_bootse (200 might be sufficient)
%
% OUTPUT
% HD is a 2D or 3D matrix with the lower CI, the median and the high CI

% -------------------------------------------
% Copyright (C) LIMO Team 2010

% See Wilcox, 2005, pages 71, 130-134, 139
% GA Rousselet, University of Glasgow, Sept 2009
% C Pernet, UoE, May 2010 - make it work in 3D



%% checking
data = varargin{1};
reduced_dim = 0;
if length(size(data)) ==2
    tmp = zeros(2,size(data,1),size(data,2));
    tmp(1,:,:) = data; tmp(2,:,:) = data;
    clear data; data = tmp; clear tmp
    reduced_dim = 1;
end

if length(size(data)) ~=3
    error('data must be a 3D matrix')
end

q = varargin{2};

if length(varargin) == 3
    option = 1;
    nboot = varargin{3};
else
    option = [];
end

%% compute Harell-Davis of q
n=size(data,3); vec=1:n;
m1=(n+1).*q; m2=(n+1).*(1-q);
w=betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
datasort=sort(data,3);
if isempty(option)
    HD = NaN(size(data,1),size(data,2));
    for i=1:size(data,1)
        for j=1:size(data,2)
            tmp_data = squeeze(datasort(i,j,:));
            HD(i,j,:)=sum(w(~isnan(tmp_data))'.*tmp_data(~isnan(tmp_data)));
        end
    end
end

%% compute HD of q and 95% confidence intervals
% The constant c was determined so that the
% probability coverage of the confidence interval is
% approximately 95% when sampling from normal and
% non-normal distributions - from Wilcox

if option == 1

    if size(data,3)<=10
        
        fprintf('***************************************************\n');
        fprintf('Error in limo_harrell_davis:\n');
        fprintf('confidence intervals of the hd estimates of the deciles cannot be computed for less than 11 observations\n')
        fprintf('***************************************************\n');
    
    else

        c = 1.96 + .5064.* (size(data,3).^-.25);
        if q<=2 || q>=8
            if size(data,3) <=21
                c = -6.23./size(data,3)+5.01;
            end
        end

        if q<=1 || q>=9
            if  size(data,3)<=40
                c = 36.2./size(data,3)+1.31;
            end
        end

        for i=1:size(data,1)
            fprintf('processing electrode %g',i); disp(' ');
            for j=1:size(data,2)
                tmp_data = squeeze(datasort(i,j,:));
                HD(i,j,2)=sum(w(~isnan(tmp_data))'.*tmp_data(~isnan(tmp_data)));
                xd_bse = limo_bootse(tmp_data(~isnan(tmp_data)),nboot,'hd',q);
                HD(i,j,1) = HD(i,j,2)+c.*xd_bse;
                HD(i,j,3) = HD(i,j,2)-c.*xd_bse;
            end
        end
    end
end

if reduced_dim == 1
   HD = HD(1,:,:);
end



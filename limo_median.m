function M = limo_median(varargin)

% Compute a median and 95% CI on the dimension 3 of a 3D matrix
%
% INPUT
% M = limo_median(varargin)
%      data is a 3D matrix
%      nboot if used determines the nb of bootsrtraps to compute the
%      standard error in limo_bootse
%
% OUTPUT
% M is a 2D or 3D matrix with the lower CI, the median and the high CI
%
% GA Rousselet, University of Glasgow, Dec 2007
% C. Pernet Rewritten for 3D matrice in integration of variance etc .. 
% see also limo_bootse - Version 1 June 2010
% -----------------------------
%  Copyright (C) LIMO Team 2010

%% checking
data = varargin{1};
if length(size(data)) ~=3
    error('data in must be 3D matrices')
end

if length(varargin) == 2
    option = 1;
    nboot = varargin{2};
else
    option = [];
end

%% compute median

if option == 1
    M = NaN(size(data,1),size(data,2),3);
    M(:,:,2) = nanmedian(data,3);
else
    M = NaN(size(data,1),size(data,2));
    M = nanmedian(data,3);
end

%% compute the 95% confidence intervals
% The constant c was determined so that the
% probability coverage of the confidence interval is
% approximately 95% when sampling from normal and
% non-normal distributions - from Wilcox
%
% An alternative method is described in Wilcox 2005, p. 132-134, but not implemented here.

if option == 1
    c = 1.96 + .5064.* (size(data,3).^-.25);
    for i=1:size(data,1)
        fprintf('processing electrode %g',i); disp(' ');
        for j=1:size(data,2)
            tmp_data = squeeze(data(i,j,:));
            xd_bse = limo_bootse(tmp_data(~isnan(tmp_data)),nboot,'median');
            M(i,j,1) = M(i,j,2)+c.*xd_bse;
            M(i,j,3) = M(i,j,2)-c.*xd_bse;
        end
    end
end

    
    

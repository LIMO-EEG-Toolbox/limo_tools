function [b w] = limo_IRLS(X,Y,tune)

% LIMO_IRLS Limo Iteratively Reweighted Least Squares (IRLS)
% IRLS is used to find the maximum likelihood estimates of a generalized
% linear model, and in robust regression to find an M-estimator, as a way 
% of mitigating the influence of outliers in an otherwise 
% normally-distributed data set. 
%
% FORMAT: [b w] = limo_IRLS(X,Y,tune)

% INPUTS:
%   X             = the design matrix 
%   Y             = 2D matrix of EEG data (dim trials x frames)
%  tune           = Tuning constant; Decreasing the tuning constant 
%                   increases the downweight assigned to large residuals; 
%                   increasing the tuning constant decreases the downweight
%                   assigned to large residuals.
%
% OUTPUTS:
%   b             = betas (dim parameters * time frames)
%   w             = weights (dim trials * time frames)
%
% References:
%   J�rgen Gro� (2003), "Linear Regression" pp 191-215
%   Street, J.O., R.J. Carroll, and D. Ruppert (1988), "A note on
%     computing robust regression estimates via iteratively
%     reweighted least squares," The American Statistician,
%   Wager TD, Keller MC, Lacey SC, Jonides J. (2005 May 15)
%    "Increased sensitivity in neuroimaging analyses using robust
%    regression".Neuroimage. 26(1):99-113
%    
% Ignacio Suay Mas
% -----------------------------
% Copyright (C) LIMO Team 2015

if  nargin < 2      
    error(message('Too Few Inputs'));      
end 

if nargin == 2    
 tune = 4.685; % tuning function for the bisquare
end 

[rows,cols] = size(X);
if (rows <= cols)
   error('IRLS cannot be computed, there is not enough trials for this design');     
end

% Find Ordinary Least Squares
b = pinv(X)*Y;

% H - Hat matrix; Leverages for each observation
H = diag(X*pinv(X'*X)*X');
% Adjustment factor
adjfactor = 1 ./ sqrt(1-H);
adjfactor(adjfactor==Inf) = 1; % when H=1 do nothing

numiter = 0; iterlim = 100; % set a 100 iteration max
oldRes=1; newRes=10;

while(max(abs(oldRes-newRes)) > (1E-4))
   
   numiter = numiter+1;
   oldRes = newRes;
   
   if (numiter>iterlim)
      warning('Limo_IRLS could not converge');
      break;
   end
   
   % Get residuals from previous fit
   res = Y - X*b;
   resadj = res .* repmat(adjfactor, 1, size(Y,2));

   %re - Robust Estimator
   % 0.6745 is the 0.75- quantile of the standard normal distribution
   % (makes the estimate unbiased)
   re = median(abs(resadj)) ./ 0.6745;
   re(find(re < 1e-5)) = 1e-5;
   r= resadj ./ repmat(tune.*re, size(Y,1),1);
   
   % Compute new weights 
   w = (abs(r)<1) .* (1 - r.^2).^2;
   yw = Y .* w;
   w = sqrt(w);
   
   for i= 1: size(Y,2)
       xw = X .* repmat(w(:,i), 1, size(X,2));
       b(:,i) = pinv(xw)*yw(:,i);
   end
   
   % newRes= sum(sum(res.^2));
   newRes= sum(res(:).^2);
   
%    % plot
%    if numiter == 1
%        figure; xx = zeros(100,1);
%    end
%    xx(numiter) = newRes;
%    plot([1:100],xx); drawnow
   
end

end

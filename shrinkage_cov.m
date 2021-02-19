function [sigma,rho] = shrinkage_cov(X,est)
% this program is distributed under BSD 2-Clause license
%
% Copyright (c) <2016>, <Okba BEKHELIFI, okba.bekhelifi@univ-usto.dz>
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%
%
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
% BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
% OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
% OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%
%compute a covariance matrix estimation using shrinkage estimators: RBLW or
% OAS describer in:
%     [1] "Shrinkage Algorithms for MMSE Covariance Estimation"
%     Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.
%Input:
%X: data NxP : N : samples P: features
%est: shrinkage estimator
%     'rblw' : Rao-Blackwell estimator
%     'oas'  : oracle approximating shrinkage estimator
%     default estimator is OAS.
%
%Output:
% sigma: estimated covariance matrix
% rho : shrinkage coefficient
%
%created: 02/05/2016
%last revised: 14/06/2016


[n,p] = size(X);
% sample covariance, formula (2) in the paper [1]
X = bsxfun(@minus,X,mean(X));
S = X'*X/n;
% structured estimator, formula (3) in the article [1]
mu = trace(S)/p;
F = mu*eye(p);


if (nargin < 2) est = 'oas';
end
switch lower(est)
    
    case 'oas'
        
%         rho = (1-(2/p)*trace(S^2)+trace(S)^2)/((n+1-2/p)*(trace(S^2)-1/p*trace(S)^2));
        c1 = 1-2/p;
        c2 = n+1-2/p;
        c3 = 1-n/p;
        rho = (c1*trace(S^2) + trace(S)^2) / (c2*trace(S^2) + c3*trace(S)^2);
        
    case 'rblw'
        
        c1 = (n-2)/n;
        c2 = n+2;
        rho = ( c1*trace(S^2)+trace(S)^2 )/( c2*( trace(S^2)-(trace(S)^2/p) ) );
        
    otherwise 'Shrinkage estimator not provided correctly';
        
end

% regularization, formula (4) in the paper [1]
sigma = (1-rho)*S + rho*F;

end


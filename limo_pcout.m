function [dist,out,w1,w2] = limo_pcout(x,varargin)

% LIMO_pcout Limo Principal Components Projection
% This method uses simple properties of Principal Components to identify
% outliers in the Multivariate space. This function use Filzmoser, Moronna
% and Werner implementation but have ommited dimensions where MAD = 0.
% 
% FORMAT [dist out] = limo_pcout(x,'figure',option)
%
% INPUTS:
%   x             = 2D matrix of EEG data (dim trials x frames)
%   'figure'      = indicate you want a figure out
%   option        = 'new' or 'on' to plot on a new figure or the current one
%
% OUTPUTS:
%   dist          = weights (distance) for each trial. Outliers have a near to zero weight.
%   out           = The out is a binary vector where False values are outliers
%   w1            = the kurtosis weighted robust Euclidian distance (location)
%   w2            = the robust Euclidian distance (scatter)
%
% References:
% Filzmoser, R. Maronna, and M. Werner. 
% Outlier identification in high dimensions. 
% Computational Statistics and Data Analysis, 
% Vol. 52, pp. 1694-1711, 2008.
%
% Ignacio Suay Mas - Cyril Pernet 
% ------------------------------
%  Copyright (C) LIMO Team 2019

if nargin == 3 
   if strcmpi(varargin{1},'figure')
       fig = varargin{2};
   else
      disp(''); fig = 'off'; 
   end
else
    fig = 'off';
end

% we have ommited dimensions with mad = 0, as this 
% corresponds to a case with all trials having the same value ! 

x = x(:,mad(x,1) > 1e-6); 
if isempty(x)
    error('WLS cannot be computed, for at least 1 frame, all trials have the same values')
end

[n,p] = size(x);
if n<p
    error('Principal Component Projection cannot be computed, more observations than variables are needed')
end

%% 1st Phase: Detect location outliers

% robustly rescale each component
madx    = mad(x,1) .* 1.4826; % median erp normally scaled
x2      = (x - repmat(median(x), n, 1)) ./ repmat(madx,n,1); % robust standardization
x3      = x2 - repmat(mean(x2),n,1); % standardization
 
% calculate the principal components (retain 99% of variance)
a       = svd(x3).^2./(n-1);
paux    = (1:p)' .* ((cumsum(a)./sum(a) > 0.99));
paux2   = paux(paux > 0);
p1      = paux2(1);
[~,~,V] = svd(x3);

% Project x in the principal components
xpc     = x2 * V(:,1:p1); 
madxpc  = mad(xpc,1) * 1.4826;

% rescale pc 
xpcsc   = (xpc - repmat(median(xpc),n,1)) ./ repmat(madxpc,n,1);

%% Location Outliers 
% calculate the absolute value of robust kurtosis measure
% kurtosis = 0 -> no outliers; kurtosis != 0 -> outliers
wp = abs(mean(xpcsc.^4) - 3);
% to scale wp between 0,1 -> wp/sum(wp)
% calculate the pc scaled by weights
xpcwsc = xpcsc * diag(wp/sum(wp));
xpcnorm = sqrt(sum(xpcwsc.^2,2));

%Transform the Robust distance (xpcnorm)
xdist1 = (xpcnorm .*  sqrt(chi2inv(0.5,p1))) ./ median(xpcnorm);

%Obtain an accurate classification between outliers and non-outliers
%translated biweight function 
M1 = quantile(xdist1, 1/3);
const1 = median(xdist1) + 2.5 * mad(xdist1,1)* 1.4826;
w1 = ( 1 - ((xdist1 - M1)./(const1 - M1)).^2).^2;
w1(xdist1 < M1) = 1;
w1(xdist1 > const1) = 0;

%% Scatter outliers. 
% Similar to the first phase but without the kurtosis weigthed function. 

xpcnorm = sqrt(sum(xpcsc.^2,2));
xdist2 = xpcnorm *  sqrt(chi2inv(0.5,p1)) ./ median(xpcnorm);

M2 = sqrt(chi2inv(0.25,p1));
const2 = sqrt(chi2inv(0.99,p1));
w2 = ( 1 - ((xdist2 - M2)./(const2 - M2)).^2).^2;
w2(xdist2 < M2) = 1;
w2(xdist2 > const2) = 0;

% Sometimes too many non-outliers receive a weight of 0, with s != 0
% helps to ensure that outliers are detected in both phases (w1,w2)
s = 0.25;
dist = (w1 + s) .* (w2 + s) ./ ((1 + s)^2);
out = round(dist + s);

if strcmpi(fig,'new') || strcmpi(fig,'on')
    if strcmpi(fig,'new')
        figure('Name','limo_pcout projection')
    end
    plot(xpc'); 
    scale = range(xpc(:))*0.01; grid on; box on
    axis([-0.5 size(xpc,2) min(xpc(:))-scale max(xpc(:))+scale]); ylabel('A.U.');
    xlabel(sprintf('%g time comp / %g time frames)',size(xpc,2),size(x,2)));
    title('Trials projected onto the PC space')
end


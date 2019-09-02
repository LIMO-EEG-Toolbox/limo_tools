function [mask,p_val,max_th] = limo_max_correction(varargin)

% correction for multiple testing using the max stat value
% since the type 1 error is the prob to make at least one error
% we can control the prob the make an error for the biggest value
% note this works for both bootstrapped data/tfce trandformed data under H0  
%
% FORMAT [mask,p_val,max_th] = limo_max_correction(M,bootM,p,fig)
%
% INPUT
% M     = 2D matrix of observed values 
% bootM = 3D matrix of T^2 or F values for data bootstrapped under H0
% p     = threshold to apply e.g. 5/100
% fig   = 1/0 to plot the maximum stat under H0
%
% OUTPUT
% mask is a binary matrix of the same size as M corresponding to a threshold
%      p corrected for multiple comparisons
% p_val are the p-values obtained via the matrix bootM (non-corrected)
% max_th is the threshold controlling the type 1 FWER
%
% outsourced from limo_stat_values
% Cyril Pernet 
% ------------------------------
%  Copyright (C) LIMO Team 2019

% check inputs 
M       = varargin{1};
bootM   = varargin{2};
if nargin == 2
    p   = 0.05;
    fig = [];
elseif nargin == 3
    p   = varargin{3};
    fig = [];
elseif nargin == 4
    p   = varargin{3};
    fig = varargin{4};
end
clear varargin

[a,b,nboot]=size(bootM);
if any(size(M)~=[a b])
    error('dimension error: matrices of observed and bootstrap values are different')
end

% collect highest value for each boot
parfor boot=1:nboot
    data = squeeze(bootM(:,:,boot));
    maxM(boot) = max(data(:)); 
end

% get threshold
maxM(maxM==Inf) = [];
sortmaxM        = sort(maxM); 
nboot           = length(sortmaxM);
U               = round((1-p).*nboot);
max_th          = sortmaxM(U);
mask            = squeeze(M) >= max_th;

% get the equivalent bootstrapped p-value
smalest_pval = 1/nboot;
for row =1:a
    for column=1:b
        tmp = sum(M(row,column) >= sortmaxM) / nboot;
        p_val(row,column) = 1-tmp;
        if p_val(row,column) == 0; p_val(row,column) = smalest_pval; end
    end
end

%% figure
if sum(mask(:)) == 0 && isempty(fig)
    fig = 1 ; 
end

if fig == 1
    figure('Name','Correction by max: results under H0')
    plot(sortmaxM,'LineWidth',3); grid on; hold on; 
    
    plot(min(find(sortmaxM==max_th)),max_th,'r*','LineWidth',5)
    txt = ['bootstrap threashold ' num2str(max_th) '\rightarrow'];
    text(min(find(sortmaxM==max_th)),max_th,txt,'FontSize',12,'HorizontalAlignment','right');
    
    [val,loc]=min(abs(sortmaxM-max(M(:)))); 
    plot(loc,max(M(:)),'r*','LineWidth',5)
    txt = ['maximum observed: ' num2str(max(M(:))) '\rightarrow'];
    text(loc,max(M(:)),txt,'FontSize',12,'HorizontalAlignment','right');
    
    title('Maxima under H0','FontSize',12)
    xlabel('sorted bootstrap iterations','FontSize',12); 
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end

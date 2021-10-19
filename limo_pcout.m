function [dist,out,rf,w1,w2] = limo_pcout(x,varargin)

% LIMO_pcout Limo Principal Components Projection
% This method uses simple properties of Principal Components to identify
% outliers in the Multivariate space. This function use Filzmoser, Moronna
% and Werner implementation but have ommited dimensions where MAD = 0.
% 
% FORMAT dist = limo_pcout(x)
%        [dist,out,rf,w1,w2] = limo_pcout(x,'downsample','on','figure',option)
%
% INPUTS:
%   x             = 2D matrix of EEG data (dim trials x frames)
%   'downsample'  = indicates option of downsampling when there are not
%                   enough trials set to 'on' or 'off' (default)
%   'figure'      = indicates you want a figure out
%   option        = 'new' or 'on' to plot on a new figure or the current one
%   'xaxis'       = label followed by vector value (e.g. 'Time' [-50:4:200])
%   'weightsas'   = 'histograms','kernel' (default)
%
% OUTPUTS:
%   dist          = weights (distance) for each trial. Outliers have a near to zero weight.
%   out           = the out is a binary vector where False values are outliers
%   rf            = the reduction factor (i.e. how many components were removed
%                   to compute w1/w2 and get dist/out)
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

downsamplex = 'on';
fig         = 'off';
weightsas   = 'Kernel';

if nargin >1
    for i=1:nargin-2
        if strcmpi(varargin{i},'downsample')
            downsamplex = varargin{i+1};
        elseif strcmpi(varargin{i},'figure')
            fig = varargin{i+1};
        elseif strcmpi(varargin{i},'weightsas')
            weightsas = varargin{i+1};
        elseif strcmpi(varargin{i},'xaxis')
            if ischar(varargin{i+1})
                label = varargin{i+1};
                vect  = varargin{i+2};
            else
                label = varargin{i+2};
                vect  = varargin{i+1};
            end
        end
    end
end

% we have ommited dimensions with mad = 0, as this 
% corresponds to a case with all trials having the same value ! 

x = x(:,mad(x,1) > 1e-6); 
if isempty(x)
    error('WLS cannot be computed, for at least 1 frame, all trials have the same values')
end

[n,p] = size(x);
if n<p && strcmpi(downsamplex,'on')
    f = p/n; 
    if ceil(f) == 2
        warning('automatically downsampling by 2 for WLS computation as there are not enough trials')
        x = downsample(x',2)'; [n,p] = size(x);
    else
        warning('%g frame/trial ratio is very low, running SVD as usual but consider using OLS rather than WLS especially if you have a low sampling rates',f)
    end
end

%% 1st Phase: Detect location outliers

% robustly rescale each component
madx    = mad(x,1) .* 1.4826;                                % median absolute deviations normally scaled (~std)
x2      = (x - repmat(median(x), n, 1)) ./ repmat(madx,n,1); % robust standardization (like x-mean / std)
x3      = x2 - repmat(mean(x2),n,1);                         % make it 0 centred

% calculate the principal components (retain 99% of variance)
a       = svd(x3).^2./(n-1);                      % take the square of eigen values
paux    = (1:length(a))' .* ((cumsum(a)./sum(a) > 0.99)); % indices to keep 99% of the variance
paux2   = paux(paux > 0);
p1      = paux2(1);                               % number of components to keep 
rf      = p-p1;                                   % number of components to remove

% Project x in the principal components
[~,~,V] = svd(x3);
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
xpcwsc  = xpcsc * diag(wp/sum(wp));
xpcnorm = sqrt(sum(xpcwsc.^2,2));

%Transform the Robust distance (xpcnorm)
xdist1 = (xpcnorm .*  sqrt(chi2inv(0.5,p1))) ./ median(xpcnorm);

%Obtain an accurate classification between outliers and non-outliers
% translated biweight function 
M1                  = quantile(xdist1, 1/3);
const1              = median(xdist1) + 2.5 * mad(xdist1,1)* 1.4826;
w1                  = ( 1 - ((xdist1 - M1)./(const1 - M1)).^2).^2;
w1(xdist1 < M1)     = 1;
w1(xdist1 > const1) = 0;

%% Scatter outliers. 
% Similar to the first phase but without the kurtosis weigthed function. 

xpcnorm             = sqrt(sum(xpcsc.^2,2));
xdist2              = xpcnorm *  sqrt(chi2inv(0.5,p1)) ./ median(xpcnorm);
M2                  = sqrt(chi2inv(0.25,p1));
const2              = sqrt(chi2inv(0.99,p1));
w2                  = ( 1 - ((xdist2 - M2)./(const2 - M2)).^2).^2;
w2(xdist2 < M2)     = 1;
w2(xdist2 > const2) = 0;

% Sometimes too many non-outliers receive a weight of 0, with s != 0
% helps to ensure that outliers are detected in both phases (w1,w2)
s    = 0.25;
dist = (w1 + s) .* (w2 + s) ./ ((1 + s)^2);
out  = round(dist + s);

if strcmpi(fig,'new') || strcmpi(fig,'on')
    if strcmpi(fig,'new')
        figure('Name','limo_pcout projection')
    end
    
    if ~exist('vect','var')
        vect = 1:size(x,2);
    end
    
    % show trials
    subplot(3,5,[1 2 3]); plot(vect,x); title('Single trials');
    grid on; box on; if exist('label','var'); xlabel(label); end
    text(5,  min(x(:))+5, sprintf('%g frames',size(x,2)));
    axis([-0.5 size(x,2) min(x(:)) max(x(:))]);
    if exist('label','var'); xlabel(label); end
    
    subplot(3,5,[4 5]); plot(1:size(xpc,2),xpc); title('Projected trials onto PC space');
    text(5,  min(xpc(:))+5, sprintf('%g components',size(xpc,2)));
    axis([-0.5 size(xpc,2) min(x(:)) max(x(:))]);
    xlabel('Eigen vectors'); grid on; box on; 
    
    
    % show weights
    if strcmpi(weightsas,'Kernel')
        [~,~,~,KDE]=data_plot([w1 w2 dist],'estimator','mean','figure','off');
        subplot(3,5,[6 7]); bar(KDE{1}(:,1),KDE{1}(:,2),'FaceColor',[0.1 0.4 0.7])
        ylabel('frequency'); axis tight; grid on; box on; title('location weights')
        subplot(3,5,[11 12]); bar(KDE{2}(:,1),KDE{2}(:,2),'FaceColor',[0.7 0.8 0.5]); xlabel('weight values')
        ylabel('frequency'); axis tight; grid on; box on; title('scatter weights')
        subplot(3,5,[8 13]); index = KDE{3}(:,1) > min(dist(out==1));
        bar(KDE{3}(~index,1),KDE{3}(~index,2),'FaceColor',[0.25 0.25 0.25],'EdgeColor',[0.25 0.25 0.25],'FaceAlpha',0.8,'EdgeAlpha',0.8);
        set(gca,'view',[90 -90]); hold on; bar(KDE{3}(index,1),KDE{3}(index,2),'FaceColor',[0.9 0.6 0.6]);
        ylabel('frequency'); axis tight; grid on; box on; title('final weights')
    else
        k = round(1 + log2(size(x,1))); % sturges binning
        subplot(3,5,[6 7]); histogram(w1,k,'FaceColor',[0.1 0.4 0.7])
        ylabel('frequency'); axis tight; grid on; box on; title('location weights')
        subplot(3,5,[11 12]); histogram(w2,k,'FaceColor',[0.7 0.8 0.5]); xlabel('weight values')
        ylabel('frequency'); axis tight; grid on; box on; title('scatter weights')
        subplot(3,5,[8 13]); h = histogram(dist,k); BinEdges = h.BinEdges;
        histogram(dist(out==1),BinEdges,'FaceColor',[0.9 0.6 0.6]); set(gca,'view',[90 -90])
        hold on; histogram(dist(out==0),BinEdges,'FaceColor',[0.25 0.25 0.25]);
        ylabel('frequency'); axis tight; grid on; box on; title('final weights')
    end

    % show averages
    [good_mean,good_ci] = limo_central_estimator(x(out==1,:)','Mean',95/100);
    [bad_mean, bad_ci]  = limo_central_estimator(x(out==0,:)','Mean',95/100);
    [weighted_mean,weighted_ci] = limo_central_estimator((x.*repmat(dist,[1 size(x,2)]))','Mean',95/100);
    subplot(3,5,[9 10]); plot(vect,bad_mean,'LineWidth',2,'color',[0.25 0.25 0.25])
    fillhandle = patch([vect fliplr(vect)], [bad_ci(1,:),fliplr(bad_ci(2,:))], [0.25 0.25 0.25]);
    set(fillhandle,'EdgeColor',[0.25 0.25 0.25],'FaceAlpha',0.2,'EdgeAlpha',0.8);
    hold on; plot(vect,good_mean,'LineWidth',2,'color',[0.9 0.6 0.6])
    fillhandle = patch([vect fliplr(vect)], [good_ci(1,:),fliplr(good_ci(2,:))], [0.9 0.6 0.6]);
    set(fillhandle,'EdgeColor',[0.9 0.6 0.6],'FaceAlpha',0.2,'EdgeAlpha',0.8);
    axis tight; grid on; box on; title('Mean of good trials vs. Mean of outliers')
    subplot(3,5,[14 15]); plot(vect,good_mean,'LineWidth',2,'color',[0.9 0.6 0.6])
    fillhandle = patch([vect fliplr(vect)], [good_ci(1,:),fliplr(good_ci(2,:))], [0.9 0.6 0.6]);
    set(fillhandle,'EdgeColor',[0.9 0.6 0.6],'FaceAlpha',0.2,'EdgeAlpha',0.8);
    hold on; plot(vect,weighted_mean,'LineWidth',2,'color',[0.2 0.6 0.2]);
    fillhandle = patch([vect fliplr(vect)], [weighted_ci(1,:),fliplr(weighted_ci(2,:))], [0.2 0.6 0.2]);
    set(fillhandle,'EdgeColor',[0.2 0.6 0.2],'FaceAlpha',0.2,'EdgeAlpha',0.8);
    axis tight; grid on; box on; title('Mean of good trials vs. Weighted Mean')
    drawnow
end


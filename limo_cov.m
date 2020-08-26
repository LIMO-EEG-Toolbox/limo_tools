function W = limo_cov(varargin)

% returns weights for the GLM - however, instead of computing weights
% based on the covariance between (residual) trials we compute the
% mimimum generalized variance based this on the smoothed residuals 
%
% 1st, compute glm residuals
% 2nd, residuals of each electrode are weighted by their direct 
%      neighbourgs as defined by the neighbouring matrix from the
%      montage (allows bad channel to recover some of the lost signal)
%      - the weights are the xcrorr values (with the lag being of 
%      +/- 10 ms max) - computed per condition
% 3rd, the residuals are smoothed in time with a span = FWHM of 
%      the autocorelation function obtained on the grand average (per
%      condition) - the smoothing is perform using the LOESS method
%      ie Locally Weighted Scatterplot Smoothing
% 4th, the covariance distances are computed calling limo_mgv that 
%      computes the minimum generalized variance. 
%
% % FORMAT:
% S = limo_cov(Y,LIMO)
% S = limo_cov(Y,X,nb_conditions)
%
% INPUTS:
%   Y             = 3D matrix of EEG data with format electrodes x frames x trials
%   X             = 2 dimensional design matrix
%   nb_conditions = a vector indicating the number of conditions per factor
%   LIMO          = structure that contains X and nb_conditions 
%
% OUTPUTS:
%   W     Weights = [electrode x trials]
%
% Refs:
%
% Cleveland, W.S. (1979). Robust Locally Weighted Regression and Smoothing
% Scatterplots. J Am Stat Ass, 74, 829-836.
%
% Pernet, C., Rousselet, G. & Bochkina, N. (2012)
% Robust weighted least squares for MEEG analyses.
% NeuroImage xxx
%
% See also
% LIMO_EEG, LIMO_MGV
%
% Cyril Pernet 20-04-2011
% -----------------------------
%  Copyright (C) LIMO Team 2010

%% data check
if nargin == 2 
    Yr = varargin{1};
    LIMO = varargin{2};
elseif nargin == 3
    Yr = varargin{1};
    LIMO.design.X = varargin{2};
    LIMO.design.nb_conditions = varargin{3};
else
    error('wrong number of arguments')
end

clear varargin
disp('Starting covariance estimation')

%% do a quick OLS to get Residuals
disp('compute OLS residuals ..')
X = LIMO.design.X;
R = eye(size(Yr,3)) - (X*pinv(X'*X)*X');
Res = NaN(size(Yr));
array = find(~isnan(Yr(:,1,1)));
for e=1:length(array)
    electrode = array(e);
    Res(electrode,:,:) = (R*squeeze(Yr(electrode,:,:))')';
end

clear Yr P R ; Yr = Res; clear Res;

%% get an average across trials / split data per condition
% ----------------------------------------------------
if LIMO.design.nb_conditions == 0
    avg{1} = mean(Yr,3); % average trials
else
    for condition = 1:LIMO.design.nb_conditions
        avg{condition} = mean(Yr(:,:,find(LIMO.design.X(:,condition))),3);
    end
end


%% Smooth data in space
% -----------------------------
CWYr = Yr;

% 1. Compute the x correlation between relevant electrodes
% -------------------------------------------------------------
disp('computing correlations between channels ...')
corr_matrix = NaN(size(LIMO.data.channeighbstructmat,1),size(LIMO.data.channeighbstructmat,1),size(avg,2));
[row,col] = ind2sub(size(LIMO.data.channeighbstructmat), find(triu(LIMO.data.channeighbstructmat)));
% force the xcorr to be within credible range ie +/- 10ms
frame_duration = 1000/LIMO.data.sampling_rate;
frame_lag = round(10/frame_duration);

for condition = 1:size(avg,2)
    for index = 1:length(col);
        a = row(index); b = col(index);
        if sum(isnan(avg{condition}([a b],1))) == 0
            r = max(xcorr(avg{condition}(a,:),avg{condition}(b,:),frame_lag,'coeff'));
            corr_matrix(a,b,condition) = r;
            corr_matrix(b,a,condition) = r;
        end
    end
end

% 2. get a new data matrix weighted across elecrodes
% --------------------------------------------------
disp('weighting channels ...')
for condition = 1:size(avg,2)
    for electrode = 1:size(Yr,1);
        data = squeeze(Yr(electrode,:,find(LIMO.design.X(:,condition))));
        other_channels = find([1:size(Yr,1)] - electrode);
        weights = corr_matrix(electrode,:,condition);
        channels = find(~isnan(weights));
        neighbours = [];
        if ~isempty(channels)
            for i = 1:length(channels)
                neighbours(:,:,i) = weights(channels(i)) * squeeze(Yr(channels(i),:,find(LIMO.design.X(:,condition))));
            end
            CWYr(electrode,:,find(LIMO.design.X(:,condition))) = ...
                (data + sum(neighbours,3)) ./ (size(neighbours,3)+1);
        end
    end
end

clear Yr a b channels col condition corr_matrix data electrode frame_lag frame_duration
clear avg i index neighbours other_channels r row weights varargin
Yr = CWYr; clear CWYr;

%% smooth data in time
% ---------------------------
disp('smoothing data ...')

if LIMO.design.nb_conditions ~= 0
    
    for electrode = 1:size(Yr,1)
        trial_index = 1;
        for condition = 1:LIMO.design.nb_conditions
            
            % 1 Estimate autocorrelation to smooth data in time
            % ---------------------------------------------------
            Data = squeeze(Yr(electrode,:,find(LIMO.design.X(:,condition))));
            ARf = xcorr(mean(Data,2),'coeff'); % average trials
            ARf = ARf(ceil(length(ARf)/2):end);
            FWHM = min(find(ARf<1/2));
            
            % 2 smooth the data
            % -----------------
            if ~isempty(FWHM)
                for trial = 1:size(Data,2)
                    Yr(electrode,:,trial_index) = (lowess(Data(:,trial),FWHM))';
                    trial_index = trial_index +1;
                end
            end
        end
    end
else % pure regression design
    for electrode = 1:size(Yr,1)
        trial_index = 1;
        
        % 1 Estimate autocorrelation to smooth data in time
        % ---------------------------------------------------
        Data = squeeze(TSYr(electrode,:,:));
        ARf = xcorr(mean(Data,2),'coeff'); % average trials
        ARf = ARf(ceil(length(ARf)/2):end);
        FWHM = min(find(ARf<1/2));
        
        % 2 smooth the data
        % -----------------
        if ~isempty(FWHM)
            for trial = 1:size(Data,2)
                Yr(electrode,:,trial_index) = (lowess(Data(:,trial),FWHM))';
                trial_index = trial_index +1;
            end
        end
    end
end

clear ARf Data FWHM conditon electrode trial trial_index

%% change sampling rate to get more trials than time frames
 
[e,f,t]=size(Yr);
if f>t
    if LIMO.design.nb_conditions == 0
        new_freq = LIMO.data.sampling_rate * ((f/t)+0.1);
        if new_freq > LIMO.data.sampling_rate
            disp('adjusting sampling rate for multivariate analysis')
            EEG.data = Yr; EEG.srate = LIMO.data.sampling_rate;
            EEG.nbchan = size(LIMO.data.expected_chanlocs,2);
            EEG.pnts = size(Yr,2); EEG.trials = size(Yr,3);
            EEG.event = []; EEG.setname = [];
            resampled_data = pop_resample(EEG, new_freq); % use EEGLab pop resample - faster if one has signal processing toolbox
        end
        Yr = resampled_data.data;
        clear f t new_freq EEG resampled_data
    else
        for condition = 1:LIMO.design.nb_conditions
            [f,t] = size(squeeze(Yr(1,:,find(LIMO.design.X(:,condition)))));
            new_freq(condition) = LIMO.data.sampling_rate * ((f/t)+0.1);
        end
        if max(new_freq) > LIMO.data.sampling_rate
            disp('adjusting sampling rate for multivariate analysis')
            EEG.data = Yr; EEG.srate = LIMO.data.sampling_rate;
            EEG.nbchan = size(LIMO.data.expected_chanlocs,2);
            EEG.pnts = size(Yr,2); EEG.trials = size(Yr,3);
            EEG.event = []; EEG.setname = [];
            resampled_data = pop_resample(EEG, max(new_freq));
        end
        Yr = resampled_data.data;
        clear f t new_freq EEG resampled_data condition
    end
end


%% Compute the weights
%---------------------
W = ones(size(Yr,1),size(Yr,3));
for e = 1:length(array)
    electrode = array(e);
    fprintf('Estimating weights at electrode %g \n',electrode);
    W(e,:) = limo_mgv(squeeze(Yr(electrode,:,:))); 
end

end % closes the function





%% lowess subfuncton

function ys = lowess(y, span)
% Smooth data using a Lowess method.

y = y(:);
n = length(y);
span = 2*floor(span/2) + 1; % make sure we have odd number
half_span = (span-1)/2; % halfspan 
x = 1:n;

% compute the smooth data
d = abs((1-half_span:half_span-1));    
dmax = half_span;  
x1 = (2:span-1)-(half_span+1);
weight = (1 - (d/dmax).^3).^1.5; % tri-cubic weight
v = [ones(length(x1),1) x1(:)];
V = v .* repmat(weight',1,size(v,2));

% Do QR decomposition and project on the middle point of QQ'
[Q,ignore] = qr(V,0);
alpha = Q(half_span,:)*Q';
alpha = alpha .* weight;
ys = filter(alpha,1,y);
ys(half_span+1:end-half_span) = ys(span-1:end-1);

% Loop over the points where the span is incomplete.  
x1 = 1:span-1;
v = [ones(length(x1),1) x1(:)];
for j=1:half_span
    d = abs((1:span-1) - j);
    weight = (1 - (d/(span-j)).^3).^1.5;
    V = v .* repmat(weight(:),1,size(v,2));
    [Q,ignore] = qr(V,0);
    alpha = Q(j,:)*Q';
    alpha = alpha .* weight;
    ys(j) = alpha * y(1:span-1);
    ys(end+1-j) = alpha * y(end:-1:end-span+2);
end
end


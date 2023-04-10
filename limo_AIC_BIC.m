function [aic, bic] = limo_AIC_BIC(varargin)

% Akaike information criterion (aic) and Bayesian information criterion (bic)
% for a GLM computed by LIMO tools
%
% FORMATS: [aic, bic] = limo_AIC_BIC(LIMO) % default
%          [aic, bic] = limo_AIC_BIC(LIMO,k)
%          [aic, bic] = limo_AIC_BIC(key,value)
%
% INPUTS  LIMO is the default LIMO.mat computed and updated after running limo_glm.m
%         k is an additional value added to the penalty term - for instance
%             is can be used to compare information resulting from the same
%             model but different pre-processing pipelines
%         key,values are pairs that can be used to call this function more
%                    flexibly, especially for low dimensional data
%          'Y', a data vector [n*1]
%          'X', the design matrix [n*p]
%          'Betas', the model parameters [p*1] (if omitted an OLS pinv(X)*Y is computed)
%          'family','gaussian' (default) or 'poisson', 'binomial' or 'none'
%
% Hector Lorenzo Mebenga & Cyril Pernet
% -------------------------------------
%  Copyright (C) LIMO Team 2023

%% defaults
k      = [];
family = 'gaussian';

%% check inputs
if nargin == 0
    help limo_AIC_BIC
    return
elseif nargin <= 2
    LIMO = varargin{1};
    if ischar(LIMO)
        LIMO = load(LIMO);
        LIMO = LIMO.LIMO;
    end

    % load the data
    Yr        = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr; % [channels x frames x trials]
    Betas     = load(fullfile(LIMO.dir,'Betas.mat')); Betas = Betas.Betas; %  [channels x frames x number of parameters]
    if strcmpi(LIMO.Analysis,'Time-Frequency')  %% add dim avoiding looking for LIMO
        Yr    = limo_tf_4d_reshape(Yr);
        Betas = limo_tf_4d_reshape(Betas);
    end
    n     = size(LIMO.design.X,1);
    p     = rank(LIMO.design.X);
    if n ~= size(Yr,3)
        error('the design matrix in LIMO and the data do not have the same number of trials')
    end
    
    X = LIMO.design.X;

    if nargin == 2
        k = varargin{2};
    end

else
    for v = 1:nargin
        if any(strcmpi(varargin{v},{'Y','Yr'}))
            Yr = varargin{v+1};
        elseif any(strcmpi(varargin{v},{'B','Betas'}))
            Betas = varargin{v+1};
        elseif any(strcmpi(varargin{v},{'X','Design'}))
            X     = varargin{v+1};

            n     =size(X,1);
            p     =rank(X);
        elseif strcmpi(varargin{v},'k')
            k = varargin{v+1};
        end
    end

    if ~exist('Yr','var')
        error('some data must be passed along');
    end

    if ~exist('X','var')
        error('The design matrix must be passed along');
    end

    if ~exist('Betas','var')
        Betas = pinv(X)*Yr;
    end

    if n ~= size(Yr,1)
        if n == size(Yr,2)
            Yr = Yr';
        else
            error('the design matrix and the data do not have the same number of trials')
        end
    end

    

end

%% compute

% if low dim data, make them high dim
if size(Yr,2) == 1
    Y = NaN(1,1,size(Yr,1));
    Y(1,1,:) = Yr; Y(1,2,:) = Yr;
    clear Yr; Yr = Y; clear Y;
end

if size(Betas,2) == 1
    B = NaN(1,1,size(Betas,1));
    B(1,1,:) = Betas; B(1,2,:) = Betas;
    clear Betas; Betas = B; clear B;
end

% run the analysis channel wise

array = find(~isnan(Yr(:,1,1)));
for channel= 1:length(array)
    residuals = squeeze(Yr(array(channel),:,:))' - X*squeeze(Betas(array(channel),:,:))';
    sigma2=sum(residuals.^2)/n;

    % compute likelihood
    if strcmpi(family,'none')
        ll =-(n/2)*log(2*pi*sigma2) - sum(residuals.^2)/(2*sigma2);
    elseif strcmpi(family,'binomial')
        ll = sum(log(binopdf(y, 1, exp(X*beta_hat)./(1+exp(X*beta_hat)))));
    elseif strcmpi(family,'poisson')
        ll = sum(log(poisspdf(squeeze(Yr(array(channel),:,:)), exp(X*beta_hat))));
    elseif strcmpi(family,'gaussian')
        ll = -(n/2)*log(2*pi*sigma2) - sum(residuals.^2)/(2*sigma2);
    end

    if ~isempty(k)
        aic = -2*ll + 2*p+sqrt(k);
        bic = -2*ll + p*log(n)+sqrt(k);
    else
        aic = -2.*ll + 2*p;
        bic = -2.*ll + p*log(n);
    end
end

if exist('LIMO','var')
    if strcmpi(LIMO.Analysis,'Time-Frequency') %% add dim avoiding looking for LIMO
        aic = limo_tf_4d_reshape(aic);
        bic = limo_tf_4d_reshape(bic);
    end
end

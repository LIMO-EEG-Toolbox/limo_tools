function result = limo_robust_rep_anova(varargin)

% result = limo_rep_anova(data,gp,factors)
%
% This function computes a robust repeated measures ANOVAs. 
% Unlike standard ANOVAs, we use a multivariate framework which accounts 
% for the correlation across measures. One advantage of this approach is 
% that one does not have to account for sphericity.
% We simply compute the T2 Hotelling test on the repeated measures.
% On top of that, means are substituted by trimmed means and variance by
% winsorized variances -- this allows to have a robust version of the test.
% The code implements equations described in:
% Rencher (2002) Methods of multivariate analysis John Wiley.
% with modifications for trimmed means and winsorized covariance
%
% INPUT
%
%   result = limo_rep_anova(data,factors)
%   result = limo_rep_anova(data,factors,C)
%   result = limo_rep_anova(data,factors,C,S)
%
% - data is a 3D matrix (f time frames * n subjects * p measures) of repeated measures
% - factors is a vector indicating the levels of each factor
%       (prod(factors)=p), e.g. [2 2] for a 2 by 2 factorial
% - C is optional and represent the contrast vector to compute (see limo_OrthogContrasts)
% - S is optional and is the covariance of the data 
%
% OUTPUT
%
%   F     = F value of the Hotelling T2 test
%   p     = corresponding level of significance
%   names = a letter per factor by default
%
%   - One sample repeated measure
%       result.F
%       result.p
% 
%   - Repeated measure with more than one factor
%       result.F
%       result.p
%       result.names
% 
% See also limo_rfep_anova, limo_random_robust, limo_OrthogContrasts
%
% Cyril Pernet & Guillaume Rousselet April 2014
% -----------------------------------------------------------------
%  Copyright (C) LIMO Team 2014


%% input stuff
% -------------
C = []; S = []; 
          
if nargin == 2
    Data = varargin{1};
    factors = varargin{2};
elseif nargin == 3 || nargin == 4
    Data = varargin{1};
    factors = varargin{2};
    C = varargin{3};
    if nargin == 4
        S = varargin{4}; % sample cov
    end
else
    error('wrong number of arguments')
end

clear varargin

%% basic info about the design
% -----------------------------
[f,n,p]       = size(Data);
nb_factors    = size(factors,2);
nb_effects    = (2^nb_factors - 1);
nb_conditions = prod(factors);

%% analyze
% ---------
type = nb_factors;

%% 
switch type
    
% One sample repeated measure
% ---------------------------
 
    % ---------------------------------------------------------------------
    case{1}  % 1 factor
        % -----------------------------------------------------------------
        if isempty(C) 
            C = [eye(p-1) ones(p-1,1).*-1];  % contrast matrix
        end
        
        if isempty(S)
            S = NaN(size(Data,1),size(Data,3),size(Data,3));
            for time = 1:size(Data,1)
                S(time,:,:) = limo_robust_cov(squeeze(Data(time,:,:))); % covariance to account for spericity
            end
        end

        df           = p-1; 
        dfe          = n-p+1; 
        y            = squeeze(limo_trimmed_mean(Data,2)); % these are the means to compare
        for time = 1:size(Data,1)
            Tsquare(time)      = n*(C*y(time,:)')'*inv(C*squeeze(S(time,:,:))*C')*(C*y(time,:)');   % Hotelling Tsquare
        end
        result.F     = ( dfe / ((n-1)*df) ) * Tsquare; 
        result.p     = 1 - fcdf(result.F, df, dfe);

        % -----------------------------------------------------------------
    case{2} % several factors
        % -----------------------------------------------------------------
        if isempty(C) 
            [C,result.names] = limo_OrthogContrasts(factors); % set of orthogonal contrasts between factors
        end
        
        if isempty(S)
            S = NaN(size(Data,1),size(Data,3),size(Data,3));
            for time = 1:size(Data,1)
                S(time,:,:) = limo_robust_cov(squeeze(Data(time,:,:))); % covariance to account for spericity
            end
        end
        
        y = squeeze(limo_trimmed_mean(Data,2)); % these are the means to compare   
        if iscell(C)
            for effect = 1:size(C,2)
                c                    = C{effect};
                df(effect)           = rank(c);
                dfe(effect)          = n-df(effect);
                for time = 1:size(Data,1)
                    Tsquare(time)    =  n*(c*y(time,:)')'*inv(c*squeeze(S(time,:,:))*c')*(c*y(time,:)'); % is also t = sqrt(n*(c*y)'*inv(c*S*c')*(c*y));
                end
                result.F(effect,:)   = ( dfe(effect) / ((n-1)*(df(effect))) ) * Tsquare;
                result.p(effect,:)   =  1 - fcdf(result.F(effect,:), df(effect), dfe(effect));
            end
        else
            df           = rank(C);
            dfe          = n-df;
            for time = 1:size(Data,1)
                Tsquare(time)    =  n*(C*y(time,:)')'*inv(C*squeeze(S(time,:,:))*C')*(C*y(time,:)');
            end
            result.F   = ( dfe / ((n-1)*(df)) ) * Tsquare;
            result.p   =  1 - fcdf(result.F, df, dfe);
        end
 end

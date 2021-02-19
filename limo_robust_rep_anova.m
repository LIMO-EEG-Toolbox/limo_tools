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
%   result = limo_rep_anova(data,gp,factors,percent)
%   result = limo_rep_anova(data,gp,factors,percent,C)
%   result = limo_rep_anova(data,gp,factors,percent,C,S)
%
% - data is a 3D matrix (f time frames * n subjects * p measures) of repeated measures
% - gp is a vector (n*1) indicating to which group subjects belong.
%       For instance enter [1 1 1 2 2 2 3 3 3] to indicate that subjects
%       1-3 belonged to group 1, subjects 4-6 belonged to group 2...
%       Enter [] if all subjects belonged to the same group.
%       There must be the same number of subjects in each group.
% - factors is a vector indicating the levels of each factor
%       (prod(factors)=p), e.g. [2 2] for a 2 by 2 factorial
% - C is optional and represent the contrast vector to compute (see limo_OrthogContrasts)
% - S is optional and is the covariance of the data (only for within factors)
% - X is optional and is the design matrix (only for within by between factors)
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
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% hard coded amount of trimming

%% input stuff
% -------------
C = []; S = []; X = []; 
          
if nargin >= 3
    Data    = varargin{1};
    gp      = varargin{2};
    factors = varargin{3};
    if nargin == 4
        percent = varargin{4};
    else
        percent = 20/100;
    end
else
    error('wrong number of arguments')
end

if nargin == 5 || nargin == 6
    C  = varargin{5};
    if nargin == 6
        if size(varargin{2},1) == size(varargin{6},1)
            X = varargin{6}; % design matrix for groups
        else
            S = varargin{6}; % sample cov if no groups
        end
    end
end

if isempty(gp)
    gp = ones(size(Data,1), 1);
end

clear varargin

%% basic info about the design
% -----------------------------
[~,n,p]       = size(Data);
nb_factors    = size(factors,2);
nb_effects    = (2^nb_factors - 1);
nb_conditions = prod(factors);
nb_gp         = size(gp,2);
if nb_gp > 1
    errordlg('Designs with more than 1 between factor are not supported')
    return
end

%% analyze
% ---------
if unique(gp) == 1
    % one sample
    if nb_factors ==1
        type = 1;
    elseif nb_factors >1
        type = 2;
    end
else
    % k samples
    if nb_factors ==1
        type = 3;
    elseif nb_factors >1
        type = 4;
    end
end


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
        
        df           = p; 
        h            = n-2*floor(percent*n); % number of items to trim
        dfe          = h-p;                
        y            = squeeze(limo_trimmed_mean(permute(Data,[1 3 2]),percent)); % these are the means to compare
        Tsquare      = NaN(1,size(Data,1));

        if isempty(S)
            for time = 1:size(Data,1)
                S = limo_robust_cov(squeeze(Data(time,:,:)),percent); % covariance to account for spericity
                Tsquare(time)  = h*(C*y(time,:)')'*inv(C*S*C')*(C*y(time,:)');   % Hotelling Tsquare
           end
        else
            for time = 1:size(Data,1)
                Tsquare(time)  = h*(C*y(time,:)')'*inv(C*squeeze(S(time,:,:))*C')*(C*y(time,:)');   % Hotelling Tsquare
            end
        end
        
        result.F     = (dfe/((n-1)*df)) * Tsquare;  
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

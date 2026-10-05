function result = limo_robust_rep_anova(varargin)

% result = limo_robust_rep_anova(data,gp,factors)
%
% This function computes a robust repeated measures ANOVA.
% Unlike standard ANOVAs, we use a multivariate framework which accounts
% for the correlation across measures. One advantage of this approach is
% that one does not have to account for sphericity.
% We simply compute the T2 Hotelling test on the repeated measures.
% On top of that, means are substituted by trimmed means and variance by
% winsorized variances -- this allows to have a robust version of the test.
% The code implements equations described in:
% Rencher (2002) Methods of multivariate analysis John Wiley.
% with modifications for trimmed means and winsorized covariance.
%
% This is the robust analogue of limo_rep_anova and MUST share its calling
% convention: it is called from limo_random_robust and limo_contrast with
% the group vector gp as the 2nd argument (see those callers).
%
% INPUT
%
%   result = limo_robust_rep_anova(data,gp,factors)
%   result = limo_robust_rep_anova(data,gp,factors,C)
%   result = limo_robust_rep_anova(data,gp,factors,C,S)   % S = robust cov (within only)
%
% - data is a 3D matrix (f time/freq frames * n subjects * p measures)
% - gp is a vector (n*1) indicating to which group subjects belong.
%       Enter [] (or a constant vector) if all subjects belong to the same group.
% - factors is a vector indicating the levels of each factor
%       (prod(factors)=p), e.g. [2 2] for a 2 by 2 factorial
% - C is optional and represents the contrast(s) to compute (see limo_OrthogContrasts)
% - S is optional and is the (robust) covariance of the data (within factors only)
%
% OUTPUT
%
%   - One sample repeated measure (within only, 1 factor)
%       result.F  result.df  result.dfe  result.p
%
%   - Repeated measure with more than one factor (within only)
%       result.F  result.df  result.dfe  result.p  result.names
%
% NOTE on between-group designs: a *robust* within-by-between repeated
% measures ANOVA (i.e. gp defines more than one group) is not implemented.
% The pooled between-group covariance required for a trimmed-mean Hotelling
% test is not a drop-in substitution of limo_robust_cov and needs a
% dedicated, validated derivation. Such designs error out here rather than
% return an unvalidated statistic; use limo_rep_anova (non-robust) instead.
%
% See also limo_rep_anova, limo_random_robust, limo_OrthogContrasts,
%          limo_trimmed_mean, limo_robust_cov
%
% Cyril Pernet & Guillaume Rousselet April 2014
% ------------------------------
%  Copyright (C) LIMO Team 2019

percent = 20; % percentage of trimming for trimmed means / winsorized cov

%% input stuff
% -------------
C = []; S = [];

if nargin == 3
    Data    = varargin{1};
    gp      = varargin{2};
    factors = varargin{3};
elseif nargin == 4 || nargin == 5
    Data    = varargin{1};
    gp      = varargin{2};
    factors = varargin{3};
    C       = varargin{4};
    if nargin == 5
        % as in limo_rep_anova: a 5th argument whose number of rows matches
        % gp is the between-group design matrix X; otherwise it is the
        % within covariance S.
        if size(varargin{2},1) == size(varargin{5},1)
            X = varargin{5}; %#ok<NASGU> % design matrix for groups (between)
        else
            S = varargin{5}; % sample covariance if no groups
        end
    end
else
    error('wrong number of arguments')
end

if isempty(gp)
    gp = ones(size(Data,2), 1);
end

clear varargin

%% basic info about the design
% -----------------------------
[f,n,p]    = size(Data); %#ok<ASGLU> % frames * subjects * measures
nb_factors = size(factors,2);

% dispatch on within (single group) vs between (>1 group)
if length(unique(gp)) == 1
    if nb_factors == 1
        type = 1;
    else
        type = 2;
    end
else
    if nb_factors == 1
        type = 3;
    else
        type = 4;
    end
end

%%
switch type

% One sample repeated measure (within factors only)
% -------------------------------------------------

    % ---------------------------------------------------------------------
    case 1  % 1 factor
        % -----------------------------------------------------------------
        if isempty(C)
            C = [eye(p-1) ones(p-1,1).*-1];  % contrast matrix
        end

        if isempty(S)
            S = NaN(size(Data,1),size(Data,3),size(Data,3));
            for frame = 1:size(Data,1)
                S(frame,:,:) = limo_robust_cov(squeeze(Data(frame,:,:))); % winsorized covariance
            end
        end

        df  = p-1;
        dfe = n-p+1;
        % trimmed mean ACROSS SUBJECTS (dim 2) -> f x p ; limo_trimmed_mean
        % reduces the 3rd dim, so bring subjects to the 3rd dim first.
        y   = limo_trimmed_mean(permute(Data,[1 3 2]),percent); % f x p (means to compare)
        for frame = 1:size(Data,1)
            Tsquare(frame) = n*(C*y(frame,:)')'*inv(C*squeeze(S(frame,:,:))*C')*(C*y(frame,:)'); %#ok<*MINV> % Hotelling Tsquare
        end
        result.F   = ( dfe / ((n-1)*df) ) * Tsquare;
        result.df  = df;
        result.dfe = dfe;
        result.p   = 1 - fcdf(result.F, df, dfe);

        % -----------------------------------------------------------------
    case 2 % several factors
        % -----------------------------------------------------------------
        if isempty(C)
            [C,result.names] = limo_OrthogContrasts(factors); % orthogonal contrasts between factors
        end

        if isempty(S)
            S = NaN(size(Data,1),size(Data,3),size(Data,3));
            for frame = 1:size(Data,1)
                S(frame,:,:) = limo_robust_cov(squeeze(Data(frame,:,:))); % winsorized covariance
            end
        end

        y = limo_trimmed_mean(permute(Data,[1 3 2]),percent); % f x p (means to compare)
        if iscell(C)
            for effect = 1:size(C,2)
                c                    = C{effect};
                df                   = rank(c);
                dfe                  = n-df;
                for frame = 1:size(Data,1)
                    Tsquare(frame)   =  n*(c*y(frame,:)')'*pinv(c*squeeze(S(frame,:,:))*c')*(c*y(frame,:)');
                end
                result.F(effect,:)   = ( dfe / ((n-1)*(df)) ) * Tsquare;
                result.p(effect,:)   =  1 - fcdf(result.F(effect,:), df, dfe);
                result.df(effect)    = df;
                result.dfe(effect)   = dfe;
            end
        else
            df  = rank(C);
            dfe = n-df;
            for frame = 1:size(Data,1)
                Tsquare(frame) =  n*(C*y(frame,:)')'*pinv(C*squeeze(S(frame,:,:))*C')*(C*y(frame,:)');
            end
            result.F   = ( dfe / ((n-1)*(df)) ) * Tsquare;
            result.p   =  1 - fcdf(result.F, df, dfe);
            result.df  = df;
            result.dfe = dfe;
        end

% k samples repeated measure (within-by-between)
% ----------------------------------------------
    case {3,4}
        error(['limo_robust_rep_anova:betweenNotSupported : robust repeated ' ...
            'measures ANOVA with a between-group factor is not implemented. ' ...
            'A validated trimmed-mean pooled-covariance derivation is required; ' ...
            'use limo_rep_anova (non-robust) for within-by-between designs.'])
end
end

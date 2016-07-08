function result = limo_rep_anova(varargin)

% result = limo_rep_anova(data,gp,factors)
%
% This function computes repeated measures ANOVAs. Unlike standard ANOVAs, 
% we use a multivariate framework which accounts for the correlation
% across measures. One advantage of this approach is that one does not 
% have to account for sphericity and thus it saves a lot of computational time.
% In short we simply either run a T2 on the repeated measures or a MANOVA 
% (generalized T2 on transformed data). The code implements equations described 
% in:bRencher (2002) Methods of multivariate analysis John Wiley.
%
% INPUT
%
%   result = limo_rep_anova(data,gp,factors)
%   result = limo_rep_anova(data,gp,factors,C,S)
%   result = limo_rep_anova(data,gp,factors,C,X)
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
%   - 1 within and 1 between factors
%       result.repeated_measure.F
%       result.repeated_measure.p
%       result.gp.F
%       result.gp.p
%       result.interaction.F
%       result.interaction.p
% 
%   - repeated measure with more than one factor and 1 between factor
%       result.repeated_measure.names
%       result.repeated_measure.F
%       result.repeated_measure.p
%       result.gp.F, result.gp.p
%       result.interaction.F
%       result.interaction.p
%       result.interaction.names
%
%
% EXAMPLE: 2 x 2 within subject design with 2 groups of subjects
%          20 times frames - here every frame contains the same data
%         data        = zeros(20,11,4);
%         data(:,:,1) = repmat([1 2 3 8 5 2 4 8 7 2 4],20,1); % A1B1
%         data(:,:,2) = repmat([1 5 7 9 4 2 6 4 8 3 5],20,1); % A1B2
%         data(:,:,3) = repmat([4 6 8 4 1 8 0 1 4 6 2],20,1); % A2B1 
%         data(:,:,4) = repmat([8 9 6 8 5 8 2 7 5 9 4],20,1); % A2B2
%         gp          = [1 1 1 1 1 1 2 2 2 2 2]'; % 2 gps with 6 and 5 subjects
%         result      = limo_rep_anova(data,gp,[2 2]);
%
%      Data organization:
%
%      GP    factor A     level 1        level 2
%            factor B level 1 level 2 level 1 level 2 
%      1                1       1       4       8
%      1                2       5       6       9
%      1                3       7       8       6
%      1                8       9       4       8
%      1                5       4       1       5
%      1                2       2       8       8
%      2                4       6       0       2
%      2                8       4       1       7
%      2                7       8       4       5
%      2                2       3       6       9
%      2                4       5       2       4
%
% See also limo_random_robust, limo_OrthogContrasts
%
% Cyril Pernet Septembre 2010
% Cyril Pernet May 2011 - update to run locally the mutlivariate 
%            analysis + added the time dimension to speed things up
% GAR April 2013: updated help & comments; added name output
% -----------------------------------------------------------------
%  Copyright (C) LIMO Team 2010


%% input stuff
% -------------
C = []; S = []; X = []; 
          
if nargin == 3
    Data = varargin{1};
    gp   = varargin{2};
    factors = varargin{3};
elseif nargin == 4 || nargin == 5
    Data = varargin{1};
    gp   = varargin{2};
    factors = varargin{3};
    C = varargin{4};
    if nargin == 5
        if size(varargin{2},1) == size(varargin{5},1)
            X = varargin{5}; % design matrix for groups
        else
            S = varargin{5}; % sample cov if no groups
        end
    end
else
    error('wrong number of arguments')
end

if isempty(gp)
    gp = ones(size(Data,1), 1);
end

clear varargin

%% basic info about the design
% -----------------------------
[f,n,p]       = size(Data);
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
        
        if isempty(S)
            S = NaN(size(Data,1),size(Data,3),size(Data,3));
            for time = 1:size(Data,1)
                S(time,:,:) = cov(squeeze(Data(time,:,:))); % covariance to account for spericity
            end
        end

        df           = p-1; 
        dfe          = n-p+1; 
        y            = squeeze(nanmean(Data,2)); % these are the means to compare
        if size(Data,1) == 1 %% no time or freq dim
            Tsquare = n*(C*y)'*inv(C*squeeze(S(1,:,:))*C')*(C*y);   % Hotelling Tsquare
        else
            for time = 1:size(Data,1)
                Tsquare(time)      = n*(C*y(time,:)')'*inv(C*squeeze(S(time,:,:))*C')*(C*y(time,:)');   % Hotelling Tsquare
            end
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
                S(time,:,:) = cov(squeeze(Data(time,:,:))); % covariance to account for spericity
            end
        end
        
        y = squeeze(nanmean(Data,2)); % these are the means to compare   
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
        
        
% k samples repeated measure
% ---------------------------

        % -----------------------------------------------------------------
    case{3} % 1 within and 1 between factors
        % -----------------------------------------------------------------
        
        % deal with group structure
        % -------------------------
        gp_values = unique(gp);
        k         = length(gp_values);
        
        % build the design matrix of groups
        if isempty(X)            
            X = NaN(size(gp,1),k+1);
            for g =1:k
                X(:,g) = gp == gp_values(g);
            end
            X(:,end) = 1;
        end
                
        % get the error matrix of the full model
        % ---------------------------------------------------
        R  = eye(size(Data,2)) - (X*pinv(X));
        for time=1:f
            E(time,:,:)  = (squeeze(Data(time,:,:))'*R*squeeze(Data(time,:,:)));                                                                              
        end
        
        % compute the repeated measure
        % ------------------------------
        if isempty(C) 
            C = [eye(p-1) ones(p-1,1).*-1]; % set of orthogonal contrasts between factors
        end

        df  = p-1; 
        dfe = (n-p+1) - (k-1); % remove from dfe nb_gp - 1
        yp  = squeeze(nanmean(Data,2))'; % average across gp
        ve  = sum(sum(X(:,1:end-1))-1); % - rank(X); % dfe for different sample sizes (gives the same as rank(X)*(sum(X(:,1))-1) for equal sample sizes            
        Spl = E/ve; % covariance of data split per gp
        for time = 1:f
            Tsquare(time)                 = n*(C*yp(:,time))'*inv(C*squeeze(Spl(time,:,:))*C')*(C*yp(:,time));
        end
        result.repeated_measure.F    = ( dfe / (ve*df) ) * Tsquare; 
        result.repeated_measure.p    = 1 - fcdf(result.repeated_measure.F, df, dfe);

        
        % compute the gp effect (=univariate stat on the mean across repeated measures)
        % ----------------------------------------------------------------------------              
        Y  = nanmean(Data,3); % average repeated measures
        [result.gp.F,result.gp.p] = local_glm(Y',X,k,sum(X(:,1:k),1),1);

        
       % compute interaction (= multivariate on differences)
       % -------------------------------------------------------
       for time=1:f
           I = (C*squeeze(Data(time,:,:))')';
           [result.interaction.F(time), result.interaction.p(time)]= local_glm(I,X,k,sum(X(:,1:k),1),2);
       end
       
        % -----------------------------------------------------------------
    case{4} % several within and 1 between factors
        % -----------------------------------------------------------------
        
        % deal with group structure
        % -------------------------
        gp_values = unique(gp);
        k = length(gp_values);
        
        % build the design matrix of gps
        if isempty(X)            
            X = NaN(size(gp,1),k+1);
            for g =1:k
                X(:,g) = gp == gp_values(g);
            end
            X(:,end) = 1;
        end
                
        % get the effect and error matrices of the full model
        % ---------------------------------------------------
        R  = eye(size(Data,2)) - (X*pinv(X));
        for time=1:f
            E(time,:,:)  = (squeeze(Data(time,:,:))'*R*squeeze(Data(time,:,:)));                                                                              
        end
        
        % compute the repeated measure
        % ------------------------------
        if isempty(C) 
            [C,result.names] = limo_OrthogContrasts(factors); % set of orthogonal contrasts between factors
        end

        y  = squeeze(nanmean(Data,2))'; % average across gp
        ve = 0; % dfe as a function of the number of subjects per gp
        for g=1:k
            v = sum(sum(X(:,g)==1));
            ve = ve + (v-1);
        end
                
        if iscell(C)
            for effect = 1:length(C)
                c   = C{effect};
                df  = rank(c);
                dfe = n-df-(k-1);
                Spl = E/ve;
                for time = 1:f
                    Tsquare(time) = n*(c*y(:,time))'*inv(c*squeeze(Spl(time,:,:))*c')*(c*y(:,time));
                end
                result.repeated_measure.F(effect,:) = ( dfe / (ve*df) ) .* Tsquare;
                result.repeated_measure.p(effect,:) = 1 - fcdf(result.repeated_measure.F(effect,:), df, dfe);
            end
        else
            df  = rank(C);
            dfe = n-df-(k-1);
            Spl = E/ve;
            for time = 1:f
                Tsquare(time) = n*(C*y(:,time))'*inv(C*squeeze(Spl(time,:,:))*C')*(C*y(:,time));
            end
            result.repeated_measure.F = ( dfe / (ve*df) ) .* Tsquare;
            result.repeated_measure.p = 1 - fcdf(result.repeated_measure.F, df, dfe);
        end
        
        % compute the gp effect (=univariate stat)
        % ---------------------------------------     
        Y  = nanmean(Data,3); % average repeated measures
        [result.gp.F, result.gp.p] = local_glm(Y',X,k,sum(X(:,1:k),1),1);
       
       % compute the interactions with gp
       % ---------------------------------
       if iscell(C)
           for effect = 1:length(C)
               c = C{effect};
               for time=1:f
                   I = (c*squeeze(Data(time,:,:))')';
                   [result.interaction.F(effect,time), result.interaction.p(effect,time)]= local_glm(I,X,k,sum(X(:,1:k),1),2);
               end
           end
       else
           for time=1:f
               I = (C*squeeze(Data(time,:,:))')';
               [result.interaction.F(time), result.interaction.p(time)]= local_glm(I,X,k,sum(X(:,1:k),1),2);
           end
       end
       
% --
end % closes the switch
end % closes the function 

% subfuncton
% --------------
function [F,p] = local_glm(Y,X,nb_gp,nb_subjects,flag)

T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
R        = eye(size(Y,1)) - (X*pinv(X));  % Residual matrix
E        = (Y'*R*Y);   % SS Error
Betas    = pinv(X)*Y;    
C = eye(size(X,2));
C(:,size(X,2)) = 0;
C0   = eye(size(X,2)) - C*pinv(C);
X0   = X*C0;                                                              
R0   = eye(size(Y,1)) - (X0*pinv(X0));
M    = R0 - R;    % M is the projection matrix onto Xc
H    = (Betas'*X'*M*X*Betas);    % SS Hypothesis (Effect)

if flag == 2
    [Eigen_vector, Eigen_values] = limo_decomp(E,H);
    p = size(Y,2); % = number of variables (dimension)
    vh = nb_gp - 1; % df = q above
    s = min(vh,p); % subspace in which mean Ys are located
    if sum(nb_subjects == nb_subjects(1)) == length(nb_subjects)
        ve = nb_gp*(nb_subjects(1)-1);     % dfe equal sample sizes
    else
        ve = sum(nb_subjects) - nb_gp;     % dfe different sample sizes
    end
    
    if s > 1
        m = (abs(vh-p)-1)/2;
        N = (ve-p-1) / 2;
        U = sum(Eigen_values);
        df_Hotelling = s*(2*m+s+1);
        dfe_Hotelling = 2*(s*N+1);
        F = (dfe_Hotelling*U) / (s^2*(2*m+s+1));
        p = 1-fcdf(F, df_Hotelling, dfe_Hotelling);
    else % = only one non zeros Eigen value s = 1 and/or vh = 1
        U = max(Eigen_values);
        df_Hotelling = p;
        dfe_Hotelling = ve-p+1;
        F = (dfe_Hotelling/df_Hotelling) * max(Eigen_values);
        p = 1-fcdf(F,df_Hotelling,dfe_Hotelling);
    end
else
    F    = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
    p    = 1 - fcdf(F, (rank(X)-1), (size(Y,1)-rank(X)));
end
end

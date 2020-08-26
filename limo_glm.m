function model = limo_glm(varargin)

% LIMO_GLM - core function of the LIMO toolbox
% Computes the basic statistics of the models using pinv
% Analyses are performed at one electrode only, but for all trials and all
% time points.
%
% FORMAT:
% model = limo_glm(Y,LIMO)
% model = limo_glm(Y,X,nb_conditions,nb_continuous,nb_items,Method)
%
% INPUTS:
%   Y             = 2D matrix of EEG data with format trials x frames 
%   X             = 2 dimensional design matrix - see LIMO_DESIGN_MATRIX
%   nb_conditions = number of conditions / groups - see LIMO_DESIGN_MATRIX
%   nb_continuous = number of covariates - see LIMO_DESIGN_MATRIX
%   nb_items      = number of trials per conditions - see LIMO_DESIGN_MATRIX
%   Method        = 1 for massive univariate or 2 for multivariate
%   LIMO          = structure that contains the above information (except Y)
%
% OUTPUTS:
%   model         = structure of results
%
%       model.R2_univariate
%       model.F
%       model.p
%       model.betas - format betas x frames
%       model.df
%
%       model.univariate - univariate results
%           model.univariate.conditions = categorical effects
%           model.univariate.continuous = continuous effects
%
%       model.multivariate - multivariate results:
%           model.R2_multivariate                       
%           model.multivariate.EV - Eigen Vectors
%
% NOTES:
%
% The adjusted data are not computed here but the effects of the
% covariates are regressed out during the matrix inversion, i.e. the 
% beta parameters reflect the unique part of variance of each regressor
% (i.e. if two regressors are correlated, the common part of variance
% goes into the error). The values of the beta parameters are not always the true
% values because one uses a pseudoinverse to invert X. However, the linear combination
% of vectors in X is optimal (gives the same Yhat as using an inverse) and their ratios
% give the same statistics. The multivariate statistic uses Pillais' V and Roy's
% test. In Roy's test one maximizes the spread of the transformed data using the
% 1st eigen value of inv(E)*H. The F approximation is an upper bound, i.e.
% results are safe if H0 is accepted (no effect) but not rejected. If data
% are highly correlated, there is only one high eigen value and Roy test is
% appropriate otherwise Pillai is better. We also comnpute the generalized
% Hotelling test, which is used only by limo_rep_anova on transformed data
%
% References
% Yandell, B.S. 1997. Practical Data Analysis For Designed Experiments. Chapman & Hall
% Friston et al. 2007. Statitical Parametric Mapping. Academic Press
% Rencher, A.C. 2002. Methods of Multivariate Analysis. 2nd Ed. Wiley.
%
% See also LIMO_DECOMP LIMO_DESIGN_MATRIX LIMO_REP_ANOVA
%
% Cyril Pernet 13-02-2009
% -----------------------------
%  Copyright (C) LIMO Team 2010

%% varagin

if length(varargin) == 2
    Y             = varargin{1};
    X             = varargin{2}.design.X;
    nb_conditions = varargin{2}.design.nb_conditions;
    nb_continuous = varargin{2}.design.nb_continuous;
    nb_items      = varargin{2}.design.nb_items;
    Method        = varargin{2}.Method;
elseif length(varargin) == 6
    Y             = varargin{1};
    X             = varargin{2};
    nb_conditions = varargin{3};
    nb_continuous = varargin{4};
    nb_items      = varargin{5};
    Method        = varargin{6};
else
    error('varargin error')
end


%% Data check

if size(Y,1)~=size(X,1)
    error('The number of events in Y and the design matrix are different, error line 86 in limo_glm')
end

if  nb_conditions ~=0
    if  size(X,1) ~= sum(nb_items)
        error('The number of items and the design matrix are different, error line 91 in limo_glm')
    end
end

%% Compute

T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
R        = eye(size(Y,1)) - (X*pinv(X));                                      % Residual matrix
E        = (Y'*R*Y);                                                          % SS Error
Betas    = pinv(X)*Y;                                                         % Parameters of the model


% --------------------
% compute model R^2
% --------------------
C = eye(size(X,2));
C(:,size(X,2)) = 0;
C0   = eye(size(X,2)) - C*pinv(C);
X0   = X*C0;                                                                  % Reduced (residual) model
R0   = eye(size(Y,1)) - (X0*pinv(X0));
M    = R0 - R;                                                                % M is the projection matrix onto Xc
H    = (Betas'*X'*M*X*Betas);                                                 % SS Hypothesis (Effect)
Rsquare   = diag(H)./diag(T);                                                 % Variances explained per Y
F_Rsquare = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
p_Rsquare = 1 - fcdf(F_Rsquare, (rank(X)-1), (size(Y,1)-rank(X)));
if Method == 2
    Eigen_values_R2 = limo_decomp(E,H);
end
  
% -----------------------------------
% compute F for categorical variables
% -----------------------------------
if nb_conditions ~= 0 && nb_continuous == 0
    F_conditions    = F_Rsquare;
    pval_conditions = p_Rsquare;
    if Method == 2
        Eigen_values_cond = Eigen_values_R2;
    end

elseif nb_conditions ~= 0 && nb_continuous ~= 0
    C = eye(size(X,2));
    C(:,(nb_conditions+1):size(X,2)) = 0;
    C0   = eye(size(X,2)) - C*pinv(C);
    X0   = X*C0;
    R0   = eye(size(Y,1)) - (X0*pinv(X0));
    M    = R0 - R;
    H    = (Betas'*X'*M*X*Betas);
    df_conditions = rank(C)-1;
    F_conditions    = (diag(H)/(rank(C)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
    pval_conditions = 1 - fcdf(F_conditions(:), rank(C)-1, (size(Y,1)-rank(X)));
    if Method == 2
        Eigen_values_cond = limo_decomp(E,H);
    end
end


% ----------------------------------
% compute F for continuous variables
% ----------------------------------
if nb_continuous ~= 0
    if nb_continuous == 1 && nb_conditions == 0
        F_continuous    = F_Rsquare;
        pval_continuous = p_Rsquare;
        if Method == 2
            Eigen_values_continuous = limo_decomp(E,H);
        end

    else

        for n = 1:nb_continuous
            C    = zeros(size(X,2));
            C(nb_conditions+n,nb_conditions+n) = 1;
            C0   = eye(size(X,2)) - C*pinv(C);
            X0   = X*C0;
            R0   = eye(size(Y,1)) - (X0*pinv(X0));
            M    = R0 - R;
            H    = Betas'*X'*M*X*Betas;
            for m = 1:size(Y,2)
                F_continuous(n,m)    = ((Betas(:,m)'*X'*M*X*Betas(:,m))/(rank(C)))/((Y(:,m)'*R*Y(:,m))/(size(Y,1)-rank(X)));
                pval_continuous(n,m) = 1 - fcdf(F_continuous(n,m), 1, (size(Y,1)-rank(X)));
            end

            if Method == 2
                Eigen_values_continuous(:,n) = limo_decomp(E,H);
            end
        end
                       
    end
end

 

% ---------------------
% Multivariate analyses
% ---------------------


if Method == 2

    % -------
    % R2
    % -------
    
    p = size(Y,2); % = number of variables (dimension)
    q = size(X,2)-1; % = number of regressors (df)
    n = size(Y,1); % nb of observations (dfe)
    s = min(p,q);
    m = (abs(q-p)-1)/2;
    N = (n-q-p-2)/2;
    d = max(p,q);
    
    theta = max(Eigen_values_R2) / (1+max(Eigen_values_R2)); % Roy
    V = sum(Eigen_values_R2 ./ (1+Eigen_values_R2)); % Pillai
    
    R2_Roy_value = theta; % = 1st canonical correlation
    R2_Roy_F     = ((n-d-1)*max(Eigen_values_R2))/d;
    R2_Roy_p     = 1-fcdf(R2_Roy_F, d, (n-d-1));
    
    R2_Pillai_value = V / s; % average of canonical correlations
    R2_Pillai_F     = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
    R2_Pillai_p     = 1-fcdf(R2_Pillai_F,(s*(2*m+s+1)),(s*(2*N+s+1)));

    %   % Generalized R2
    %   % variance covariance matrix
    %   S = cov([Y X(:,1:size(X,2)-1)]);
    %   Syy = S(1:size(Y,2),1:size(Y,2));
    %   Sxy = S(size(Y,2)+1:size(S,1),1:size(Y,2));
    %   Syx = S(1:size(Y,2),size(Y,2)+1:size(S,2));
    %   Sxx = S(size(Y,2)+1:size(S,1),size(Y,2)+1:size(S,2));
    %   Rsquare_multi = trace(Sxy*Syx) / sqrt(trace(Sxx.^2)*trace(Syy.^2)); % Robert and Escoufier, J.Royal Stat Soc, C - 1976
    %   Rsquare_multi = det(Syx*inv(Sxx)*Sxy) / det(Syy); % same as univariate but uses the determinant

    
    % ------
    % MANOVA
    % ------
    if nb_conditions ~= 0 && nb_continuous == 0 
        
        vh = nb_conditions - 1; % df = q above
        s = min(vh,p); % subspace in which mean Ys are located
        if sum(nb_items == nb_items(1)) == length(nb_items)
            ve = nb_conditions*(nb_items(1)-1);     % dfe equal sample sizes
        else
            ve = sum(nb_items) - nb_conditions;     % dfe different sample sizes
        end
        
        if s > 1

            m = (abs(vh-p)-1)/2;
            N = (ve-p-1) / 2;
            
            % Pillai
            V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
            df_Pillai = s*(2*m+s+1);
            dfe_Pillai = s*(2*N+s+1);
            F_conditions_Pillai = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
            pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_Pillai,dfe_Pillai);

            % Roy's test
            theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond));
            df_Roy = max(p,vh);
            dfe_Roy = ve - 1; % in Renchner it is proposed to use ve - max(p,vh) -1 while in Statistica it is ve -1
            F_conditions_Roy = (dfe_Roy*max(Eigen_values_cond))/df_Roy;
            pval_conditions_Roy = 1-fcdf(F_conditions_Roy, df_Roy, dfe_Roy);
            
            % Lawley-Hotelling's generalized T2
            U = sum(Eigen_values_cond);
            df_Hotelling = s*(2*m+s+1); dfe_Hotelling = 2*(s*N+1);
            F_conditions_Hotelling = (dfe_Hotelling*U) / (s^2*(2*m+s+1));
            pval_conditions_Hotelling = 1-fcdf(F_conditions_Hotelling, df_Hotelling, dfe_Hotelling);
          

        else % = only one non zeros Eigen value s = 1 and/or vh = 1
            
            V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
            U = max(Eigen_values_cond);
            theta = U;
            
            df_Pillai = p; % number of frames
            dfe_Pillai = ve-p+1; 
            df_Roy = df_Pillai; dfe_Roy = dfe_Pillai;
            df_Hotelling = df_Roy; dfe_Hotelling = dfe_Roy;
            
            F_conditions_Pillai = (dfe_Pillai/df_Pillai) * max(Eigen_values_cond);
            pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_Pillai,dfe_Pillai);
            F_conditions_Roy = F_conditions_Pillai;
            pval_conditions_Roy = pval_conditions_Pillai;
            F_conditions_Hotelling = F_conditions_Pillai;
            pval_conditions_Hotelling = pval_conditions_Pillai;


        end
    end % closes MANOVA case

     
    % --------------------------
    % Regressions and/or MAnCOVA
    % --------------------------
    if nb_continuous ~= 0 
                
        p = size(Y,2);
        vh = nb_conditions - 1;
        s = min(vh,p);

        % MAnCOVA
        % -------
        if nb_conditions ~= 0

            if sum(nb_items == nb_items(1)) == length(nb_items)
                ve = nb_conditions*(nb_items(1)-1);
            else
                ve = sum(nb_items) - nb_conditions;
            end

            if s > 1

                m = (abs(vh-p)-1)/2;
                N = (ve-p-1) / 2;
                
                % Pillai
                V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
                df_Pillai = s*(2*m+s+1);
                dfe_Pillai = s*((2*N+s+1)-nb_continuous); 
                F_conditions_Pillai = (((2*N+s+1)-nb_continuous)*V) / ((2*m+s+1)*(s-V)');
                pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_Pillai,dfe_Pillai);

                % Roy's test
                theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond));
                df_Roy = max(p,vh);
                dfe_Roy = ve - 1 - nb_continuous;
                F_conditions_Roy = (dfe_Roy*max(Eigen_values_cond))/df_Roy;
                pval_conditions_Roy = 1-fcdf(F_conditions_Roy, df_Roy, dfe_Roy);
                
                % Lawley-Hotelling's generalized T2
                U = sum(Eigen_values_cond);
                df_Hotelling = s*(2*m+s+1); dfe_Hotelling = 2*(s*N+1);
                F_conditions_Hotelling = (dfe_Hotelling*U) / (s^2*(2*m+s+1));
                pval_conditions_Hotelling = 1-fcdf(F_conditions_Hotelling, df_Hotelling, dfe_Hotelling);

            else % = only one non zeros Eigen value s = 1 and/or vh = 1
                
                V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
                U = max(Eigen_values_cond);
                theta = U;
                
                df_Pillai = p; df_Pillai = p; % number of frames
                dfe_Pillai = (ve-p+1) - nb_continuous;
                df_Roy = df_Pillai; dfe_Roy = dfe_Pillai;
                df_Hotelling = df_Roy; dfe_Hotelling = dfe_Roy;
                
                F_conditions_Pillai = (dfe_Pillai/df_Pillai) * max(Eigen_values_cond);
                pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_Pillai,dfe_Pillai);
                F_conditions_Roy = F_conditions_Pillai;
                pval_conditions_Roy = pval_conditions_Pillai;
                F_conditions_Hotelling = F_conditions_Pillai;
                pval_conditions_Hotelling = pval_conditions_Pillai;
            
            end
        end
        
        % ---------------------
        % Continous regressors 
        % --------------------
       
        df_continuous = size(Y,2);
        dfe_continuous = size(Y,1)-nb_continuous-size(Y,2);
        if nb_conditions ~= 0
            dfe_continuous = dfe_continuous - 1;
        end

        s = min(size(Y,2),size(X,2));
        m = (abs(size(X,2)-size(Y,2))-1) / 2;

        for i=1:nb_continuous

            % Roy's test
            theta_continuous(i) = max(Eigen_values_continuous(:,i)) / (1+max(Eigen_values_continuous(:,i)));
            F_continuous_Roy(i) = (dfe_continuous*max(Eigen_values_continuous(:,i)))/df_continuous;
            pval_continuous_Roy(i) = 1-fcdf(F_continuous_Roy(i), df_continuous, dfe_continuous);

            % Pillai
            V_continuous(i) = sum(Eigen_values_continuous(:,i) ./ (1+Eigen_values_continuous(:,i)));
            F_continuous_Pillai(i) = F_continuous_Roy(i);
            pval_continuous_Pillai(i) = pval_continuous_Roy(i);
            
            % Lawley-Hotelling's generalized T2
            U_continuous(i) = sum(Eigen_values_continuous(:,i));
            F_continuous_Hotelling(i) = F_continuous_Roy(i);
            pval_continuous_Hotelling(i) = pval_continuous_Roy(i);
            
            % only for multivariate multiple regressions
            if nb_conditions == 0
                R2_Roy(i)    = theta_continuous(i);
                R2_Pillai(i) = V_continuous(i) / s;
                R2_F_Pillai(i) = (R2_Pillai(i) / p-1) / ((1-R2_Pillai(i))/ size(Y,1)-p-1);
                R2_pval_Pillai(i) = 1-fcdf(R2_F_Pillai(i),p-1,size(Y,1)-p-1);
            end
            
        end      
    end % closes if nb_continuous ~= 0 
end % closes if Method == 2


%% Save into a structure

% Overall model
model.R2_univariate   = Rsquare;
model.F               = F_Rsquare;
model.p               = p_Rsquare;
model.betas           = Betas;
model.df              = [rank(X)-1 (size(Y,1)-rank(X))];

% univariate results
if nb_conditions > 0 && nb_continuous == 0
    model.univariate.conditions.F             = F_conditions;
    model.univariate.conditions.p             = pval_conditions;
    model.univariate.conditions.df            = [rank(X)-1 (size(Y,1)-rank(X))];

elseif nb_conditions == 0 && nb_continuous > 0
    model.univariate.continuous.F             = F_continuous;
    model.univariate.continuous.p             = pval_continuous;
    model.univariate.continuous.df            = [1 (size(Y,1)-rank(X))];

elseif nb_conditions > 0 && nb_continuous > 0
    model.univariate.conditions.F             = F_conditions;
    model.univariate.conditions.p             = pval_conditions;
    model.univariate.conditions.df            = [df_conditions; (size(Y,1)-rank(X))];
    model.univariate.continuous.F             = F_continuous;
    model.univariate.continuous.p             = pval_continuous;
    model.univariate.continuous.df            = [1 (size(Y,1)-rank(X))];
end


if Method == 2

    model.multivariate.R2.Roy.V = R2_Roy_value;
    model.multivariate.R2.Roy.F = R2_Roy_F;
    model.multivariate.R2.Roy.p = R2_Roy_p;   
    model.multivariate.R2.pillai.V = R2_Pillai_value;
    model.multivariate.R2.pillai.F = R2_Pillai_F;
    model.multivariate.R2.pillai.p = R2_Pillai_p;
    
    
    if nb_conditions > 0 && nb_continuous == 0
        
        model.multivariate.EV                      = [Eigen_values_R2'; Eigen_values_cond'];
        model.multivariate.conditions.Pillai.V     = V;
        model.multivariate.conditions.Pillai.F     = F_conditions_Pillai;
        model.multivariate.conditions.Pillai.p     = pval_conditions_Pillai;
        model.multivariate.conditions.Pillai.df    = df_Pillai;
        model.multivariate.conditions.Pillai.dfe   = dfe_Pillai;
        model.multivariate.conditions.Pillai.R2.R2 = R2_Pillai_value;
        model.multivariate.conditions.Pillai.R2.F  = R2_Pillai_F;
        model.multivariate.conditions.Pillai.R2.p  = R2_Pillai_p;
        model.multivariate.conditions.Roy.theta    = theta;
        model.multivariate.conditions.Roy.F        = F_conditions_Roy;
        model.multivariate.conditions.Roy.p        = pval_conditions_Roy;
        model.multivariate.conditions.Roy.df       = df_Roy;
        model.multivariate.conditions.Roy.dfe      = dfe_Roy;
        model.multivariate.conditions.Roy.R2.R2    = R2_Roy_value;
        model.multivariate.conditions.Roy.R2.F     = F_conditions_Roy;
        model.multivariate.conditions.Roy.R2.p     = pval_conditions_Roy;
        model.multivariate.conditions.Hotelling.U  = U;
        model.multivariate.conditions.Hotelling.F  = F_conditions_Hotelling;
        model.multivariate.conditions.Hotelling.p  = pval_conditions_Hotelling;
        model.multivariate.conditions.Hotelling.df = df_Hotelling;
        model.multivariate.conditions.Hotelling.dfe= dfe_Hotelling;

    elseif nb_conditions == 0 && nb_continuous > 0
        
        model.multivariate.EV                      = [Eigen_values_R2'; Eigen_values_continuous'];
        model.multivariate.continuous.df           = df_continuous;
        model.multivariate.continuous.dfe          = dfe_continuous;

        model.multivariate.continuous.Pillai.V     = V_continuous;
        model.multivariate.continuous.Pillai.F     = F_continuous_Pillai;
        model.multivariate.continuous.Pillai.p     = pval_continuous_Pillai;
        model.multivariate.continuous.Pillai.R2.R2 = R2_Pillai;
        model.multivariate.continuous.Pillai.R2.F  = R2_F_Pillai;
        model.multivariate.continuous.Pillai.R2.p  = R2_pval_Pillai;
        model.multivariate.continuous.Roy.theta    = theta_continuous;
        model.multivariate.continuous.Roy.F        = F_continuous_Roy;
        model.multivariate.continuous.Roy.p        = pval_continuous_Roy;
        model.multivariate.continuous.Roy.R2.R2    = R2_Roy;
        model.multivariate.continuous.Roy.R2.F     = F_continuous_Roy;
        model.multivariate.continuous.Roy.R2.p     = pval_continuous_Roy;
        model.multivariate.continuous.Hotelling.U  = U_continuous;
        model.multivariate.continuous.Hotelling.F  = F_continuous_Hotelling;
        model.multivariate.continuous.Hotelling.p  = pval_continuous_Hotelling;

    elseif nb_conditions > 0 && nb_continuous > 0
        
        model.multivariate.EV                     = [Eigen_values_R2'; Eigen_values_cond'; Eigen_values_continuous'];
        model.multivariate.conditions.Pillai.V    = V;
        model.multivariate.conditions.Pillai.F    = F_conditions_Pillai;
        model.multivariate.conditions.Pillai.p    = pval_conditions_Pillai;
        model.multivariate.conditions.Pillai.df   = df_Pillai;
        model.multivariate.conditions.Pillai.dfe  = dfe_Pillai;
        model.multivariate.conditions.Roy.theta   = theta;
        model.multivariate.conditions.Roy.F       = F_conditions_Roy;
        model.multivariate.conditions.Roy.p       = pval_conditions_Roy;
        model.multivariate.conditions.Roy.df      = df_Roy;
        model.multivariate.conditions.Roy.dfe     = dfe_Roy;
        model.multivariate.conditions.Hotelling.U  = U;
        model.multivariate.conditions.Hotelling.F  = F_conditions_Hotelling;
        model.multivariate.conditions.Hotelling.p  = pval_conditions_Hotelling;
        model.multivariate.conditions.Hotelling.df = df_Hotelling;
        model.multivariate.conditions.Hotelling.dfe= dfe_Hotelling;

        model.multivariate.continuous.df          = df_continuous;
        model.multivariate.continuous.dfe         = dfe_continuous;
        model.multivariate.continuous.Pillai.V    = V_continuous;
        model.multivariate.continuous.Pillai.F    = F_continuous_Pillai;
        model.multivariate.continuous.Pillai.p    = pval_continuous_Pillai;
        model.multivariate.continuous.Roy.theta   = theta_continuous;
        model.multivariate.continuous.Roy.F       = F_continuous_Roy;
        model.multivariate.continuous.Roy.p       = pval_continuous_Roy;
        model.multivariate.continuous.Hotelling.U = U_continuous;
        model.multivariate.continuous.Hotelling.F = F_continuous_Hotelling;
        model.multivariate.continuous.Hotelling.p = pval_continuous_Hotelling;
    end
end 
 

function model = limo_mglm(varargin)

% Multivariate General Linear Model for EEG data
% The model consider subjects / trials as independent observations
% i.e. this is similar as running a N-way MANOVA or MANCOVA
% Analyses are performed at one time frame only, but for all
% subjects/ trials and considering the pattern of effect across all electrodes.
% In addition, we compute the multivariate Discriminable value D which allows
% finding if conditions/regressors (classes) are different. This is similar
% as a pattern classification (LDA, SVM etc) without to have to reply on
% particular classifier algorithm.
%
% FORMAT:
% model = limo_mglm(Y,LIMO)
% model = limo_mglm(Y,X,nb_conditions,nb_interactions,nb_continuous,method, cov_method,W)
%
% INPUTS:
%   Y             = 2D matrix of EEG data with format trials/subjects x electrodes
%   LIMO          = structure that contains the information below (except Y)
%   or input info that are in LIMO.mat
%   X             = 2 dimensional design matrix
%   nb_conditions = a vector indicating the number of conditions per factor
%   nb_interactions = a vector indicating number of columns per interactions
%   nb_continuous = number of covariates
%   method        = 'OLS', 'WLS', 'IRLS' (bisquare)
%   cov_mehtod    = either pseudo inverse of cov or inverse on regularized cov
%   W             = optional - a matrix a trial weights ; ie method 'WLS'
%                   using these weights rather than something else
%
% OUTPUTS:
%     model.R2.V
%     model.R2.EV = eigenvalues
%     model.R2.Roy.F
%     model.R2.Roy.p
%     model.R2.Pillai.F
%     model.R2.Pillai.p
%
%     model.betas  = the beta parameters (dimensions nb of paramters x electrodes)
%     model.conditions    = categorical effects
%          --> F/p in rows are the factors
%          --> df row 1 = df, row 2 = dfe, columns are factors
%     model.continuous = continuous effects
%          --> F/p in rows are the variables
%          --> df column 1 = df, column2 2 = dfe (same for all covariates)
%
%     model.Discriminant.Z = discriminant functions (for plotting observations)
%          --> [trials/subjects x nb discriminant functions]
%     model.Discriminant.Z_cent = centered discriminant functions (for plotting observations)
%     model.Discriminant.importanceZ = relative importance of each discrim function
%          --> considers its eigenvalue as a proportion of the total
%     model.Discriminant.cc_squared = squared canonical correlation
%          --> corresponding to each discriminant function
%     model.Discriminant.D = Discriminable values (standardized eigenvector coefficients) for plotting topography
%          --> Identifying the relative contribution of the electrodes to seperation of the classes.
%          In other words, find spatial patterns that distinguish classes
%          [nb of electrodes x nb of discriminant functions]
%
%     model.Classification.Linear.Acc = decoding accuracy
%     model.Classification.Linear.cvAcc = 10-fold cross validated decoding accuracy
%     model.Classification.Linear.CvSD = sd error of cross validated accuracy
%     model.Classification.Quadratic.Acc = decoding accuracy
%     model.Classification.Quadratic.cvAcc = 10-fold cross validated decoding accuracy
%     model.Classification.Quadratic.CvSD = sd error of cross validated accuracy
%
% NOTES:
%
% The MANOVA uses Pillais' V and Roy's tests.
% In Roy's test one maximizes the spread of the transformed data using the
% 1st eigen value of inv(E)*H. The F approximation is an upper bound, i.e.
% results are safe if H0 is accepted (no effect) but not rejected. If data
% are highly correlated, there is only one high eigen value and Roy test is
% appropriate otherwise Pillai is better.
%
% The parameters can be computed using 3 methods: ordinary least squares,
% weighted least squares and iterative reweighted least squares.
% For OLS, W=ones(size(Y)), for IRLS and WLS the weights can be provided.
% For IRLS, each cell (time and space) is weighted so it doesn't matter if
% precomputed by limo_eeg.m *unless you want to pass some special weights)
% For WLS, the default from limo_eeg.m is that weights are computed
% trial-wise, seeing a trial as multivariate in time, but then the MANOVA is computed
% across electrodes (ie all electrode have the same weights). If W is not provided,
% the weights will be computed across electrodes, seeing a trial as a multivariate
% measures in space.
%
% Within the MGLM, each effect is accounted for given the other effects; this means
% that one can have a different number of trials/subjects per conditions/factors
% provided there is no interactions - for interactions this is not possible
% and no correction is provided - a design created by limo_design_matrix
% would have sampled trials to make sure the number of trials is identical
% across interaction terms.
%
% References
% - Christensen, R. 2002. Plane answers to complex questions. 3rd Ed. Springer-Verlag
% - Friston et al. 2007. Statitical Parametric Mapping. Academic Press
% - Rencher, A.C. 2002. Methods of Multivariate Analysis. 2nd Ed. Wiley.
% - Robert and Escoufier, 1976. A Unifying Tool for Linear Multivariate Statistical Methods:
%   The RV- Coefficient J.Royal Stat Soc, C
% - Bassez et al. Spatial Multivariate Analyses of EEG data: theory and application. In prep.
%
% See also
% LIMO_DESIGN_MATRIX, LIMO_PCOUT, LIMO_IRLS, LIMO_EEG, LIMO_DECOMP
%
% Cyril Pernet v1 May 2013
% Iege Bassez  v2 June 2018
% ----------------------------
% Copyright (C) LIMO Team 2018

%% defaults
decomp_method = 'pseudo'; % default decomposition method ie eig((pinv(E)*H))
folds = 10; % 10 folds cross validated classification

%% varagin

W = [];
Y                  = varargin{1}; % the data
if nargin == 2
    X               = varargin{2}.design.X;
    nb_conditions   = varargin{2}.design.nb_conditions;
    nb_interactions = varargin{2}.design.nb_interactions;
    nb_continuous   = varargin{2}.design.nb_continuous;
    method          = varargin{2}.design.method;
    cov_method      = varargin{2}.design.cov_method;
    if isfield(varargin{2}.design,'weigths')
        W           = varargin{2}.design.weigths;
    end
elseif nargin >= 7
    X               = varargin{2};
    nb_conditions   = varargin{3};
    nb_interactions = varargin{4};
    nb_continuous   = varargin{5};
    method          = varargin{6};
    cov_method      = varargin{7};
    if nargin == 8
        W           = varargin(8);
    end
else
    error('varargin error')
end

nb_factors = numel(nb_conditions);
if nb_factors == 1 && nb_conditions == 0
    nb_factors = 0;
end

% -----------
%% Data check
% -----------

if size(Y,1)~=size(X,1)
    error('The number of events in Y and the design matrix are different')
elseif size(Y,1)>=size(X,1)
    error('The number of events is bigger or equal to the number of channels, the analysis cannot be done')
end

if nb_interactions == 0
    nb_interactions = [];
end

% ------------------------------
%% Compute F for dummy variables
% ------------------------------

% -------------------------
if nb_factors == 1   %  1-way MANOVA/MANCOVA
    % -------------------------
    
    % total sum of squares, projection matrix for errors, residuals and betas
    % -----------------------------------------------------------------------
    T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
    R     = eye(size(Y,1)) - (X*pinv(X));                                      % Projection on E
    E     = (Y'*R*Y);                                                          % SS Error
    
    % compute Beta parameters and weights
    if strcmp(method,'OLS')
        W = ones(size(Y,1),1);
        Betas = pinv(X)*Y;
    elseif strcmp(method,'WLS')
        if isempty(W)
            [Betas,W] = limo_WLS(X,Y);
        else
            Betas = pinv(WX)*WY;
        end
    elseif strcmp(method,'IRLS')
        [Betas,W] = limo_IRLS(X,Y);
    end
    
    % compute model R^2
    % -----------------
    C = eye(size(X,2));
    C(:,size(X,2)) = 0;                  % contrast matrix for effect
    C0 = eye(size(X,2)) - C*pinv(C);     % contrast matrix for error
    X0 = X*C0;                           % reduced model
    R0 = eye(size(Y,1)) - (X0*pinv(X0));
    M  = R0 - R;                         % projection matrix onto Xc
    H = (Betas'*X'*M*X*Betas);           % SS Effect
                                         % only works with rank deficient matrix
                                         % intercept column with ones as last column of X
    if round(H,6) ~= round(T - E, 6)     % since T = H + E
        H = T - E;                       % 'temporary' solution ; is there a better solution?
    end
    
    
    % Generalized R2 - Robert and Escoufier, J.Royal Stat Soc, C - 1976
    % -----------------------------------------------------------------
    % variance covariance matrix
    S   = cov([Y X(:,1:size(X,2)-1)]);
    Syy = S(1:size(Y,2),1:size(Y,2));
    Sxy = S(size(Y,2)+1:size(S,1),1:size(Y,2));
    Syx = S(1:size(Y,2),size(Y,2)+1:size(S,2));
    Sxx = S(size(Y,2)+1:size(S,1),size(Y,2)+1:size(S,2));
    Rsquare_multi = trace(Sxy*Syx) / sqrt(trace(Sxx.^2)*trace(Syy.^2));
    
    % compute multivariate effects
    % ----------------------------
    if strcmp(cov_method, 'regularized')
        [RegularizedCovariance, ~] = cov1para(Y); % shrinkage
        E = RegularizedCovariance .* (size(Y,1)-nb_conditions);
    end
    
    [~,pd] = chol(E); % check if E is positive definite
    if pd == 0
        decomp_method = 'choldesky';
    end
    [Eigen_vectors_R2, Eigen_values_R2] = limo_decomp(E,H, decomp_method);
    
    p = size(Y,2);      % = number of variables (dimension)
    q = rank(X)-1;      % -1 because of intercept column
    s = min(p,q);       % df
    n = size(Y,1);      % nb of observations (dfe)
    m = (abs(q-p)-1)/2;
    N = (n-q-p-2)/2;
    ve = n - rank(X);
    
    % Roy
    theta        = max(Eigen_values_R2) / (1+max(Eigen_values_R2));
    R2_Roy_value = theta; % = 1st canonical correlation
    R2_Roy_F     = ((ve-max(p,q)+q)*max(Eigen_values_R2))/max(p,q);
    R2_Roy_p     = 1-fcdf(R2_Roy_F, max(p,q), ve-max(p,q)+q);
    
    % Pillai
    V               = sum(Eigen_values_R2 ./ (1+Eigen_values_R2));
    R2_Pillai_value = V / s; % average of canonical correlations
    R2_Pillai_F     = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
    R2_Pillai_p     = 1-fcdf(R2_Pillai_F,(s*(2*m+s+1)),(s*(2*N+s+1)));
    
    % compute F for categorical variables
    % -----------------------------------
    if nb_conditions ~= 0 && nb_continuous == 0
        Eigen_values_cond  = Eigen_values_R2;
        Eigen_vectors_cond = Eigen_vectors_R2;
        
    elseif nb_conditions ~= 0 && nb_continuous ~= 0
        C = eye(size(X,2));
        C(:,(nb_conditions+1):size(X,2)) = 0;
        C0 = eye(size(X,2)) - C*pinv(C);
        X0 = X*C0; % Here the reduced model includes the covariates
        R0 = eye(size(Y,1)) - (X0*pinv(X0));
        M  = R0 - R;
        H  = (Betas'*X'*M*X*Betas);
        [Eigen_vectors_cond, Eigen_values_cond] = limo_decomp(E,H,decomp_method);
    end
    
    vh = nb_conditions - 1; % df = q above
    s = min(vh,p); % subspace in which mean Ys are located
    for c=1:nb_conditions
        nb_items(c) = numel(find(X(:,c)));
    end
    
    if sum(nb_items == nb_items(1)) == length(nb_items)
        ve = nb_conditions*(nb_items(1)-1);     % dfe equal sample sizes
    else
        ve = sum(nb_items) - nb_conditions;     % dfe different sample sizes
    end
    
    if s > 1 % (see page 165 and p.166 Rencher)
        m = (abs(vh-p)-1)/2;
        N = (ve-p-1) / 2;
        
        % Pillai
        V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
        df_conditions_Pillai   = s*(2*m+s+1);
        dfe_conditions_Pillai  = s*(2*N+s+1);
        F_conditions_Pillai    = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
        pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_conditions_Pillai,dfe_conditions_Pillai);
        
        % Roy's test
        theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond));
        df_conditions_Roy   = max(p,vh);
        dfe_conditions_Roy  = ve - max(p,q) + q;
        F_conditions_Roy    = (dfe_conditions_Roy*max(Eigen_values_cond))/df_conditions_Roy;
        pval_conditions_Roy = 1-fcdf(F_conditions_Roy, df_conditions_Roy, dfe_conditions_Roy);
        
    else % = only one non zeros Eigen value s = 1 and/or vh = 1 (see p. 169 Rencher)
        
        theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond)); % Roy
        V = theta; % Pillai(1) equals theta
        
        df_conditions_Pillai  = p;
        dfe_conditions_Pillai = ve-p+1;
        df_conditions_Roy     = df_conditions_Pillai;
        dfe_conditions_Roy    = dfe_conditions_Pillai;
        
        F_conditions_Pillai    = (dfe_conditions_Pillai/df_conditions_Pillai) * max(Eigen_values_cond);
        pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_conditions_Pillai,dfe_conditions_Pillai);
        F_conditions_Roy       = F_conditions_Pillai;
        pval_conditions_Roy    = pval_conditions_Pillai;
    end
    
    %% Discriminant Analysis
    % ------------------------
    if length(Y)-nb_conditions <= nb_conditions
        errordlg('there is not enough data point to run a discriminant analysis')
    else
        
        % rescale eigenvectors so the within-group variance is 1
        n = size(Y,1); % nb of observations (dfe)
        q = rank(X); % number of groups (df)
        a = Eigen_vectors_cond;
        vs = diag((a' * E * a))' ./ (n-q);
        vs(vs<=0) = 1;
        a = a ./ repmat(sqrt(vs), size(a,1), 1);
        scaled_eigenvectors = a;
        
        % validate if correct eigenvectors
        if round((pinv(E)*H) * scaled_eigenvectors(:,1), 4) ~= round(Eigen_values_cond(1) * scaled_eigenvectors(:,1), 4);
            errordlg('something went wrong with scaling the eigenvectors')
        end
        
        % get the discriminant function(s)
        scaled_eigenvectors = scaled_eigenvectors(:,1:s); % s: nb nonzero eigenvalues
        Eigen_values_cond   = Eigen_values_cond(1:s); % s: nb nonzero eigenvalues
        z                   = Y * scaled_eigenvectors; % the discriminant functions themself
        z_importance        = Eigen_values_cond ./ sum(Eigen_values_cond); % variance corresponding to each eigenvalue
        cc_squared          = Eigen_values_cond ./ (1 + Eigen_values_cond); % canonical correlation
        centered_Y          = bsxfun(@minus, Y,mean(Y));
        centered_z          = centered_Y * scaled_eigenvectors; % centered discriminant functions
        
        % now, get standardized eigenvectors
        Spooled = E./ve;
        if strcmp(cov_method, 'regularized')
            Spooled = RegularizedCovariance;
        end
        standardized_eigenvectors = bsxfun(@times, scaled_eigenvectors, sqrt(diag(Spooled)));
    end
    
    %% do the classification
    %--------------------------------------------------------
    % get training linear decoding accuracy:
    [class,~]            = find(X(:,1:nb_conditions)');
    [~, ~,training_Acc]  = limo_LDA(Y, class, Y, cov_method);
    % confmat      = confusionmat(class, predicted); clear predicted;
    % training_Acc = sum(diag(confmat/sum(sum(confmat))));
    
    % get 10-fold CV linear decoding accuracy:
    cvp = cvpartition(class,'k',folds);
    accuracyvector = NaN(1,folds);
    for foldi=1:folds
        trainIdx              = cvp.training(foldi);
        testIdx               = cvp.test(foldi);
        % make a training and test set based on indexes:
        trainingset           = Y(trainIdx,:);
        testset               = Y(testIdx,:);
        % class labels:
        traininglabel         = class(trainIdx,1);
        testlabel             = class(testIdx,1);
        % Classification:
        predicted             = limo_LDA(trainingset,traininglabel,testset, cov_method);
        confmat               = confusionmat(testlabel, predicted); clear predicted;
        accuracyvector(foldi) = sum(diag(confmat/sum(sum(confmat))));
    end
    cvAcc  = mean(accuracyvector);
    cvSD   = std(accuracyvector);
    
    %     % get training linear decoding accuracy:
    %     LinearModel = fitcdiscr(Y,class, 'DiscrimType', 'pseudolinear'); %The software inverts the covariance matrix using the pseudo inverse.
    %     training_Acc = sum(diag(confusionmat(LinearModel.Y, predict(LinearModel, Y))))/sum(sum(confusionmat(LinearModel.Y, predict(LinearModel, Y))));
    %
    %     % get CV linear decoding accuracy:
    %     LinearModel_CV = fitcdiscr(Y,class, 'DiscrimType', 'pseudolinear','CrossVal', 'on');
    %     acc = NaN(1,LinearModel_CV.KFold);
    %     for k=1:LinearModel_CV.KFold
    %         acc(k) = sum(diag(confusionmat(LinearModel_CV.Y, predict(LinearModel_CV.Trained{k,1}, Y))))/sum(sum(confusionmat(LinearModel_CV.Y, predict(LinearModel_CV.Trained{k,1}, Y))));
    %     end
    %     cvAcc = mean(acc);
    %     cvSD  = sqrt(sum((acc - mean(acc)).^2)/(LinearModel_CV.KFold-1));
    %
    
    % get training quadratic decoding accuracy:
    QuadraticModel = fitcdiscr(Y, class, 'DiscrimType', 'pseudoquadratic');
    q_training_Acc = sum(diag(confusionmat(QuadraticModel.Y, predict(QuadraticModel, Y))))/sum(sum(confusionmat(QuadraticModel.Y, predict(QuadraticModel, Y))));
    
    % get CV quadratic decoding accuracy:
    QuadraticModel_CV = fitcdiscr(Y, class, 'DiscrimType', 'pseudoquadratic', 'CrossVal', 'on');
    q_acc = NaN(1,QuadraticModel_CV.KFold);
    for k=1:QuadraticModel_CV.KFold
        q_acc(k) = sum(diag(confusionmat(QuadraticModel_CV.Y, predict(QuadraticModel_CV.Trained{k,1}, Y))))/sum(sum(confusionmat(QuadraticModel_CV.Y, predict(QuadraticModel_CV.Trained{k,1}, Y))));
    end
    q_cvAcc = mean(q_acc);
    q_cvSD  = sqrt(sum((q_acc - mean(q_acc)).^2)/(QuadraticModel_CV.KFold-1));
    
    %%
    % ------------------------------------------------
elseif nb_factors > 1  && isempty(nb_interactions) % N-ways MANOVA without interactions
    % ------------------------------------------------
    
    % compute basic SS total, projection matrices and parameters
    T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
    R        = eye(size(Y,1)) - (X*pinv(X));
    E        = (Y'*R*Y);

    % compute Beta parameters and weights
    if strcmp(method,'OLS')
        W = ones(size(Y,1),1);
        Betas = pinv(X)*Y;
    elseif strcmp(method,'WLS')
        if isempty(W)
            [Betas,W] = limo_WLS(X,Y);
        else
            Betas = pinv(WX)*WY;
        end
    elseif strcmp(method,'IRLS')
        [Betas,W] = limo_IRLS(X,Y);
    end
    
    % --------------------
    % compute model R^2
    % --------------------
    C = eye(size(X,2));
    C(:,size(X,2)) = 0;
    C0   = eye(size(X,2)) - C*pinv(C);
    X0   = X*C0; % Reduced model (i.e. only intercept)
    R0   = eye(size(Y,1)) - (X0*pinv(X0));
    M    = R0 - R;      % M is the projection matrix onto Xc
    H    = (Betas'*X'*M*X*Betas);   % SSCP Hypothesis (Effect)
    [Eigen_vectors_R2, Eigen_values_R2] = limo_decomp(E,H);
    
    % Generalized R2
    % variance covariance matrix
    S = cov([Y X(:,1:size(X,2)-1)]);
    Syy = S(1:size(Y,2),1:size(Y,2));
    Sxy = S(size(Y,2)+1:size(S,1),1:size(Y,2));
    Syx = S(1:size(Y,2),size(Y,2)+1:size(S,2));
    Sxx = S(size(Y,2)+1:size(S,1),size(Y,2)+1:size(S,2));
    Rsquare_multi = trace(Sxy*Syx) / sqrt(trace(Sxx.^2)*trace(Syy.^2)); % Robert and Escoufier, J.Royal Stat Soc, C - 1976
    
    p = size(Y,2); % = number of variables (dimension)
    q = rank(X); % = number of regressors (df)
    s = min(p,q); % df
    n = size(Y,1); % nb of observations (dfe)
    m = (abs(q-p)-1)/2;
    N = (n-q-p-2)/2;
    d = max(p,q);
    
    theta = max(Eigen_values_R2) / (1+max(Eigen_values_R2)); % Roy
    R2_Roy_value = theta; % = 1st canonical correlation
    R2_Roy_F     = ((n-d-1)*max(Eigen_values_R2))/d;
    R2_Roy_p     = 1-fcdf(R2_Roy_F, d, (n-d-1));
    
    V = sum(Eigen_values_R2 ./ (1+Eigen_values_R2)); % Pillai
    R2_Pillai_value = V / s; % average of canonical correlations
    R2_Pillai_F     = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
    R2_Pillai_p     = 1-fcdf(R2_Pillai_F,(s*(2*m+s+1)),(s*(2*N+s+1)));
    
    % --------------------------------------
    % compute F and p values of each factor
    % --------------------------------------
    
    eoi = zeros(1,size(X,2));
    eoi(1:nb_conditions(1)) = 1:nb_conditions(1);
    eoni = [1:size(X,2)];
    eoni = find(eoni - eoi);
    model.conditions.EV = [];
    for f = 1:length(nb_conditions)
        C = eye(size(X,2));
        C(:,eoni) = 0;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0;
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;
        H    = (Betas'*X'*M*X*Betas);
        [Eigen_vectors_cond, Eigen_values_cond] = limo_decomp(E,H);
        model.conditions.EV(f,:) = Eigen_values_cond';
        
        vh = nb_conditions(f) - 1; % df = q above
        s = min(vh,p); % subspace in which mean Ys are located
        clear nb_items; x = X(:,find(eoi));
        for c=1:nb_conditions(f)
            nb_items(c) = numel(find(x(:,c)));
        end
        
        if sum(nb_items == nb_items(1)) == length(nb_items)
            ve = nb_conditions(f)*(nb_items(1)-1);     % dfe equal sample sizes
        else
            ve = sum(nb_items) - nb_conditions(f);     % dfe different sample sizes
        end
        
        if s > 1
            m = (abs(vh-p)-1)/2;
            N = (ve-p-1) / 2;
            
            % Pillai
            V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
            df_conditions_Pillai(f) = s*(2*m+s+1);
            dfe_conditions_Pillai(f) = s*(2*N+s+1);
            F_conditions_Pillai(f) = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
            pval_conditions_Pillai(f) = 1-fcdf(F_conditions_Pillai(f),df_conditions_Pillai(f),dfe_conditions_Pillai(f));
            
            % Roy's test
            theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond));
            df_conditions_Roy(f) = max(p,vh);
            dfe_conditions_Roy(f) = ve - 1; % in Renchner it is proposed to use ve - max(p,vh) -1 while in Statistica it is ve -1
            F_conditions_Roy(f) = (dfe_conditions_Roy(f)*max(Eigen_values_cond))/df_conditions_Roy(f);
            pval_conditions_Roy(f) = 1-fcdf(F_conditions_Roy(f), df_conditions_Roy(f), dfe_conditions_Roy(f));
            
        else % = only one non zeros Eigen value s = 1 and/or vh = 1
            
            V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
            U = max(Eigen_values_cond);
            theta = U;
            
            df_conditions_Pillai(f) = p; % number of electrodes
            dfe_conditions_Pillai(f) = ve-p+1;
            df_conditions_Roy(f) = df_conditions_Pillai(f);
            
            F_conditions_Pillai(f) = (dfe_conditions_Pillai(f)/df_conditions_Pillai(f)) * max(Eigen_values_cond);
            pval_conditions_Pillai(f) = 1-fcdf(F_conditions_Pillai(f),df_conditions_Pillai(f),dfe_conditions_Pillai(f));
            F_conditions_Roy(f) = F_conditions_Pillai(f);
            pval_conditions_Roy(f) = pval_conditions_Pillai(f);
        end
                
        % update factors
        if f<length(nb_conditions)
            update = max(find(eoi));
            eoi = zeros(1,size(X,2));
            eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
            eoni = [1:size(X,2)];
            eoni = find(eoni - eoi);
        end
    end
    
    % ------------------------------------------------
elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways MANOVA with interactions
    % ------------------------------------------------
    
    % compute basic SS total, projection matrices and parameters
    T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
    R        = eye(size(Y,1)) - (X*pinv(X));
    E        = (Y'*R*Y);
    if strcmp(method,'OLS')
        W = ones(size(Y,1),1);
        Betas = pinv(X)*Y;
    elseif strcmp(method,'WLS')
        if isempty(W)
            [Betas,W] = limo_WLS(X,Y);
        else
            Betas = pinv(WX)*WY;
        end
    elseif strcmp(method,'IRLS')
        [Betas,W] = limo_IRLS(X,Y);
    end
    
    % --------------------
    % compute model R^2
    % --------------------
    C = eye(size(X,2));
    C(:,size(X,2)) = 0;
    C0   = eye(size(X,2)) - C*pinv(C);
    X0   = X*C0; % Reduced model (i.e. only intercept)
    R0   = eye(size(Y,1)) - (X0*pinv(X0));
    M    = R0 - R;      % M is the projection matrix onto Xc
    H    = (Betas'*X'*M*X*Betas);   % SSCP Hypothesis (Effect)
    [Eigen_vectors_cond, Eigen_values_cond] = limo_decomp(E,H);
    
    % Generalized R2
    % variance covariance matrix
    S = cov([Y X(:,1:size(X,2)-1)]);
    Syy = S(1:size(Y,2),1:size(Y,2));
    Sxy = S(size(Y,2)+1:size(S,1),1:size(Y,2));
    Syx = S(1:size(Y,2),size(Y,2)+1:size(S,2));
    Sxx = S(size(Y,2)+1:size(S,1),size(Y,2)+1:size(S,2));
    Rsquare_multi = trace(Sxy*Syx) / sqrt(trace(Sxx.^2)*trace(Syy.^2)); % Robert and Escoufier, J.Royal Stat Soc, C - 1976
    
    [Eigen_vectors_R2, Eigen_values_R2] = limo_decomp(E,H);
    p = size(Y,2); % = number of variables (dimension)
    q = rank(X); % = number of regressors (df)
    s = min(p,q); % df
    n = size(Y,1); % nb of observations (dfe)
    m = (abs(q-p)-1)/2;
    N = (n-q-p-2)/2;
    d = max(p,q);
    
    theta = max(Eigen_values_R2) / (1+max(Eigen_values_R2)); % Roy
    R2_Roy_value = theta; % = 1st canonical correlation
    R2_Roy_F     = ((n-d-1)*max(Eigen_values_R2))/d;
    R2_Roy_p     = 1-fcdf(R2_Roy_F, d, (n-d-1));
    
    V = sum(Eigen_values_R2 ./ (1+Eigen_values_R2)); % Pillai
    R2_Pillai_value = V / s; % average of canonical correlations
    R2_Pillai_F     = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
    R2_Pillai_p     = 1-fcdf(R2_Pillai_F,(s*(2*m+s+1)),(s*(2*N+s+1)));
    
    % ---------------------------------------------------
    % start by ANOVA without interaction for main effects
    % ---------------------------------------------------
    
    % covariates
    covariate_columns = [(sum(nb_conditions)+sum(nb_interactions)+1):(size(X,2)-1)];
    
    % main effects
    dummy_columns = 1:sum(nb_conditions);
    
    % re-define X
    x = [X(:,dummy_columns) X(:,covariate_columns) ones(size(X,1),1)];
    
    % run same model as above
    R        = eye(size(Y,1)) - (x*pinv(x));
    if strcmp(method,'IRLS')
        betas = pinv(Wx)*WY;
    else
        w = repamt(W,1,size(Y,2));
        betas = pinv(wx)*wY;
    end
    
    eoi = zeros(1,size(x,2));
    eoi(1:nb_conditions(1)) = 1:nb_conditions(1);
    eoni = [1:size(x,2)];
    eoni = find(eoni - eoi);
    model.conditions.EV = [];
    
    for f = 1:length(nb_conditions)
        C = eye(size(x,2));
        C(:,eoni) = 0;
        C0   = eye(size(x,2)) - C*pinv(C);
        X0   = x*C0;
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;
        H(f,:) = diag((betas'*x'*M*x*betas));
        [Eigen_vectors_cond, Eigen_values_cond] = limo_decomp(E,H);
        model.conditions.EV(f,:) = Eigen_values_cond';
        
        vh = nb_conditions(f) - 1; % df = q above
        s = min(vh,p); % subspace in which mean Ys are located
        clear nb_items; x = X(:,find(eoi));
        for c=1:nb_conditions(f)
            nb_items(c) = numel(find(x(:,c)));
        end
        
        if sum(nb_items == nb_items(1)) == length(nb_items)
            ve = nb_conditions(f)*(nb_items(1)-1);     % dfe equal sample sizes
        else
            ve = sum(nb_items) - nb_conditions(f);     % dfe different sample sizes
        end
        
        if s > 1
            m = (abs(vh-p)-1)/2;
            N = (ve-p-1) / 2;
            
            % Pillai
            V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
            df_conditions_Pillai(f) = s*(2*m+s+1);
            dfe_conditions_Pillai(f) = s*(2*N+s+1);
            F_conditions_Pillai(f) = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
            pval_conditions_Pillai(f) = 1-fcdf(F_conditions_Pillai(f),df_conditions_Pillai(f),dfe_conditions_Pillai(f));
            
            % Roy's test
            theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond));
            df_conditions_Roy(f) = max(p,vh);
            dfe_conditions_Roy(f) = ve - 1; % in Renchner it is proposed to use ve - max(p,vh) -1 while in Statistica it is ve -1
            F_conditions_Roy(f) = (dfe_conditions_Roy(f)*max(Eigen_values_cond))/df_conditions_Roy(f);
            pval_conditions_Roy(f) = 1-fcdf(F_conditions_Roy(f), df_conditions_Roy(f), dfe_conditions_Roy(f));
            
        else % = only one non zeros Eigen value s = 1 and/or vh = 1
            
            V = sum(Eigen_values_cond ./ (1+Eigen_values_cond));
            U = max(Eigen_values_cond);
            theta = U;
            
            df_conditions_Pillai(f) = p; % number of electrodes
            dfe_conditions_Pillai(f) = ve-p+1;
            df_conditions_Roy(f) = df_conditions_Pillai;
            dfe_conditions_Roy(f) = dfe_conditions_Pillai;
            
            F_conditions_Pillai(f) = (dfe_conditions_Pillai/df_conditions_Pillai) * max(Eigen_values_cond);
            pval_conditions_Pillai(f) = 1-fcdf(F_conditions_Pillai(f),df_conditions_Pillai(f),dfe_conditions_Pillai(f));
            F_conditions_Roy(f) = F_conditions_Pillai(f);
            pval_conditions_Roy(f) = pval_conditions_Pillai(f);
        end
        
        % update factors
        if f<length(nb_conditions)
            update = max(find(eoi));
            eoi = zeros(1,size(x,2));
            eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
            eoni = [1:size(x,2)];
            eoni = find(eoni - eoi);
        end
    end
    
    % ---------------------------
    % now deal with interactions
    % ---------------------------
    
    if nb_factors == 2 && nb_continuous == 0 % the quick way with only one interaction
        HI = diag(T)' - H(1,:) - H(2,:) - diag(E)';
        [Eigen_vectors_inter, Eigen_values_inter] = limo_decomp(E,HI);
        model.interactions.EV = [Eigen_values_inter'];
        
        vh = nb_interactions - 1; % df = q above
        s = min(vh,p); % subspace in which mean Ys are located
        clear nb_items; x = X(:,sum(nb_conditions)+1:end-1);
        for c=1:nb_interactions
            nb_items(c) = numel(find(x(:,c)));
        end
        
        if sum(nb_items == nb_items(1)) == length(nb_items)
            ve = nb_interactions*(nb_items(1)-1);     % dfe equal sample sizes
        else
            ve = sum(nb_items) - nb_interactions;     % dfe different sample sizes
        end
        
        if s > 1
            m = (abs(vh-p)-1)/2;
            N = (ve-p-1) / 2;
            
            % Pillai
            V = sum(Eigen_values_inter ./ (1+Eigen_values_inter));
            df_interaction_Pillai = s*(2*m+s+1);
            dfe_interaction_Pillai = s*(2*N+s+1);
            F_interaction_Pillai = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
            pval_interaction_Pillai = 1-fcdf(F_interaction_Pillai,df_interaction_Pillai,dfe_interaction_Pillai);
            
            % Roy's test
            theta = max(Eigen_values_inter) / (1+max(Eigen_values_inter));
            df_interaction_Roy = max(p,vh);
            dfe_interaction_Roy = ve - 1; % in Renchner it is proposed to use ve - max(p,vh) -1 while in Statistica it is ve -1
            F_interaction_Roy = (dfe_interaction_Roy*max(Eigen_values_cond))/df_interaction_Roy;
            pval_interaction_Roy = 1-fcdf(F_interaction_Roy, df_interaction_Roy, dfe_interaction_Roy);
            
        else % = only one non zeros Eigen value s = 1 and/or vh = 1
            
            V = sum(Eigen_values_inter ./ (1+Eigen_values_inter));
            U = max(Eigen_values_inter);
            theta = U;
            
            df_interaction_Pillai = p; % number of electrodes
            dfe_interaction_Pillai = ve-p+1;
            df_interaction_Roy = df_interaction_Pillai;
            dfe_interaction_Roy = dfe_interaction_Pillai;
            
            F_interaction_Pillai = (dfe_interaction_Pillai/df_interaction_Pillai) * max(Eigen_values_cond);
            pval_interaction_Pillai = 1-fcdf(F_interaction_Pillai,df_interaction_Pillai,dfe_interaction_Pillai);
            F_interaction_Roy = F_interaction_Pillai;
            pval_interaction_Roy = pval_interaction_Pillai;
        end
        
    else % run through each interaction
        
        % part of X unchanged
        Main_effects = [X(:,dummy_columns)];
        Cov_and_Mean = [X(:,covariate_columns) ones(size(Y,1),1)];
        
        % get interactions
        start = size(Main_effects,2)+1;
        for i=1:length(nb_interactions)
            I{i} = X(:,start:(start+nb_interactions(i)-1));
            start = start+nb_interactions(i);
        end
        start = size(Main_effects,2)+1;
        
        % check interaction levels
        index = 1;
        for n=2:nb_factors
            combinations = nchoosek([1:nb_factors],n); % note it matches I above because computed with nchoosek the same way in limo_design_matrix
            for c = 1:size(combinations,1)
                interaction{index} = combinations(c,:);
                index = index + 1;
            end
        end
        
        add = 0; start_at_I = 1;
        model.interactions.EV = [];
        % run substituting and/or incrementing parts of X
        for f = 1:length(nb_interactions)
            
            % re-define X with interactions
            test = size(interaction{f},2);
            if test == 2
                x = [Main_effects I{f} Cov_and_Mean]; % havbing the same nb of trials simply additive model
                add = add+1;
            else
                if add == test
                    for a = start_at_I:add
                        Main_effects = [Main_effects I{a}];
                    end
                    start = size(Main_effects,2)+1;
                    start_at_I = add+1;
                end
                x = [Main_effects I{f} Cov_and_Mean];
            end
            
            % run same model as above
            R  = eye(size(Y,1)) - (x*pinv(x));
            if strcmp(method,'IRLS')
                betas = pinv(Wx)*WY;
            else
                w = repamt(W,1,size(Y,2));
                betas = pinv(wx)*wY;
            end
            
            eoi = zeros(1,size(x,2));
            eoi(start:(start-1+nb_interactions(f))) = start:(start-1+nb_interactions(f));
            eoni = [1:size(x,2)];
            eoni = find(eoni - eoi);
            
            C = eye(size(x,2));
            C(:,eoni) = 0;
            C0   = eye(size(x,2)) - C*pinv(C);
            X0   = x*C0;
            R0   = eye(size(Y,1)) - (X0*pinv(X0));
            M    = R0 - R;
            HI(f,:) = diag((betas'*x'*M*x*betas))';
            [Eigen_vectors_inter, Eigen_values_inter] = limo_decomp(E,HI(f,:));
            model.interactions.EV(f,:) = Eigen_values_inter';
            
            vh = nb_interactions(f) - 1; % df = q above
            s = min(vh,p); % subspace in which mean Ys are located
            clear nb_items; x = I{f};
            for c=1:nb_interactions(f)
                nb_items(c) = numel(find(x(:,c)));
            end
            
            if sum(nb_items == nb_items(1)) == length(nb_items)
                ve = nb_interactions(f)*(nb_items(1)-1);     % dfe equal sample sizes
            else
                ve = sum(nb_items) - nb_interactions(f);     % dfe different sample sizes
            end
            
            if s > 1
                m = (abs(vh-p)-1)/2;
                N = (ve-p-1) / 2;
                
                % Pillai
                V = sum(Eigen_values_inter ./ (1+Eigen_values_inter));
                df_interaction_Pillai(f) = s*(2*m+s+1);
                dfe_interaction_Pillai(f) = s*(2*N+s+1);
                F_interaction_Pillai(f) = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
                pval_interaction_Pillai(f) = 1-fcdf(F_interaction_Pillai(f),df_interaction_Pillai(f),dfe_interaction_Pillai(f));
                
                % Roy's test
                theta = max(Eigen_values_inter) / (1+max(Eigen_values_inter));
                df_interaction_Roy(f) = max(p,vh);
                dfe_interaction_Roy(f) = ve - 1; % in Renchner it is proposed to use ve - max(p,vh) -1 while in Statistica it is ve -1
                F_interaction_Roy(f) = (dfe_interaction_Roy(f)*max(Eigen_values_cond))/df_interaction_Roy(f);
                pval_interaction_Roy(f) = 1-fcdf(F_interaction_Roy(f), df_interaction_Roy(f), dfe_interaction_Roy(f));
                
            else % = only one non zeros Eigen value s = 1 and/or vh = 1
                
                V = sum(Eigen_values_inter ./ (1+Eigen_values_inter));
                U = max(Eigen_values_inter);
                theta = U;
                
                df_interaction_Pillai(f) = p; % number of electrodes
                dfe_interaction_Pillai(f) = ve-p+1;
                df_interaction_Roy(f) = df_interaction_Pillai(f);
                dfe_interaction_Roy(f) = dfe_interaction_Pillai(f);
                
                F_interaction_Pillai(f) = (dfe_interaction_Pillai(f)/df_interaction_Pillai(f)) * max(Eigen_values_cond);
                pval_interaction_Pillai(f) = 1-fcdf(F_interaction_Pillai(f),df_interaction_Pillai(f),dfe_interaction_Pillai(f));
                F_interaction_Roy(f) = F_interaction_Pillai(f);
                pval_interaction_Roy(f) = pval_interaction_Pillai(f);
            end
            
        end
    end
end


% -----------------------------------
%% compute F for continuous variables
% -----------------------------------

if nb_continuous ~=0
    
    if nb_factors == 0
        T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        R     = eye(size(Y,1)) - (X*pinv(X));
        E     = (Y'*R*Y);
        if strcmp(method,'OLS')
            W = ones(size(Y,1),1);
            Betas = X\Y; % numerically more stable than pinv
        elseif strcmp(method,'WLS')
            if isempty(W)
                [Betas,W] = limo_WLS(X,Y);
            else
                Betas = pinv(WX)*WY;
            end
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X,Y);
        end
        
        % compute model R^2
        % -----------------
        C = eye(size(X,2));
        C(:,size(X,2)) = 0;
        C0 = eye(size(X,2)) - C*pinv(C);
        X0 = X*C0;
        R0 = eye(size(Y,1)) - (X0*pinv(X0));
        M  = R0 - R;
        H  = (Betas'*X'*M*X*Betas);
        [Eigen_vectors_R2, Eigen_values_R2] = limo_decomp(E,H);
        
        % Generalized R2
        % variance covariance matrix
        S = cov([Y X(:,1:size(X,2)-1)]);
        Syy = S(1:size(Y,2),1:size(Y,2));
        Sxy = S(size(Y,2)+1:size(S,1),1:size(Y,2));
        Syx = S(1:size(Y,2),size(Y,2)+1:size(S,2));
        Sxx = S(size(Y,2)+1:size(S,1),size(Y,2)+1:size(S,2));
        Rsquare_multi = trace(Sxy*Syx) / sqrt(trace(Sxx.^2)*trace(Syy.^2)); % Robert and Escoufier, J.Royal Stat Soc, C - 1976
        
        [Eigen_vectors_R2, Eigen_values_R2] = limo_decomp(E,H);
        p = size(Y,2); % = number of variables (dimension)
        q = rank(X); % = number of regressors (df)
        s = min(p,q); % df
        n = size(Y,1); % nb of observations (dfe)
        m = (abs(q-p)-1)/2;
        N = (n-q-p-2)/2;
        d = max(p,q);
        
        theta = max(Eigen_values_R2) / (1+max(Eigen_values_R2)); % Roy
        R2_Roy_value = theta; % = 1st canonical correlation
        R2_Roy_F     = ((n-d-1)*max(Eigen_values_R2))/d;
        R2_Roy_p     = 1-fcdf(R2_Roy_F, d, (n-d-1));
        
        V = sum(Eigen_values_R2 ./ (1+Eigen_values_R2)); % Pillai
        R2_Pillai_value = V / s; % average of canonical correlations
        R2_Pillai_F     = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
        R2_Pillai_p     = 1-fcdf(R2_Pillai_F,(s*(2*m+s+1)),(s*(2*N+s+1)));
        
    else
        
        % compute
        nb_conditions = sum(nb_conditions) + sum(nb_interactions);
        model.continuous.EV = [];
        for n = 1:nb_continuous
            C    = zeros(size(X,2));
            C(nb_conditions+n,nb_conditions+n) = 1;
            C0   = eye(size(X,2)) - C*pinv(C);
            X0   = X*C0;
            R0   = eye(size(Y,1)) - (X0*pinv(X0));
            M    = R0 - R;
            H    = Betas'*X'*M*X*Betas;
            [Eigen_vectors_continuous, Eigen_values_continuous] = limo_decomp(E,H);
            model.continuous.EV = [model.continuous.EV Eigen_values_continuous'];
            
            df_continuous = size(Y,2);
            dfe_continuous = size(Y,1)-nb_continuous-size(Y,2);
            if nb_conditions ~= 0
                dfe_continuous = dfe_continuous - length(nb_conditions);
            end
            s = min(size(Y,2),size(X,2));
            m = (abs(size(X,2)-size(Y,2))-1) / 2;
            
            % Roy's test
            theta_continuous(n) = max(Eigen_values_continuous(:,n)) / (1+max(Eigen_values_continuous(:,n)));
            F_continuous_Roy(n) = (dfe_continuous*max(Eigen_values_continuous(:,n)))/df_continuous;
            pval_continuous_Roy(n) = 1-fcdf(F_continuous_Roy(n), df_continuous, dfe_continuous);
            
            % Pillai
            V_continuous(n) = sum(Eigen_values_continuous(:,n) ./ (1+Eigen_values_continuous(:,n)));
            F_continuous_Pillai(n) = F_continuous_Roy(n);
            pval_continuous_Pillai(n) = pval_continuous_Roy(n);
            
        end
    end
end

% ----------------------------
%% update the model structure
% ----------------------------

model.R2.V = Rsquare_multi;
model.R2.EV = Eigen_values_R2(1:s);
model.R2.Roy.F = R2_Roy_F;
model.R2.Roy.p = R2_Roy_p;
model.R2.Pillai.F = R2_Pillai_F;
model.R2.Pillai.p = R2_Pillai_p;
model.R2.varEV = Eigen_values_R2(1:s) ./ sum(Eigen_values_R2(1:s))*100;
model.betas = Betas;

if nb_conditions ~= 0
    model.conditions.EV           = [Eigen_values_cond'];
    model.conditions.Pillai.F     = F_conditions_Pillai;
    model.conditions.Pillai.p     = pval_conditions_Pillai;
    model.conditions.Pillai.df    = df_conditions_Pillai;
    model.conditions.Pillai.dfe   = dfe_conditions_Pillai;
    model.conditions.Roy.F        = F_conditions_Roy;
    model.conditions.Roy.p        = pval_conditions_Roy;
    model.conditions.Roy.df       = df_conditions_Roy;
    model.conditions.Roy.dfe      = dfe_conditions_Roy;
end

if nb_factors == 1
    model.conditions.varEV  = z_importance * 100; % to see in how many dimension the mean vectors lie in (guide using Roy or Pillai)
    
    % outcomes from Discriminant Analysis:
    model.Discriminant.D           = standardized_eigenvectors;
    model.Discriminant.scaledEV_coef = scaled_eigenvectors;
    model.Discriminant.importanceZ = z_importance * 100; % so also shows the relative importance of each discriminant function:
    model.Discriminant.cc_squared  = cc_squared;
    model.Discriminant.Z           = z;
    model.Discriminant.Z_cent     = centered_z;
    model.Discriminant.cc_squared  = cc_squared;
    
    % outcome from Classification:
    model.Classification.Linear.Acc       = training_Acc;
    model.Classification.Linear.cvAcc     = cvAcc;
    model.Classification.Linear.cvSD      = cvSD;
    
    model.Classification.Quadratic.Acc       = q_training_Acc;
    model.Classification.Quadratic.cvAcc     = q_cvAcc;
    model.Classification.Quadratic.cvSD      = q_cvSD;
    
end

if nb_interactions ~= 0
    model.interactions.Pillai.F     = F_interactions_Pillai;
    model.interactions.Pillai.p     = pval_interactions_Pillai;
    model.interactions.Pillai.df    = df_interactions_Pillai;
    model.interactions.Pillai.dfe   = dfe_interactions_Pillai;
    model.interactions.Roy.F        = F_interactions_Roy;
    model.interactions.Roy.p        = pval_interactions_Roy;
    model.interactions.Roy.df       = df_interactions_Roy;
    model.interactions.Roy.dfe      = dfe_interactions_Roy;
end

if nb_continuous > 0
    model.continuous.Pillai.F     = F_continuous_Pillai;
    model.continuous.Pillai.p     = pval_continuous_Pillai;
    model.continuous.Pillai.df    = df_continuous;
    model.continuous.Pillai.dfe   = dfe_continuous;
    model.continuous.Roy.F        = F_continuous_Roy;
    model.continuous.Roy.p        = pval_continuous_Roy;
    model.continuous.Roy.df       = df_continuous;
    model.continuous.Roy.dfe      = dfe_continuous;
end


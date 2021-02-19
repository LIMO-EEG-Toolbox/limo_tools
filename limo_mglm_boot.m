function model = limo_mglm_boot(varargin)

% Bootstrapped version of limo_mglm
% Importantly it also runs per time frame (under same resampling) - but do B bootstraps to obtain
% the distributon of F (and associated p values) under H0
% H0 is obtained by either by resampling from centered data (categorical designs)
% or sampling Y but leaving X intact, i.e. breaking the link between Y and X
% 
%
% FORMAT:
% model = limo_mglm(Y,LIMO, boot_table)
% model = limo_mglm(Y,X,nb_conditions,nb_interactions,nb_continuous,method,W, boot_table)
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
%   boot_table    = an optional argument - this is the resampling table
%                    if one calls limo_mglm_boot to loop throughout time frames,
%                    this might a good idea to provide such table so that
%                    the same resampling applies to each time frame.
% See also
% LIMO_DESIGN_MATRIX, LIMO_PCOUT, LIMO_IRLS, LIMO_EEG, LIMO_DECOMP
%
% Cyril Pernet & Iege Bassez v1 April 2018
% ----------------------------
% Copyright (C) LIMO Team 2018

%% varagin
nboot = 800; %

if nargin == 2 || nargin == 3
    y               = varargin{1};
    X               = varargin{2}.design.X;
    nb_conditions   = varargin{2}.design.nb_conditions;
    nb_interactions = varargin{2}.design.nb_interactions;
    nb_continuous   = varargin{2}.design.nb_continuous;
    method          = varargin{2}.design.method;
    cov_method      = varargin{2}.design.cov_method;
    
    if nargin == 2
        boot_table = randi(size(y,1), size(y,1), nboot); % resample observations (trials/subjects)
    elseif nargin == 3
        boot_table = varargin{3};
        nboot = size(boot_table,2);
    end
    
elseif nargin >= 7
    y               = varargin{1};
    X               = varargin{2};
    nb_conditions   = varargin{3};
    nb_interactions = varargin{4};
    nb_continuous   = varargin{5};
    method          = varargin{6};
    cov_method      = varargin{7};
    if nargin == 8
        boot_table  = varargin(8);
        nboot = size(boot_table); 
    elseif nargin == 7
        boot_table = randi(size(y,1), size(y,1),nboot); % resample observations (trials/subjects)
    end
    
else
    error('varargin error in limo_mglm_boot')
end

nb_factors = numel(nb_conditions);
if nb_factors == 1 && nb_conditions == 0
    nb_factors = 0;
end

% -----------
%% Data check
% -----------

if size(y,1)~=size(X,1)
    error('The number of events in Y and the design matrix are different')
end

if nb_interactions == 0
    nb_interactions = [];
end

% ----------
%% Bootstrap
% -----------
design = X;

% if categorical design, center data 1st
% ---------------------------------------
if nb_continuous == 0
    centered_y = NaN(size(y,1),size(y,2));
    if ~isempty(nb_interactions)
        % look up the last interaction to get unique groups
        if length(nb_interactions) == 1
            start_at = sum(nb_conditions);
        else
            start_at = sum(nb_conditions)+sum(nb_interactions(1:end-1));
        end
        
        for cel=(start_at+1):(start_at+nb_interactions(end))
            index = find(X(:,cel));
            centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),length(index),1);
        end
        
    elseif size(nb_conditions,2) == 1
        % no interactions because just 1 factor
        %figure; % as test everything went well with null distributions
        for cel=1:nb_conditions
            index = find(X(:,cel));
            centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),length(index),1);
            fprintf(['average centered y (for every electrode) for class' num2str(cel) ': ' num2str(round(mean(centered_y(index,:), 1),10)) '\n'])
            %subplot(1,nb_conditions,cel);histogram(centered_y(index,1)); title(['class ' num2str(cel)]); % quick check
        end
        
    else
        % create fake interaction to get groups
        [tmpX, interactions] = limo_make_interactions(X(:,1:(end-1)), nb_conditions);
        if length(interactions) == 1
            start_at = sum(nb_conditions);
        else
            start_at = sum(nb_conditions)+sum(interactions(1:end-1));
        end
        
        for cel=(start_at+1):(start_at+interactions(end))
            index = find(tmpX(:,cel));
            centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),[size(y(index,:),1)],1);
        end
        
    end
    
else
    centered_y = y;
    design = X;
    
end
clear y 

% workout the interaction increment (blocks on interactions terms to add-up)
if nb_factors > 1  && ~isempty(nb_interactions) % N-ways ANOVA with interactions
    for n=2:nb_factors
        increment(n-1) = size(nchoosek([1:nb_factors],n),1);
    end
else
    increment = [];
end

% compute for each bootstrap
% ---------------------------
parfor B = 1:nboot
    fprintf('boot n %g\n',B)
 
    % create data under H0
    if nb_continuous == 0
        % if just categorical variables, sample from the centered data and
        % the design simultaneously - rezscore if needed
        Y = centered_y(boot_table(:,B),:); % resample Y
        X = design(boot_table(:,B),:); % resample X
        
    else
        % sample and break the link between Y and X (regression and AnCOVA designs)
        Y = centered_y(boot_table(:,B),:); % resample (not really centerd y, see above)
        X = design; % stays the same
    end
    
    nb_items = NaN(1,nb_conditions);
    for c=1:nb_conditions
        nb_items(c) = sum(X(:,c));
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
            [Betas,W] = limo_WLS(X,Y);
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X,Y);
        end

        % compute model R^2
        % -----------------
        C = eye(size(X,2));
        C(:,size(X,2)) = 0;
        C0 = eye(size(X,2)) - C*pinv(C);
        X0 = X*C0;  % Reduced model
        R0 = eye(size(Y,1)) - (X0*pinv(X0));
        M  = R0 - R;  % Projection matrix onto Xc
        H = (Betas'*X'*M*X*Betas);  % SS Effect % only works
        % with rank deficient matrix (intercept column with ones as last column of X)
        if round(H,6) ~= round(T - E, 6) % if H is not equal to T - E problem!
            H = T - E; % temporary solution
        end 

        % Generalized R2
        % variance covariance matrix
        S = cov([Y X(:,1:size(X,2)-1)]);
        Syy = S(1:size(Y,2),1:size(Y,2));
        Sxy = S(size(Y,2)+1:size(S,1),1:size(Y,2));
        Syx = S(1:size(Y,2),size(Y,2)+1:size(S,2));
        Sxx = S(size(Y,2)+1:size(S,1),size(Y,2)+1:size(S,2));
        Rsquare_multi = trace(Sxy*Syx) / sqrt(trace(Sxx.^2)*trace(Syy.^2)); % Robert and Escoufier, J.Royal Stat Soc, C - 1976
        
        if strcmp(cov_method, 'regularized')
            [RegularizedCovariance, ~] = cov1para(Y);
            E = RegularizedCovariance .* (size(Y,1)-nb_conditions);
        end  
        [Eigen_vectors_R2, Eigen_values_R2] = limo_decomp(E,H,cov_method); 
        p = size(Y,2); % = number of variables (dimension)
        q = rank(X)-1; % -1 because of intercept column
        s = min(p,q); % df
        n = size(Y,1); % nb of observations (dfe)
        m = (abs(q-p)-1)/2;
        N = (n-q-p-2)/2;
        ve = n - rank(X);

        % Roy
        theta = max(Eigen_values_R2) / (1+max(Eigen_values_R2)); % Roy
        R2_Roy_value = theta; % = 1st canonical correlation
        R2_Roy_F     = ((ve-max(p,q)+q)*max(Eigen_values_R2))/max(p,q);
        R2_Roy_p     = 1-fcdf(R2_Roy_F, max(p,q), ve-max(p,q)+q);

        % Pillai
        V = sum(Eigen_values_R2 ./ (1+Eigen_values_R2)); % Pillai
        R2_Pillai_value = V / s; % average of canonical correlations
        R2_Pillai_F     = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
        R2_Pillai_p     = 1-fcdf(R2_Pillai_F,(s*(2*m+s+1)),(s*(2*N+s+1)));

        % compute F for categorical variables
        % -----------------------------------
        if nb_conditions ~= 0 && nb_continuous == 0
            Eigen_values_cond = Eigen_values_R2;
            Eigen_vectors_cond = Eigen_vectors_R2;

        elseif nb_conditions ~= 0 && nb_continuous ~= 0
            C = eye(size(X,2));
            C(:,(nb_conditions+1):size(X,2)) = 0;
            C0 = eye(size(X,2)) - C*pinv(C);
            X0 = X*C0; % Here the reduced model includes the covariates
            R0 = eye(size(Y,1)) - (X0*pinv(X0));
            M  = R0 - R;
            H  = (Betas'*X'*M*X*Betas); 
            [Eigen_vectors_cond, Eigen_values_cond] = limo_decomp(E,H, cov_method);
        end

        vh = nb_conditions - 1; % df = q above
        s = min(vh,p); % subspace in which mean Ys are located

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
            df_conditions_Pillai = s*(2*m+s+1);
            dfe_conditions_Pillai = s*(2*N+s+1);
            F_conditions_Pillai = ((2*N+s+1)*V) / ((2*m+s+1)*(s-V)');
            pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_conditions_Pillai,dfe_conditions_Pillai);

            % Roy's test
            theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond));
            df_conditions_Roy = max(p,vh);
            %dfe_conditions_Roy = ve - 1; % in Renchner it is proposed to use ve - max(p,vh) -1  while in Statistica it is ve -1. 
            dfe_conditions_Roy = ve - max(p,q) + q; % SAS site
            F_conditions_Roy = (dfe_conditions_Roy*max(Eigen_values_cond))/df_conditions_Roy;
            pval_conditions_Roy = 1-fcdf(F_conditions_Roy, df_conditions_Roy, dfe_conditions_Roy);

        else % = only one non zeros Eigen value s = 1 and/or vh = 1 (see p. 169 Rencher)

            theta = max(Eigen_values_cond) / (1+max(Eigen_values_cond)); % Roy
            V = theta; % Pillai(1) equals theta

            df_conditions_Pillai = p; 
            dfe_conditions_Pillai = ve-p+1;
            df_conditions_Roy = df_conditions_Pillai;
            dfe_conditions_Roy = dfe_conditions_Pillai;

            F_conditions_Pillai = (dfe_conditions_Pillai/df_conditions_Pillai) * max(Eigen_values_cond);
            pval_conditions_Pillai = 1-fcdf(F_conditions_Pillai,df_conditions_Pillai,dfe_conditions_Pillai);
            F_conditions_Roy = F_conditions_Pillai;
            pval_conditions_Roy = pval_conditions_Pillai;
        end

        %% 
        % Discriminant Analysis
        % ----------------------------------------
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

%             % get the discriminant function(s) 
%             scaled_eigenvectors = scaled_eigenvectors(:,1:s); % s: nb nonzero eigenvalues
%             Eigen_values_cond = Eigen_values_cond(1:s); % s: nb nonzero eigenvalues
%             z = Y * scaled_eigenvectors; % the discriminant functions themself 
%             z_importance = Eigen_values_cond ./ sum(Eigen_values_cond); % variance corresponding to each eigenvalue
%             cc_squared = Eigen_values_cond ./ (1 + Eigen_values_cond); % canonical correlation
%             centered_Y = bsxfun(@minus, Y,mean(Y));     
%             centered_z = centered_Y * scaled_eigenvectors; % centered discriminant functions

            % now, get standardized eigenvectors
            Spooled = E./ve;
%             standardized_eigenvectors = bsxfun(@times, scaled_eigenvectors, sqrt(diag(Spooled)));                
        end      
        %%
        % do the classification 
        %--------------------------------------------------------
        % get training linear decoding accuracy:
        [class,~] = find(X(:,1:nb_conditions)');
        [predicted] = limo_LDA(Y, class, Y, cov_method);
        confmat = confusionmat(class, predicted); 
        training_Acc = sum(diag(confmat/sum(sum(confmat))));

        % get 10-fold CV linear decoding accuracy:
        folds = 10;
        cvp = cvpartition(class,'k',folds); 
        accuracyvector = NaN(1,folds);
        for foldi=1:folds
            trainIdx = cvp.training(foldi);
            testIdx = cvp.test(foldi);
            % make a training and test set based on indexes:
            trainingset = Y(trainIdx,:);
            testset = Y(testIdx,:);
            % class labels:
            traininglabel = class(trainIdx,1);
            testlabel = class(testIdx,1);
            % Classification:
            predicted = limo_LDA(trainingset,traininglabel,testset, cov_method);
            confmat = confusionmat(testlabel, predicted); 
            accuracyvector(foldi) = sum(diag(confmat/sum(sum(confmat))));
        end
        cvAcc  = mean(accuracyvector);
        cvSD = std(accuracyvector);
        
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
    % to be added
        % ------------------------------------------------
    elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways MANOVA with interactions
        % ------------------------------------------------
    % to be added
    end    

    % ----------------------------
    %% update the model structure
    % ----------------------------
    BETASB(:,:,B) = Betas';
    MODELR2{B} = Rsquare_multi;
    
    MODELF_Pillai{B} = R2_Pillai_F;
    MODELp_Pillai{B} = R2_Pillai_p;
    F_CONDVALUES_Pillai{B}  = F_conditions_Pillai;
    p_CONDVALUES_Pillai{B}  = pval_conditions_Pillai;    

    MODELF_Roy{B} = R2_Roy_F;
    MODELp_Roy{B} = R2_Roy_p;
    F_CONDVALUES_Roy{B}  = F_conditions_Roy;
    p_CONDVALUES_Roy{B}  = pval_conditions_Roy;  
    
    Linear_Classification{B} = training_Acc; %overestimate but good under H0
    Quadratic_Classification{B} = q_training_Acc; %overestimate but good under H0
    
end

model.R2 = MODELR2;
model.F.Pillai = MODELF_Pillai;
model.p.Pillai = MODELp_Pillai;
model.F.Roy = MODELF_Roy;
model.p.Roy = MODELp_Roy;
model.Betas = BETASB;

try
    %MANOVA
    model.conditions.F.Pillai = F_CONDVALUES_Pillai;
    model.conditions.p.Pillai  = p_CONDVALUES_Pillai;
    model.conditions.F.Roy = F_CONDVALUES_Roy;
    model.conditions.p.Roy  = p_CONDVALUES_Roy;
    %DISCRIMINANT - is descriptive
    %CLASSIFICATION
    model.Classification.Linear = Linear_Classification;
    model.Classification.Quadratic = Quadratic_Classification;    
end
end

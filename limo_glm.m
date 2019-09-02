function model = limo_glm(varargin)

% General Linear model analysis of EEG data.
% The model consider trials/subjects as independent observations and analyses
% are performed at one channel/IC/source only, but for all frames (time and/or
% freq) points. All multiple comparisions correction are performed post-hoc
% using bootstrap under H0.
%
% see the <https://www.hindawi.com/journals/cin/2011/831409/ LIMO EEG paper>
%
% FORMATS: model = limo_glm(Y,LIMO)
%          model = limo_glm(Y, X, nb_conditions, nb_interactions, ...
%                  nb_continuous, method, Analysis type, n_freqs, n_times)
% INPUTS:
%   Y               = 2D matrix of EEG data with format trials x frames
%   LIMO            = structure that contains the above information (except Y)
%                     or input info that are in LIMO.mat
%   X               = 2 dimensional design matrix
%   nb_conditions   = a vector indicating the number of conditions per factor
%   nb_interactions = a vector indicating number of columns per interactions
%   nb_continuous   = number of covariates
%   method          = 'OLS', 'WLS', 'IRLS' (bisquare)
%   analysis type   = 'Time', 'Frequency' or 'Time-Frequency'
%   n_freqs         = the nb of frequency bins
%   n_times         = the nb of time bins
%
% OUTPUTS:
%   model.R2_univariate = the R2 of the model
%   model.F             = the F value of the model
%   model.df            = the df associated to the model F
%   model.p             = the p value of the model
%   model.betas         = the beta parameters (dimensions nb of paramters x frames)
%   model.W             = the weights for ea ch trial/subject (and frames if IRLS)
%   model.conditions    = main categorical effects
%          --> F/p in rows are the factors, in columns time frames
%          --> df row 1 = df, row 2 = dfe, columns are factors
%   model.interactions  = interaction effects
%          --> F/p in rows are the factors, in columns time frames
%          --> df row 1 = df, row 2 = dfe, columns are interaction levels
%   model.continuous = continuous effects
%          --> F/p in rows are the variables, in columns time frames
%          --> df column 1 = df, column2 2 = dfe (same for all covariates)
%
% NOTES:
%
% - The parameters can be computed using 3 methods: ordinary least squares
% (OLS), weighted least squares (WLS) which attribute unique weights per trial
% and iterative reweighted least squares (IRLS) which attribute weights for
% each observations
% - For WLS, since this weights are derived for the entire trial, there is a
% depencies similar to a multivariate approach, for which we have no analytical
% solution and the p-values must then be derived post-hoc via bootstrap.
% see < WLS paper>
% - For IRLS, dfe are computed using the Satterthwaite approximation and
% a MSE correction is used
% - Each effect is accounted for given the other effects; this means that one
% can have a different number of trials per conditions/factors
% provided there is no interactions. For interaction models (typically for
% single subject anayses), this is not directly possible and no correction
% is provided. Instead, with a design created using limo_design_matrix, data
% would have been sampled to make sure the number of trials or subjects is
% identical across interaction terms.
% - For time*frequency analyses, limo_eeg_tf sends a vector of length freq*time
% that is unwrapped after running limo_glm. Weights cannot be computed that
% way and thus one calls LIMO.Analysis to unwrap and then rewarp to get the
% right weights and betas.
%
% References
% Christensen, R. 2002. Plane answers to complex questions. 3rd Ed. Springer-Verlag
% Friston et al. 2007. Statitical Parametric Mapping. Academic Press
% Yandell, B.S. 1997. Practical Data Analysis For Designed Experiments. Chapman & Hall
%
% See also
% LIMO_GLM_HANDLING, LIMO_DESIGN_MATRIX, LIMO_WLS, LIMO_IRLS
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% varagin

if nargin == 2
    Y               = varargin{1};
    X               = varargin{2}.design.X;
    nb_conditions   = varargin{2}.design.nb_conditions;
    nb_interactions = varargin{2}.design.nb_interactions;
    nb_continuous   = varargin{2}.design.nb_continuous;
    method          = varargin{2}.design.method;
    Analysis        = varargin{2}.Analysis;
    if strcmp(Analysis,'Time-Frequency')
        if strcmp(method,'WLS')
            method = 'WLS-TF'; % run weights per freq band
        end
        n_freqs = varargin{2}.data.size4D(2);
        n_times = varargin{2}.data.size4D(3);
    end
    clear varargin
elseif nargin == 9
    Y               = varargin{1};
    X               = varargin{2};
    nb_conditions   = varargin{3};
    nb_interactions = varargin{4};
    nb_continuous   = varargin{5};
    method          = varargin{6};
    Analysis        = varargin{7};
    n_freqs         = varargin{8};
    n_times         = varargin{9};
else
    error('varargin error')
end

if isempty(nb_conditions);   nb_conditions = 0; end
if isempty(nb_interactions); nb_interactions = 0; end
if isempty(nb_continuous);   nb_continuous = 0; end
    
nb_factors = numel(nb_conditions);
if nb_factors == 1 && nb_conditions == 0
    nb_factors = 0;
end

% -----------
%% Data check
% -----------
if ~isreal(Y)
    Y = abs(Y);
end

if size(Y,1)~=size(X,1)
    error('The number of events in Y and the design matrix are different')
end

if nb_interactions == 0
    nb_interactions = [];
end

%% Compute model parameters
% ------------------------------

% compute Beta parameters and weights
% -----------------------------------
if strcmp(method,'OLS')
    
    if strcmp(Analysis,'Time-Frequency')
        W = ones(n_freqs,size(X,1));
    else
        W = ones(size(Y,1),1);
    end
    WX = X;
    
    if nb_continuous ~=0 && nb_factors == 0
        Betas = WX\Y; % numerically more stable than pinv
    else
        Betas = pinv(WX)*Y;
    end
    
elseif strcmp(method,'WLS')
    [Betas,W] = limo_WLS(X,Y);
    WX        = [X(:,1:end-1).*repmat(W,1,size(X,2)-1) X(:,end)];
    
elseif strcmp(method,'WLS-TF')
    % unpack the data
    [n_freq_times, N] = size(Y');
    if n_freq_times ~= n_freqs*n_times
        error('dimensions disagreement to reshape freq*time')
    else
        reshaped = nan(n_freqs, n_times, N);
    end
    
    for tr = 1:N
        eft_3d = NaN(n_freqs,n_times);
        for tm = 1:n_times
            this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 2D tf
            eft_3d(:,tm) =Y(tr,this_freq_start_index:(this_freq_start_index+n_freqs-1))';
        end
        reshaped(:,:,tr) = eft_3d; clear eft_3d
    end
    
    % get estimates per frequency band
    index1 = 1;
    Betas  = NaN(size(X,2),n_freqs*n_times);
    W      = NaN(n_freqs,size(X,1));
    for f=1:n_freqs
        [Betas(:,index1:6:(n_freqs*n_times)),W(f,:)] = limo_WLS(X,squeeze(reshaped(f,:,:))');
        index1=index1+1;
    end
    clear reshaped
    WX = X .* repmat(W,1,size(X,2));
    
elseif strcmp(method,'IRLS')
    [Betas,W] = limo_IRLS(X,Y);
    % WX = X.*W; per frame! switch method
end

%% ------------------------------------

switch method
    
    case {'OLS','WLS'}
        % -----------------------------------------------------------------
        
        %% Compute model statistics
        % ------------------------------
        % total sum of squares, projection matrix for errors, residuals
        % --------------------------------------------------------------
        T   = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total (the data)
        % HM  = WX*pinv(WX);                                                     % Hat matrix, projection onto X
        R   = eye(size(Y,1)) - WX*pinv(WX);                                      % Projection onto E
        E   = Y'*R*Y;                                                            % SS Error
        
        % degrees of freedom
        % -------------------
        df = rank(WX)-1;
        if strcmp(method,'OLS')
            dfe = size(Y,1)-rank(WX);
        else
            %             Satterthwaite approximation: dfe = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM));
            %             Multivariate df: dfe = size(Y,1)-size(Y,2)+rank(WX);
            %             Cheverud 2001: EV  = eig(corr((R*Y)));
            %                            M   = length(EV);
            %                            V   = sum((EV-1).^2) / (M-1);
            %                            dfe = (1 + (M-1)*(1-V/M)) - df;
            %           Li and Ji 2005
            EV  = abs(eig(corr(R*Y)));
            x   = single(EV>=1) + (EV - floor(EV));
            dfe = size(Y,1) - sum(x) + rank(WX) + 1;
        end
        
        % model R^2
        % -----------
        C              = eye(size(X,2));
        C(:,size(X,2)) = 0;                              % all columns but the constant
        C0             = eye(size(X,2)) - C*pinv(C);     % only the constant
        X0             = WX*C0;                          % Reduced model design matrix
        R0             = eye(size(Y,1)) - (X0*pinv(X0)); % Projection onto error
        M              = R0 - R;                         % Projection matrix onto Xc
        H              = (Betas'*X'*M*X*Betas);          % SS Effects
        Rsquare        = diag(H)./diag(T);               % Variance explained
        F_Rsquare      = (diag(H)./df) ./ (diag(E)/dfe);
        p_Rsquare      = 1 - fcdf(F_Rsquare, df, dfe);
        
        % update the model structure
        % ----------------------------
        
        model.W                 = W;
        model.betas             = Betas;
        model.betas_se          = Betas;
        for b=1:size(Y,2)
            model.betas_se(:,b) = sqrt(diag((E(b,b)/dfe)*pinv(WX'*WX)));
        end
        model.R2_univariate     = Rsquare;
        model.F                 = F_Rsquare;
        model.df                = [df dfe];
        if strcmp(method,'WLS')
            model.p             = NaN(size(p_Rsquare)); % p values aren't valid in this scheme
        else
            model.p             = p_Rsquare;
        end
        
        %% Compute effects
        % ------------------
        
        % ---------------------------------
        if nb_factors == 1   %  1-way ANOVA
            % ---------------------------------
            
            % compute F for categorical variables
            % -----------------------------------
            if nb_conditions ~= 0 && nb_continuous == 0
                df_conditions                    = df;
                F_conditions                     = F_Rsquare;
                pval_conditions                  = p_Rsquare;
                
            elseif nb_conditions ~= 0 && nb_continuous ~= 0
                C                                = eye(size(X,2));
                C(:,(nb_conditions+1):size(X,2)) = 0;
                C0                               = eye(size(X,2)) - C*pinv(C);
                X0                               = WX*C0; % here the reduced model includes the covariates
                R0                               = eye(size(Y,1)) - (X0*pinv(X0));
                M                                = R0 - R; % hat matrix for all categorical regressors (1 factor)
                H                                = (Betas'*X'*M*X*Betas);
                df_conditions                    = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C)-1 if OLS; same as tr(M)?
                F_conditions                     = (diag(H)/df_conditions) ./ (diag(E)/dfe);
                pval_conditions                  = 1 - fcdf(F_conditions(:), df_conditions, dfe);
            end
            
            model.conditions.F                   = F_conditions;
            model.conditions.df                  = [df_conditions ; dfe];
            if strcmp(method,'WLS')
                model.conditions.p               = NaN(size(pval_conditions)); % p values aren't valid in this scheme
            else
                model.conditions.p               = pval_conditions;
            end
            
            
            % ------------------------------------------------
        elseif nb_factors > 1  && isempty(nb_interactions) % N-ways ANOVA without interactions
            % ------------------------------------------------
            
            % --------------------------------------
            % compute F and p values of each factor
            % --------------------------------------
            
            df_conditions            = zeros(1,length(nb_conditions));
            F_conditions             = zeros(length(nb_conditions),size(Y,2));
            pval_conditions          = zeros(length(nb_conditions),size(Y,2));
            
            % define the effect of interest (eoi)
            eoi                      = zeros(1,size(X,2));
            eoi(1:nb_conditions(1))  = 1:nb_conditions(1);
            eoni                     = 1:size(X,2);
            eoni                     = find(eoni - eoi);
            
            for f = 1:length(nb_conditions)
                C                    = eye(size(X,2));
                C(:,eoni)            = 0; % set all but factor of interest to 0
                C0                   = eye(size(X,2)) - C*pinv(C);
                X0                   = WX*C0; % the reduced model include all but the factor f
                R0                   = eye(size(Y,1)) - (X0*pinv(X0));
                M                    = R0 - R; % hat matrix for factor f
                H                    = (Betas'*X'*M*X*Betas);
                df_conditions(f)     = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C)-1 if OLS;
                F_conditions(f,:)    = (diag(H)/df_conditions(f)) ./ (diag(E)/dfe);
                pval_conditions(f,:) = 1 - fcdf(F_conditions(f,:), df_conditions(f), dfe);
                
                % update factors
                if f<length(nb_conditions)
                    update           = find(eoi,1,'last'); % max(find(eoi));
                    eoi              = zeros(1,size(X,2));
                    eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
                    eoni             = 1:size(X,2);
                    eoni             = find(eoni - eoi);
                end
            end
            
            model.conditions.F      = F_conditions;
            model.conditions.df     = [df_conditions ; repmat(dfe,1,numel(df_conditions))]';
            if strcmp(method,'WLS')
                model.conditions.p  = NaN(size(pval_conditions)); % p values aren't valid in this scheme
            else
                model.conditions.p  = pval_conditions;
            end
            
            % ------------------------------------------------
        elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways ANOVA with interactions
            % ------------------------------------------------
            
            % ---------------------------------------------------
            % start by ANOVA without interaction for main effects
            % ---------------------------------------------------
            
            H                 = NaN(length(nb_conditions),size(Y,2));
            df_conditions     = NaN(1,length(nb_conditions));
            F_conditions      = NaN(length(nb_conditions),size(Y,2));
            pval_conditions   = NaN(length(nb_conditions),size(Y,2));
            HI                = NaN(length(nb_interactions),size(Y,2));
            df_interactions   = NaN(1,length(nb_interactions));
            F_interactions    = NaN(length(nb_interactions),size(Y,2));
            pval_interactions = NaN(length(nb_interactions),size(Y,2));
            
            % covariates
            covariate_columns = (sum(nb_conditions)+sum(nb_interactions)+1):(size(X,2)-1);
            
            % main effects
            dummy_columns = 1:sum(nb_conditions);
            
            % re-define X for main effects
            x = [X(:,dummy_columns) X(:,covariate_columns) ones(size(X,1),1)];
            
            % run same model as above with re-defined model x and
            % using the weights from the full model
            wx                       = x.*repmat(W,1,size(x,2));
            betas                    = pinv(wx)*(Y.*repmat(W,1,size(Y,2)));
            R                        = eye(size(Y,1)) - wx*pinv(wx);
            eoi                      = zeros(1,size(x,2));
            eoi(1:nb_conditions(1))  = 1:nb_conditions(1);
            eoni                     = 1:size(x,2);
            eoni                     = find(eoni - eoi);
            
            for f = 1:length(nb_conditions)
                C                    = eye(size(x,2));
                C(:,eoni)            = 0;
                C0                   = eye(size(x,2)) - C*pinv(C);
                X0                   = wx*C0;
                R0                   = eye(size(Y,1)) - (X0*pinv(X0));
                M                    = R0 - R;
                H(f,:)               = diag((betas'*x'*M*x*betas));
                df_conditions(f)     = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C)-1 if OLS;
                F_conditions(f,:)    = (H(f,:)./df_conditions(f)) ./ (diag(E)./dfe)'; % note dfe from full model
                pval_conditions(f,:) = 1 - fcdf(F_conditions(f,:), df_conditions(f), dfe);
                
                % update factors
                if f<length(nb_conditions)
                    update           = find(eoi,1,'last'); % max(find(eoi));
                    eoi              = zeros(1,size(x,2));
                    eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
                    eoni             = 1:size(x,2);
                    eoni             = find(eoni - eoi);
                end
            end
            
            model.conditions.F       = F_conditions;
            model.conditions.df      = [df_conditions ; repmat(dfe,1,numel(df_conditions))]';
            if strcmp(method,'WLS')
                model.conditions.p   = NaN(size(pval_conditions)); % p values aren't valid in this scheme
            else
                model.conditions.p   = pval_conditions;
            end
            
            % ---------------------------
            % now deal with interactions
            % ---------------------------
            
            if nb_factors == 2 && nb_continuous == 0 % the quick way with only one interaction
                HI                 = diag(T)' - H(1,:) - H(2,:) - diag(E)';
                df_interactions    = prod(df_conditions);
                F_interactions     = (HI./df_interactions) ./ (diag(E)/dfe)';
                pval_interactions  = 1 - fcdf(F_interactions, df_interactions, dfe);
                
            else % run through each interaction
                
                % part of X unchanged
                Main_effects = X(:,dummy_columns);
                Cov_and_Mean = [X(:,covariate_columns) ones(size(Y,1),1)];
                
                % check interaction level sizes in X
                index = 1;
                Ifactors = NaN(1,length(nb_interactions));
                interaction = cell(1,length(nb_interactions));
                for n=2:nb_factors
                    combinations = nchoosek(1:nb_factors,n); % note it matches X below because computed the same way in limo_design_matrix
                    for c = 1:size(combinations,1)
                        Ifactors(index) = length(combinations(c,:));
                        interaction{index} = combinations(c,:);
                        index = index + 1;
                    end
                end
                
                % loop through interactions
                % substituting and/or incrementing parts of X
                Istart     = size(Main_effects,2)+1; % where we start interaction in X
                Ilowbound  = size(Main_effects,2)+1;
                for f=1:length(nb_interactions)
                    I              = X(:,Istart:(Istart+nb_interactions(f)-1));
                    if length(interaction{f}) == 2 % 1st oder interaction is main + I
                        x          = [Main_effects I Cov_and_Mean];
                    else % higher oder inteaction includes lower levels
                        Isize      = sum(nb_interactions(1:find(Ifactors == Ifactors(f),1) - 1));
                        Ihighbound = size(Main_effects,2)+Isize;
                        x          = [Main_effects X(:,Ilowbound:Ihighbound) I Cov_and_Mean];
                    end
                    eoibound = size(x,2) - size(I,2) - size(Cov_and_Mean,2);
                    
                    % run same model as above
                    wx                     = x.*repmat(W,1,size(x,2));
                    betas                  = pinv(wx)*(Y.*repmat(W,1,size(Y,2)));
                    R                      = eye(size(Y,1)) - (wx*pinv(wx));
                    eoi                    = zeros(1,size(x,2));
                    eoi(eoibound+1:(eoibound+nb_interactions(f))) = eoibound+1:(eoibound+nb_interactions(f));
                    eoni                   = 1:size(x,2);
                    eoni                   = find(eoni - eoi);
                    C                      = eye(size(x,2));
                    C(:,eoni)              = 0;
                    C0                     = eye(size(x,2)) - C*pinv(C);
                    X0                     = wx*C0;
                    R0                     = eye(size(Y,1)) - (X0*pinv(X0));
                    M                      = R0 - R;
                    HI(f,:)                = diag((betas'*x'*M*x*betas))';
                    df_interactions(f)     = prod(df_conditions(interaction{f}));
                    F_interactions(f,:)    = (HI(f,:)./df_interactions(f)) ./ (diag(E)/dfe)';
                    pval_interactions(f,:) = 1 - fcdf(F_interactions(f,:), df_interactions(f), dfe);
                    Istart                 = Istart+nb_interactions(f);
                end
            end
            
            model.interactions.F       = F_interactions;
            model.interactions.df      = [df_interactions ; repmat(dfe,1,numel(df_interactions))]';
            if strcmp(method,'WLS')
                model.interactions.p   = NaN(size(pval_interactions)); % p values aren't valid in this scheme
            else
                model.interactions.p   = pval_interactions;
            end
        end
        
        % -----------------------------------
        %% compute F for continuous variables
        % -----------------------------------
        
        if nb_continuous ~=0
            
            if nb_factors == 0 && nb_continuous == 1 % simple regression
                model.continuous.F  = F_Rsquare;
                model.continuous.df = [1 (size(Y,1)-rank(X))];
                if strcmp(method,'WLS')
                    model.continuous.p         = NaN(size(p_Rsquare)); % p values aren't valid in this scheme
                else
                    model.continuous.p       = p_Rsquare;
                end
                
            else % ANCOVA type of designs
                
                % pre-allocate space
                df_continuous   = zeros(nb_continuous,size(Y,2));
                F_continuous    = zeros(nb_continuous,size(Y,2));
                pval_continuous = zeros(nb_continuous,size(Y,2));
                
                % compute
                N_conditions = sum(nb_conditions) + sum(nb_interactions);
                for n = 1:nb_continuous
                    C                                = zeros(size(X,2));
                    C(N_conditions+n,N_conditions+n) = 1; % pick up one regressor at a time
                    C0                               = eye(size(X,2)) - C*pinv(C);
                    X0                               = WX*C0; % all but rehressor of interest
                    R0                               = eye(size(Y,1)) - (X0*pinv(X0));
                    M                                = R0 - R; % hat matrix for regressor of interest
                    H                                = Betas'*X'*M*X*Betas;
                    df_continuous(n)                 = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C) if OLS;
                    F_continuous(n,:)                = (diag(H)./(df_continuous(n))) ./ (diag(E)/dfe);
                    pval_continuous(n,:)             = 1 - fcdf(F_continuous(n,:), 1, dfe); % dfe same as size(Y,1)-rank(X) if OLS
                end
                
                model.continuous.F                   = F_continuous';
                model.continuous.df                  = [1 dfe];
                if strcmp(method,'WLS')
                    model.continuous.p               = NaN(size(pval_continuous')); % p values aren't valid in this scheme
                else
                    model.continuous.p               = pval_continuous';
                end
            end
        end
        
        
        % ---------------------------------------------------------------------
    case 'IRLS'
        % -----------------------------------------------------------------
        
        % pre-allocate memory space
        Rsquare   = NaN(1,size(Y,2));
        F_Rsquare = NaN(1,size(Y,2));
        p_Rsquare = NaN(1,size(Y,2));
        betas_se  = NaN(size(Betas,1),size(Y,2));
        df        = NaN(1,size(Y,2));
        dfe       = NaN(1,size(Y,2));
        dof       = NaN(2,size(Y,2));
        
        if nb_factors ~=0
            F_conditions    = NaN(length(nb_conditions),size(Y,2));
            pval_conditions = NaN(length(nb_conditions),size(Y,2));
            df_conditions   = NaN(length(nb_conditions),size(Y,2));
        end
        
        if nb_interactions ~=0
            HI                = NaN(length(nb_interactions),size(Y,2));
            F_interactions    = NaN(length(nb_interactions),size(Y,2));
            pval_interactions = NaN(length(nb_interactions),size(Y,2));
            df_interactions   = NaN(length(nb_interactions),size(Y,2));
            
            % check interaction level sizes in X
            index = 1;
            Ifactors = NaN(1,length(nb_interactions));
            interaction = cell(1,length(nb_interactions));
            for n=2:nb_factors
                combinations = nchoosek(1:nb_factors,n); % note it matches X below because computed the same way in limo_design_matrix
                for c = 1:size(combinations,1)
                    Ifactors(index) = length(combinations(c,:));
                    interaction{index} = combinations(c,:);
                    index = index + 1;
                end
            end
        end
        
        if nb_continuous >1
            F_continuous    = NaN(nb_continuous,size(Y,2));
            pval_continuous = NaN(nb_continuous,size(Y,2));
            df_continuous   = NaN(nb_continuous,size(Y,2));
        end
        
        % start computing
        T   = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        for frame = 1:size(Y,2)
            % model stats
            % -------------------------------------------------------------
            WX                     = X.*W(:,frame);
            HM                     = WX*pinv(WX);
            R                      = eye(size(Y,1)) - WX*pinv(WX);
            E                      = Y(:,frame)'*R*Y(:,frame);
            % The number of degrees of freedom can be defined as the minimum number of
            % independent coordinates that can specify the position of the system completely.
            % This gives the same as [rank(X)-1 (size(Y,1)-rank(X))] if OLS, here we
            % use the Satterthwaite approximation
            df(frame)              = trace(HM'*HM)^2/trace((HM'*HM)*(HM'*HM))-1;
            dfe(frame)             = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM));
            R_ols                  = eye(size(Y,1)) - X*pinv(X);
            E_ols                  = Y(:,frame)'*R_ols*Y(:,frame);
            % MSE adjustment, E cannot be smaller than OLS since the
            % hyperplane we fit is farther away from some observations
            if E < E_ols
                n = size(X,1); p = rank(X);
                sigmar = E/(n-p); sigmals = E_ols/(n-p);
                MSE = (n*sigmar + p^2*sigmals) / (n+p^2);
                E = MSE * dfe(frame);
            end
            C                 = eye(size(X,2));
            C(:,size(X,2))    = 0;
            C0                = eye(size(X,2)) - C*pinv(C);
            X0                = WX*C0;
            R0                = eye(size(Y,1)) - (X0*pinv(X0));
            M                 = R0 - R;
            H                 = (Betas(:,frame)'*X'*M*X*Betas(:,frame));
            Rsquare(frame)    = H./T(frame,frame);
            F_Rsquare(frame)  = (H/df(frame))/(E/dfe(frame));
            p_Rsquare(frame)  = 1 - fcdf(F_Rsquare(frame), df(frame), dfe(frame));
            betas_se(:,frame) = sqrt(diag((E/dfe(frame))*pinv(WX'*WX)));
            dof(:,frame)      = [df(frame) dfe(frame)];
            
            %% Compute effects
            % ------------------
            
            % ---------------------------------
            if nb_factors == 1   %  1-way ANOVA
                % ---------------------------------
                
                % compute F for categorical variables
                % -----------------------------------
                if nb_conditions ~= 0 && nb_continuous == 0
                    df_conditions(frame)             = df(frame);
                    F_conditions(frame)              = F_Rsquare(frame);
                    pval_conditions(frame)           = p_Rsquare(frame);
                    
                elseif nb_conditions ~= 0 && nb_continuous ~= 0
                    C                                = eye(size(X,2));
                    C(:,(nb_conditions+1):size(X,2)) = 0;
                    C0                               = eye(size(X,2)) - C*pinv(C);
                    X0                               = WX*C0; % here the reduced model includes the covariates
                    R0                               = eye(size(Y,1)) - (X0*pinv(X0));
                    M                                = R0 - R; % hat matrix for all categorical regressors (1 factor)
                    H                                = (Betas(:,frame)'*X'*M*X*Betas(:,frame));
                    df_conditions(frame)             = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C)-1 if OLS; same as tr(M)?
                    F_conditions(frame)              = (H/df_conditions(frame)) ./ (E/dfe(frame));
                    pval_conditions(frame)           = 1 - fcdf(F_conditions(frame), df_conditions(frame), dfe(frame));
                end
                
                % ------------------------------------------------
            elseif nb_factors > 1  && isempty(nb_interactions) % N-ways ANOVA without interactions
                % ------------------------------------------------
                
                % --------------------------------------
                % compute F and p values of each factor
                % --------------------------------------
                
                % define the effect of interest (eoi)
                eoi                      = zeros(1,size(X,2));
                eoi(1:nb_conditions(1))  = 1:nb_conditions(1);
                eoni                     = 1:size(X,2);
                eoni                     = find(eoni - eoi);
                
                for f = 1:length(nb_conditions)
                    C                        = eye(size(X,2));
                    C(:,eoni)                = 0; % set all but factor of interest to 0
                    C0                       = eye(size(X,2)) - C*pinv(C);
                    X0                       = WX*C0; % the reduced model include all but the factor f
                    R0                       = eye(size(Y,1)) - (X0*pinv(X0));
                    M                        = R0 - R; % hat matrix for factor f
                    H                        = (Betas(:,frame)'*X'*M*X*Betas(:,frame));
                    df_conditions(f,frame)   = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C)-1 if OLS;
                    F_conditions(f,frame)    = (H/df_conditions(f,frame)) ./ (E/dfe(frame));
                    pval_conditions(f,frame) = 1 - fcdf(F_conditions(f,frame), df_conditions(f,frame), dfe(frame));
                    
                    % update factors
                    if f<length(nb_conditions)
                        update           = find(eoi,1,'last'); % max(find(eoi));
                        eoi              = zeros(1,size(X,2));
                        eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
                        eoni             = 1:size(X,2);
                        eoni             = find(eoni - eoi);
                    end
                end
                
                % ------------------------------------------------
            elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways ANOVA with interactions
                % ------------------------------------------------
                
                % ---------------------------------------------------
                % start by ANOVA without interaction for main effects
                % ---------------------------------------------------
                
                % covariates
                covariate_columns = (sum(nb_conditions)+sum(nb_interactions)+1):(size(X,2)-1);
                
                % main effects
                dummy_columns = 1:sum(nb_conditions);
                
                % re-define X
                x = [X(:,dummy_columns) X(:,covariate_columns) ones(size(X,1),1)];
                
                % run same model as above with re-defined model x and
                % using the weights from the full model
                wx                       = x.*repmat(W(:,frame),1,size(x,2));
                betas                    = pinv(wx)*(Y(:,frame).*W(:,frame));
                R                        = eye(size(Y,1)) - wx*pinv(wx);
                eoi                      = zeros(1,size(x,2));
                eoi(1:nb_conditions(1))  = 1:nb_conditions(1);
                eoni                     = 1:size(x,2);
                eoni                     = find(eoni - eoi);
                
                for f = 1:length(nb_conditions)
                    C                        = eye(size(x,2));
                    C(:,eoni)                = 0;
                    C0                       = eye(size(x,2)) - C*pinv(C);
                    X0                       = wx*C0;
                    R0                       = eye(size(Y,1)) - (X0*pinv(X0));
                    M                        = R0 - R;
                    H(f,frame)               = betas'*x'*M*x*betas;
                    df_conditions(f,frame)   = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C)-1 if OLS;
                    F_conditions(f,frame)    = (H(f,frame)./df_conditions(f,frame)) ./ (E./dfe(frame))'; % note dfe from full model
                    pval_conditions(f,frame) = 1 - fcdf(F_conditions(f,frame), df_conditions(f,frame), dfe(frame));
                    
                    % update factors
                    if f<length(nb_conditions)
                        update           = find(eoi,1,'last'); % max(find(eoi));
                        eoi              = zeros(1,size(x,2));
                        eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
                        eoni             = 1:size(x,2);
                        eoni             = find(eoni - eoi);
                    end
                end
                
                % ---------------------------
                % now deal with interactions
                % ---------------------------
                
                if nb_factors == 2 && nb_continuous == 0 % the quick way with only one interaction
                    HI(frame)                = T(frame,frame) - sum(H(:,frame)) - E;
                    df_interactions(frame)   = prod(df_conditions(frame));
                    F_interactions(frame)    = (HI(frame)./df_interactions(frame)) ./ (E/dfe(frame))';
                    pval_interactions(frame) = 1 - fcdf(F_interactions(frame), df_interactions(frame), dfe(frame));
                    
                else % run through each interaction
                    
                    % part of X unchanged
                    Main_effects = X(:,dummy_columns);
                    Cov_and_Mean = [X(:,covariate_columns) ones(size(Y,1),1)];
                    
                    % loop through interactions
                    % substituting and/or incrementing parts of X
                    Istart     = size(Main_effects,2)+1; % where we start interaction in X
                    Ilowbound  = size(Main_effects,2)+1;
                    for f=1:length(nb_interactions)
                        I              = X(:,Istart:(Istart+nb_interactions(f)-1));
                        if length(interaction{f}) == 2 % 1st oder interaction is main + I
                            x          = [Main_effects I Cov_and_Mean];
                        else % higher oder inteaction includes lower levels
                            Isize      = sum(nb_interactions(1:find(Ifactors == Ifactors(f),1) - 1));
                            Ihighbound = size(Main_effects,2)+Isize;
                            x          = [Main_effects X(:,Ilowbound:Ihighbound) I Cov_and_Mean];
                        end
                        eoibound = size(x,2) - size(I,2) - size(Cov_and_Mean,2);
                        
                        % run same model as above
                        wx                     = x.*repmat(W(:,frame),1,size(x,2));
                        betas                  = pinv(wx)*(Y(:,frame).*W(:,frame));
                        R                      = eye(size(Y,1)) - (wx*pinv(wx));
                        eoi                    = zeros(1,size(x,2));
                        eoi(eoibound+1:(eoibound+nb_interactions(f))) = eoibound+1:(eoibound+nb_interactions(f));
                        eoni                   = 1:size(x,2);
                        eoni                   = find(eoni - eoi);
                        C                      = eye(size(x,2));
                        C(:,eoni)              = 0;
                        C0                     = eye(size(x,2)) - C*pinv(C);
                        X0                     = wx*C0;
                        R0                     = eye(size(Y,1)) - (X0*pinv(X0));
                        M                      = R0 - R;
                        HI(f,frame)              = betas'*x'*M*x*betas;
                        df_interactions(f,frame) = prod(df_conditions(interaction{f},frame));
                        F_interactions(f,frame)  = (HI(f,frame)./df_interactions(f,frame)) ./ (E/dfe(frame))';
                        pval_interactions(f,:)   = 1 - fcdf(F_interactions(f,frame), df_interactions(f,frame), dfe(frame));
                        Istart                   = Istart+nb_interactions(f);
                    end
                end
            end
                    
            % -----------------------------------
            %% compute F for continuous variables
            % -----------------------------------
            
            if nb_continuous > 1
                
                N_conditions = sum(nb_conditions) + sum(nb_interactions);
                for n = 1:nb_continuous
                    C                                = zeros(size(X,2));
                    C(N_conditions+n,N_conditions+n) = 1; % pick up one regressor at a time
                    C0                               = eye(size(X,2)) - C*pinv(C);
                    X0                               = WX*C0; % all but rehressor of interest
                    R0                               = eye(size(Y,1)) - (X0*pinv(X0));
                    M                                = R0 - R; % hat matrix for regressor of interest
                    H                                = Betas(:,frame)'*X'*M*X*Betas(:,frame);
                    df_continuous(n,frame)           = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C) if OLS;
                    F_continuous(n,frame)            = (H./(df_continuous(n,frame))) ./ (E/dfe(frame));
                    pval_continuous(n,frame)         = 1 - fcdf(F_continuous(n,frame), 1, dfe(frame)); % dfe same as size(Y,1)-rank(X) if OLS
                end
            end
            
        end
        
        model.R2_univariate     = Rsquare';
        model.F                 = F_Rsquare';
        model.p                 = p_Rsquare';
        model.betas             = Betas;
        model.df                = mean(dof,2)'; % not the true dof but a summary
        model.W                 = W;
        model.betas_se          = betas_se;
        
        % record all in model
        % transpose to get the same format as OLS and WLS
        if nb_conditions > 0
            if nb_factors == 1
                model.conditions.F    = F_conditions';
                model.conditions.p    = pval_conditions';
                model.conditions.df   = [mean(df_conditions,2) ; repmat(mean(dfe,2),nb_factors,1)];
            else
                model.conditions.F    = F_conditions;
                model.conditions.p    = pval_conditions;
                model.conditions.df   = [mean(df_conditions,2) repmat(mean(dfe,2),nb_factors,1)];
            end
            
            if nb_interactions > 0
                model.interactions.F  = F_interactions;
                model.interactions.p  = pval_interactions;
                model.interactions.df = [mean(df_interactions,2) repmat(mean(dfe,2),numel(nb_interactions),1)];
            end
        end
        
        if nb_continuous ~=0
            if nb_factors == 0 && nb_continuous == 1
                model.continuous.F   = F_Rsquare';
                model.continuous.p   = p_Rsquare';
                model.continuous.df  = mean(dof,2)';
            else
                model.continuous.F   = F_continuous';
                model.continuous.p   = pval_continuous';
                model.continuous.df  = [mean(df_continuous,2) repmat(mean(dfe,2),nb_continuous,1)];
            end
        end
end

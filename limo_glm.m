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
%
% INPUTS:
%   Y               = 2D matrix of EEG data with format trials x frames
%                     (frame can be the concatenation for freq*time)
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
% - For WLS, since the weights are derived for the entire trial, using a
% dimention reduction the df are adjusted for the subspace spanned by the
% parameters see < WLS paper>
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
elseif nargin >= 7
    Y               = varargin{1};
    X               = varargin{2};
    nb_conditions   = varargin{3};
    nb_interactions = varargin{4};
    nb_continuous   = varargin{5};
    method          = varargin{6};
    Analysis        = varargin{7};
    if nargin == 9
        n_freqs     = varargin{8};
        n_times     = varargin{9};
    end
else
    error('varargin error')
end

if isempty(nb_conditions);   nb_conditions   = 0; end
if isempty(nb_interactions); nb_interactions = 0; end
if isempty(nb_continuous);   nb_continuous   = 0; end

nb_factors = numel(nb_conditions);
if nb_factors == 1 && nb_conditions == 0
    nb_factors = 0;
end

% -----------
%% Data check
% -----------
if ~isreal(Y)
    Y = abs(Y).^2;
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
    [Betas,W,rf] = limo_WLS(X,Y);
    WX           = X.*repmat(W,1,size(X,2));
    
elseif strcmp(method,'WLS-TF')
    % unpack the data
    [n_freq_times, N] = size(Y');
    if n_freq_times ~= n_freqs*n_times
        error('dimensions disagreement to reshape freq*time')
    else
        reshaped = nan(n_freqs, N, n_times);
        for tr = 1:N
            eft_3d = NaN(n_freqs,n_times);
            for tm = 1:n_times
                this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 2D tf
                eft_3d(:,tm) =Y(tr,this_freq_start_index:(this_freq_start_index+n_freqs-1))';
            end
            reshaped(:,tr,:) = eft_3d;
        end
        Y = reshaped;
        clear reshaped;
    end
    
    % get estimates per frequency band
    index1 = 1;
    Betas  = NaN(size(X,2),n_freqs,n_times);
    W      = NaN(size(X,1),n_freqs);
    WX     = cell(1,n_freqs);
    rf     = NaN(1,n_freqs);
    for f=1:n_freqs
        [Betas(:,f,:),W(:,f),rf(f)] = limo_WLS(X,squeeze(Y(f,:,:)));
        WX{f} = X .* repmat(W(:,f),1,size(X,2)); 
        index1=index1+1;
    end
    
elseif strcmp(method,'IRLS')
    [Betas,W] = limo_IRLS(X,Y);
    warning off
    % WX = X.*W per frame =  switch method
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
        R   = eye(size(Y,1)) - WX*pinv(WX);                                      % Projection onto E
        E   = Y'*R*Y;                                                            % SS Error
        
        % degrees of freedom
        % -------------------
        df = rank(WX)-1;
        if strcmp(method,'OLS')
            dfe = size(Y,1)-rank(WX);
        else
            % Satterthwaite approximation minus the number of dimensions removed by pcout to get W
            HM  = WX*pinv(WX); % Hat matrix, projection onto X
            dfe = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM)) - (rf-1);   
        end
        
        % model R^2
        % -----------
        if all(nb_conditions==0) && nb_continuous == 0       % just the mean
            Yhat           = X*Betas;
            H              = (Yhat-repmat(mean(Yhat),size(Y,1),1))'*(Y-repmat(mean(Yhat),size(Y,1),1));
        else
            C              = eye(size(X,2));
            C(:,size(X,2)) = 0;                              % all columns but the constant
            C0             = eye(size(X,2)) - C*pinv(C);     % only the constant
            X0             = WX*C0;                          % Reduced model design matrix
            R0             = eye(size(Y,1)) - (X0*pinv(X0)); % Projection onto error
            M              = R0 - R;                         % Projection matrix onto WXc
            H              = (Betas'*X'*M*X*Betas);          % SS Effects (X'*M*X is weighted)
        end
        Rsquare            = diag(H)./diag(T);               % Variance explained
        F_Rsquare          = (diag(H)./df) ./ (diag(E)/dfe);
        p_Rsquare          = 1 - fcdf(F_Rsquare, df, dfe);
        
        % update the model structure
        % ----------------------------
        
        model.W                 = W;
        model.betas             = Betas;
        model.betas_se          = Betas;
        for t=1:size(Y,2)
            model.betas_se(:,t) = sqrt(diag((E(t,t)/dfe)*pinv(WX'*WX)));
                                  % same as sqrt(E(t,t)/dfe)./ sqrt(sum(sum((WX-mean(WX)).^2)));
        end
        model.R2_univariate     = Rsquare;
        model.F                 = F_Rsquare;
        model.df                = [df dfe];
        model.p                 = p_Rsquare;
        
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
            model.conditions.p                   = pval_conditions;           
            
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
            model.conditions.p      = pval_conditions;
            
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
            model.conditions.p       = pval_conditions;
            
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
                Cov_and_Mean = [X(:,covariate_columns) ones(size(X,1),1)];
                
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
            model.interactions.p       = pval_interactions;
        end
        
        % -----------------------------------
        %% compute F for continuous variables
        % -----------------------------------
        
        if nb_continuous ~=0
            
            if nb_factors == 0 && nb_continuous == 1 % simple regression
                model.continuous.F  = F_Rsquare;
                model.continuous.df = [1 (size(Y,1)-rank(X))];
                model.continuous.p  = p_Rsquare;
                
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
                    X0                               = WX*C0; % all but regressor of interest
                    R0                               = eye(size(Y,1)) - (X0*pinv(X0));
                    M                                = R0 - R; % hat matrix for regressor of interest
                    H                                = Betas'*X'*M*X*Betas;
                    df_continuous(n)                 = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C) if OLS;
                    F_continuous(n,:)                = (diag(H)./(df_continuous(n))) ./ (diag(E)/dfe);
                    pval_continuous(n,:)             = 1 - fcdf(F_continuous(n,:), 1, dfe); % dfe same as size(Y,1)-rank(X) if OLS
                end
                
                model.continuous.F                   = F_continuous';
                model.continuous.df                  = [1 dfe];
                model.continuous.p                   = pval_continuous';
            end
        end
        
        % ---------------------------------------------------------------------
    case 'WLS-TF'
        % ---------------------------------------------------------------------
        
        model.W             = W;
        model.betas         = Betas;
        model.betas_se      = Betas;
        model.dfe           = NaN(n_freqs,2);
        model.conditions.df = squeeze(NaN(nb_factors,n_freqs,2));
        
        % iterate per freqency band
        % ---------------------------
        for freq=n_freqs:-1:1
            
            T   = (squeeze(Y(freq,:,:))-repmat(mean(squeeze(Y(freq,:,:))),size(Y,2),1))'*(squeeze(Y(freq,:,:))-repmat(mean(squeeze(Y(freq,:,:))),size(Y,2),1));  % SS Total (the data)
            R   = eye(size(Y,2)) - WX{freq}*pinv(WX{freq});
            E   = squeeze(Y(freq,:,:))'*R*squeeze(Y(freq,:,:));
            
            % degrees of freedom
            % -------------------
            df  = rank(WX{freq})-1;
            HM  = WX{freq}*pinv(WX{freq}); % Hat matrix, projection onto X
            dfe = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM)) - (rf(freq)-1);

            % model R^2
            % -----------
            C              = eye(size(X,2));
            C(:,size(X,2)) = 0;
            C0             = eye(size(X,2)) - C*pinv(C);
            X0             = WX{freq}*C0;
            R0             = eye(size(Y,2)) - (X0*pinv(X0));
            M              = R0 - R;
            H              = (squeeze(Betas(:,freq,:))'*X'*M*X*squeeze(Betas(:,freq,:)));
            Rsquare        = diag(H)./diag(T);
            F_Rsquare      = (diag(H)./df) ./ (diag(E)/dfe);
            p_Rsquare      = 1 - fcdf(F_Rsquare, df, dfe);
           
            % update the model structure
            % ----------------------------
            
            for t=1:size(Y,3)
                model.betas_se(:,freq,t) = sqrt(diag((E(t,t)/dfe)*pinv(WX{freq}'*WX{freq})));
            end
            model.R2_univariate(freq,:)  = Rsquare;
            model.F(freq,:)              = F_Rsquare;
            model.df(freq,:)             = [df dfe];
            model.p(freq,:)              = p_Rsquare;
            
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
                    X0                               = WX{freq}*C0; % here the reduced model includes the covariates
                    R0                               = eye(size(Y,2)) - (X0*pinv(X0));
                    M                                = R0 - R; % hat matrix for all categorical regressors (1 factor)
                    H                                = (squeeze(Betas(:,freq,:))'*X'*M*X*squeeze(Betas(:,freq,:)));
                    df_conditions                    = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C)-1 if OLS; same as tr(M)?
                    F_conditions                     = (diag(H)/df_conditions) ./ (diag(E)/dfe);
                    pval_conditions                  = 1 - fcdf(F_conditions(:), df_conditions, dfe);
                end
                
                model.conditions.F(freq,:)              = F_conditions';
                model.conditions.df(freq,:)             = [df_conditions ; dfe]';
                model.conditions.p(freq,:)              = pval_conditions;
                
                % ------------------------------------------------
            elseif nb_factors > 1  && isempty(nb_interactions) % N-ways ANOVA without interactions
                % ------------------------------------------------
                
                % --------------------------------------
                % compute F and p values of each factor
                % --------------------------------------
                
                df_conditions            = zeros(1,length(nb_conditions));
                F_conditions             = zeros(length(nb_conditions),size(Y,3));
                pval_conditions          = zeros(length(nb_conditions),size(Y,3));
                
                % define the effect of interest (eoi)
                eoi                      = zeros(1,size(X,2));
                eoi(1:nb_conditions(1))  = 1:nb_conditions(1);
                eoni                     = 1:size(X,2);
                eoni                     = find(eoni - eoi);
                
                for f = 1:length(nb_conditions)
                    C                    = eye(size(X,2));
                    C(:,eoni)            = 0; % set all but factor of interest to 0
                    C0                   = eye(size(X,2)) - C*pinv(C);
                    X0                   = WX{freq}*C0; % the reduced model include all but the factor f
                    R0                   = eye(size(Y,2)) - (X0*pinv(X0));
                    M                    = R0 - R; % hat matrix for factor f
                    H                    = (squeeze(Betas(:,freq,:))'*X'*M*X*squeeze(Betas(:,freq,:)));
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
                
                model.conditions.F(:,freq,:)      = F_conditions';
                model.conditions.df(:,freq,:)     = [df_conditions ; repmat(dfe,1,numel(df_conditions))];
                model.conditions.p(:,freq,:)      = pval_conditions;
            
                % ------------------------------------------------
            elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways ANOVA with interactions
                % ------------------------------------------------
                
                % ---------------------------------------------------
                % start by ANOVA without interaction for main effects
                % ---------------------------------------------------
                
                H                 = NaN(length(nb_conditions),size(Y,3));
                df_conditions     = NaN(1,length(nb_conditions));
                F_conditions      = NaN(length(nb_conditions),size(Y,3));
                pval_conditions   = NaN(length(nb_conditions),size(Y,3));
                HI                = NaN(length(nb_interactions),size(Y,3));
                df_interactions   = NaN(1,length(nb_interactions));
                F_interactions    = NaN(length(nb_interactions),size(Y,3));
                pval_interactions = NaN(length(nb_interactions),size(Y,3));
                
                % covariates
                covariate_columns = (sum(nb_conditions)+sum(nb_interactions)+1):(size(X,2)-1);
                
                % main effects
                dummy_columns = 1:sum(nb_conditions);
                
                % re-define X for main effects
                x = [X(:,dummy_columns) X(:,covariate_columns) ones(size(X,1),1)];
                
                % run same model as above with re-defined model x and
                % using the weights from the full model
                wx                       = x.*repmat(W(:,freq),1,size(x,2));
                betas                    = pinv(wx)*(Y.*repmat(W(:,freq),1,size(Y,2)));
                R                        = eye(size(Y,2)) - wx*pinv(wx);
                eoi                      = zeros(1,size(x,2));
                eoi(1:nb_conditions(1))  = 1:nb_conditions(1);
                eoni                     = 1:size(x,2);
                eoni                     = find(eoni - eoi);
                
                for f = 1:length(nb_conditions)
                    C                    = eye(size(x,2));
                    C(:,eoni)            = 0;
                    C0                   = eye(size(x,2)) - C*pinv(C);
                    X0                   = wx*C0;
                    R0                   = eye(size(Y,2)) - (X0*pinv(X0));
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
                
                model.conditions.F(:,freq,:)      = F_conditions';
                model.conditions.df(:,freq,:)     = [df_conditions ; repmat(dfe,1,numel(df_conditions))];
                model.conditions.p(:,freq,:)      = pval_conditions;

                % ---------------------------
                % now deal with interactions
                % ---------------------------
                
                if nb_factors == 2 && nb_continuous == 0 % the quick way with only one interaction
                    HI                 = diag(T)' - H(1,:) - H(2,:) - diag(E)';
                    df_interactions    = prod(df_conditions);
                    F_interactions     = (HI./df_interactions) ./ (diag(E)/dfe)';
                    
                else % run through each interaction
                    
                    % part of X unchanged
                    Main_effects = X(:,dummy_columns);
                    Cov_and_Mean = [X(:,covariate_columns) ones(size(X,1),1)];
                    
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
                        wx                     = x.*repmat(W(:,freq),1,size(x,2));
                        betas                  = pinv(wx)*(Y.*repmat(W(:,freq),1,size(Y,2)));
                        R                      = eye(size(Y,2)) - (wx*pinv(wx));
                        eoi                    = zeros(1,size(x,2));
                        eoi(eoibound+1:(eoibound+nb_interactions(f))) = eoibound+1:(eoibound+nb_interactions(f));
                        eoni                   = 1:size(x,2);
                        eoni                   = find(eoni - eoi);
                        C                      = eye(size(x,2));
                        C(:,eoni)              = 0;
                        C0                     = eye(size(x,2)) - C*pinv(C);
                        X0                     = wx*C0;
                        R0                     = eye(size(Y,2)) - (X0*pinv(X0));
                        M                      = R0 - R;
                        HI(f,:)                = diag((betas'*x'*M*x*betas))';
                        df_interactions(f)     = prod(df_conditions(interaction{f}));
                        F_interactions(f,:)    = (HI(f,:)./df_interactions(f)) ./ (diag(E)/dfe)';
                        pval_interactions(f,:) = 1 - fcdf(F_interactions(f,:), df_interactions(f), dfe);
                        Istart                 = Istart+nb_interactions(f);
                    end
                end
                
                model.interactions.F(:,freq,:)  = F_interactions;
                model.interactions.df(:,freq,:) = [df_interactions ; repmat(dfe,1,numel(df_interactions))]';
                model.interactions.p(:,freq,:)  = pval_interactions;
            end
            
            % -----------------------------------
            %% compute F for continuous variables
            % -----------------------------------
            
            if nb_continuous ~=0
                
                if nb_factors == 0 && nb_continuous == 1 % simple regression
                    model.continuous.F(:,freq)  = F_Rsquare;
                    model.continuous.df(:,freq) = [1 (size(Y,2)-rank(X))];
                    
                else % ANCOVA type of designs
                    
                    % pre-allocate space
                    df_continuous   = zeros(nb_continuous,size(Y,3));
                    F_continuous    = zeros(nb_continuous,size(Y,3));
                    
                    % compute
                    N_conditions = sum(nb_conditions) + sum(nb_interactions);
                    for n = 1:nb_continuous
                        C                                = zeros(size(X,2));
                        C(N_conditions+n,N_conditions+n) = 1; % pick up one regressor at a time
                        C0                               = eye(size(X,2)) - C*pinv(C);
                        X0                               = WX{freq}*C0; % all but regressor of interest
                        R0                               = eye(size(Y,2)) - (X0*pinv(X0));
                        M                                = R0 - R; % hat matrix for regressor of interest
                        H                                = squeeze(Betas(:,freq,:))'*X'*M*X*squeeze(Betas(:,freq,:));
                        df_continuous(n)                 = trace(M'*M)^2/trace((M'*M)*(M'*M)); % same as rank(C) if OLS;
                        F_continuous(n,:)                = (diag(H)./(df_continuous(n))) ./ (diag(E)/dfe);
                        pval_continuous(n,:)             = 1 - fcdf(F_continuous(n,:), 1, dfe); % dfe same as size(Y,1)-rank(X) if OLS
                   end
                    
                    model.continuous.F(:,freq,:)         = F_continuous';
                    model.continuous.df(:,freq,:)        = [1 dfe];
                    model.continuous.p(:,freq,:)         = pval_continuous';
                end
            end
        end
        
        % reshape to output the same size as input
        % except betas - keeping the parameters for modeling
        model.R2_univariate = reshape(model.R2_univariate, [n_freqs*n_times,1]);
        model.F             = reshape(model.F, [n_freqs*n_times,1]);
        model.p             = reshape(model.p, [n_freqs*n_times,1]);
        
        if isfield(model,'conditions')
            if nb_factors == 1
                model.conditions.F  = reshape(model.conditions.F, [n_freqs*n_times,1]);
                model.conditions.p  = reshape(model.conditions.p, [n_freqs*n_times,1]);
            else
                for f = length(nb_conditions):-1:1
                    tmpF(f,:)  = reshape(model.conditions.F(f,:,:), [n_freqs*n_times,1]);
                    tmpp(f,:)  = reshape(model.conditions.p(f,:,:), [n_freqs*n_times,1]);
                end
                model.conditions.F  = tmpF;
                model.conditions.p  = tmpp;
                clear tmpF tmpp
            end
        end
        
        if isfield(model,'interactions')
            if nb_factors == 2 && nb_continuous == 0
                model.interactions.F  = reshape(model.conditions.F, [n_freqs*n_times,1]);
                model.interactions.p  = reshape(model.conditions.p, [n_freqs*n_times,1]);
            else
                for f = length(nb_interactions)
                    tmpF(f,:)  = reshape(model.interactions.F(f,:,:), [n_freqs*n_times,1]);
                    tmpp(f,:)  = reshape(model.interactions.p(f,:,:), [n_freqs*n_times,1]);
                end
                model.conditions.F  = tmpF;
                model.conditions.p  = tmpp;
                clear tmpF tmpp
            end
        end
 
        if isfield(model,'continuous')
            if nb_factors == 0 && nb_continuous == 1
                model.continuous.F  = reshape(model.continuous.F, [n_freqs*n_times,1]);
                model.continuous.p  = reshape(model.continuous.p, [n_freqs*n_times,1]);
            else
                for f = 1:nb_continuous
                    tmpF(:,f)  = reshape(model.continuous.F(:,:,f), [n_freqs*n_times,1]);
                    tmpp(:,f)  = reshape(model.continuous.p(:,:,f), [n_freqs*n_times,1]);
                end
                model.continuous.F  = tmpF;
                model.continuous.p  = tmpp;
                clear tmpF tmpp
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
        
        if nb_continuous ~=0
            F_continuous    = NaN(nb_continuous,size(Y,2));
            pval_continuous = NaN(nb_continuous,size(Y,2));
            df_continuous   = NaN(nb_continuous,size(Y,2));
        end
        
        % start computing
        T   = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        for frame = 1:size(Y,2)
            % model stats
            % -------------------------------------------------------------
            WX                     = X.*repmat(W(:,frame),1,size(X,2));
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
            
            if nb_continuous ~=0
                
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
        model.df                = dof'; 
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

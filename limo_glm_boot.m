function model = limo_glm_boot(varargin)

% Boostrapped version for limo_glm
% Importantly it runs per electrodes but do N bootstraps to obtain
% the distributon of F (and associated p values) under H0.
% H0 is obtained by either by resampling from centered data (categorical designs)
% or sampling Y but leaving X intact, i.e. breaking the link between Y and X
%
% FORMAT:
% model = limo_glm_boot(Y,LIMO,boot_table)
% model = limo_glm_boot(Y,X,nb_conditions,nb_interactions,nb_continuous,method,analysis type,n_freqs,n_times,boot_table)
%
% INPUTS
%         Y = 2D matrix of EEG data with format trials x frames
%         LIMO is a structure that contains information below
%         X = 2 dimensional design matrix
%         nb_conditions = a vector indicating the number of conditions per factor
%         nb_interactions = a vector indicating number of columns per interactions
%         nb_continuous = number of covariates
%         method = 'OLS', 'WLS', 'IRLS' (bisquare)
%         analysis type =  'Time', 'Frequency' or 'Time-Frequency'
%         n_freqs is the nb of frequency bins
%         n_times is the nb of time bins
%         boot_table is an optional argument - this is the resampling table
%                    if one calls limo_glm_boot to loop throughout channels,
%                    this might a good idea to provide such table so that
%                    the same resampling applies to each channel
%
% See also
% LIMO_GLM_HANDLING, LIMO_GLM, LIMO_WLS, LIMO_IRLS
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% varagin
nboot = 599; %

if nargin == 2 || nargin == 3
    y               = varargin{1};
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
    else
        n_freqs = []; n_times =[];
    end
    
    if nargin == 2
        boot_table = randi(size(y,1),size(y,1),nboot);
    elseif nargin == 3
        boot_table = varargin{3};
        nboot = size(boot_table,2);
    end
    
elseif nargin == 10 || nargin == 11
    y               = varargin{1};
    X               = varargin{2};
    nb_conditions   = varargin{3};
    nb_interactions = varargin{4};
    nb_continuous   = varargin{5};
    method          = varargin{6};
    Analysis        = varargin{7};
    n_freqs         = varargin{8};
    n_times         = varargin{9};
    if nargin == 10
        boot_table = randi(size(y,1),size(y,1),nboot);
    elseif nargin == 11
        boot_table = varargin{11};
        nboot = size(boot_table,2);
    end
else
    error('varargin error in limo_glm_boot')
end

clear varargin
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

% ---------------
%% Make null data
% ---------------

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
        for cel=1:nb_conditions
            index = find(X(:,cel));
            centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),length(index),1);
        end
        
    elseif size(nb_conditions,2) > 1
        % create fake interaction to get groups
        [tmpX, interactions] = limo_make_interactions(X(:,1:(end-1)), nb_conditions);
        if length(interactions) == 1
            start_at = sum(nb_conditions);
        else
            start_at = sum(nb_conditions)+sum(interactions(1:end-1));
        end
        
        for cel=(start_at+1):(start_at+interactions(end))
            index = find(tmpX(:,cel));
            centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),size(y(index,:),1),1);
        end
    end
else 
    centered_y = y; % actually not centered
end

clear y
design = X;

% compute for each bootstrap
% ---------------------------
BETASB  = cell(1,nboot);
MODELR2 = cell(1,nboot);
MODELF  = cell(1,nboot);
MODELp  = cell(1,nboot);

if nb_factors ~= 0
    F_CONDVALUES = cell(1,nboot);
    p_CONDVALUES = cell(1,nboot);
end

if ~isempty(nb_interactions)
    F_INTERVALUES  = cell(1,nboot);
    p_INTERVALUES  = cell(1,nboot);
end

if nb_continuous ~=0
    F_CONTVALUES  = cell(1,nboot);
    p_CONTVALUES  = cell(1,nboot);
end

switch method
    
        % -----------------------------------------------------------------
    case {'OLS','WLS'}
        parfor B = 1:nboot
            
            % create data under H0
            Y = centered_y(boot_table(:,B),:); % resample Y
            if nb_continuous == 0
                % if just categorical variables, sample from the centered data and
                % the design simultaneously
                X = design(boot_table(:,B),:);     % resample X
            else
                % sample and break the link between Y and X (regression and AnCOVA designs)
                X = design;                        % stays the same
            end
            
            %% Compute model parameters
            % ------------------------------
            
            % compute Beta parameters
            if strcmp(method,'OLS')
                if strcmp(Analysis,'Time-Frequency')
                    W = ones(n_freqs,size(X,1));
                else
                    W = ones(size(Y,1),1);
                end
                WX = X;
                
                if nb_continuous ~=0 && nb_factors == 0
                    Betas = X\Y; % numerically more stable than pinv
                else
                    Betas = pinv(X)*Y;
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
                    eft_3d = nan(n_freqs,n_times);
                    for tm = 1:n_times
                        this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 2D tf
                        eft_3d(:,tm) =Y(tr,this_freq_start_index:(this_freq_start_index+n_freqs-1))';
                    end
                    reshaped(:,:,tr) = eft_3d;
                end
                
                % get estimates per freq band
                index1 = 1;
                Betas  = NaN(size(X,2),n_freqs*n_times);
                W      = NaN(n_freqs,size(X,1));
                for f=1:n_freqs
                    [Betas(:,index1:6:(n_freqs*n_times)),W(f,:)] = limo_WLS(X,squeeze(reshaped(f,:,:))');
                    index1=index1+1;
                end
                WX = X .* repmat(W,1,size(X,2));
            end
            
            % Betas bootstap
            BETASB{B} = Betas';
            
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
            
            
            % ----------------------------
            %% update the model structure
            % ----------------------------
            MODELR2{B} = Rsquare;
            MODELF{B}  = F_Rsquare;
            MODELp{B}  = p_Rsquare;
            
            %% Compute effects
            % ------------------
            
            % -------------------------
            if nb_factors == 1   %  1-way ANOVA
                % -------------------------
                
                % compute F for categorical variables
                % -----------------------------------
                if nb_conditions ~= 0 && nb_continuous == 0
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
                    F_conditions                     = (diag(H)/df) ./ (diag(E)/dfe);
                    pval_conditions                  = 1 - fcdf(F_conditions(:), df_conditions, dfe);
                end
                
                F_CONDVALUES{B}  = F_conditions;
                p_CONDVALUES{B}  = pval_conditions;
                
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
                
                F_CONDVALUES{B}  = F_conditions;
                p_CONDVALUES{B}  = pval_conditions;
                
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
                
                % re-define X
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
                
                F_CONDVALUES{B}  = F_conditions;
                p_CONDVALUES{B}  = pval_conditions;
                
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
                
                F_INTERVALUES{B}  = F_interactions;
                p_INTERVALUES{B}  = pval_interactions;
                
            end
            
            % -----------------------------------
            %% compute F for continuous variables
            % -----------------------------------
            
            if nb_continuous ~=0
                
                if nb_factors == 0 && nb_continuous == 1 % simple regression
                    F_CONTVALUES{B}  = F_Rsquare;
                    p_CONTVALUES{B}  = p_Rsquare;
                    
                else  % ANCOVA
                    
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
                    
                    F_CONTVALUES{B}  = F_continuous';
                    p_CONTVALUES{B}  = pval_continuous';
                end
            end
        end
        
        %% ---------------------------------------------------------------------
    case 'IRLS'
        parfor B = 1:nboot
            
            % create data under H0
            Y = centered_y(boot_table(:,B),:); % resample Y
            if nb_continuous == 0
                % if just categorical variables, sample from the centered data and
                % the design simultaneously
                X = design(boot_table(:,B),:);     % resample X
            else
                % sample and break the link between Y and X (regression and AnCOVA designs)
                X = design;                        % stays the same
            end
            
            tmp        = limo_glm(Y, X, nb_conditions, nb_interactions, nb_continuous, method, Analysis, n_freqs, n_times);
            BETASB{B}  = tmp.betas;
            MODELR2{B} = tmp.R2_univariate;
            MODELF{B}  = tmp.F;
            MODELp{B}  = tmp.p;
            
            if nb_factors ~= 0
                F_CONDVALUES{B} = tmp.conditions.F;
                p_CONDVALUES{B} = tmp.conditions.p;
            end
            
            if ~isempty(nb_interactions)
                F_INTERVALUES{B}  = tmp.interactions.F;
                p_INTERVALUES{B}  = tmp.interactions.F;
            end
            
            if nb_continuous ~=0
                F_CONTVALUES{B}  = tmp.continuous.F;
                p_CONTVALUES{B}  = tmp.continuous.F;
            end
        end
end

model.R2_univariate  = MODELR2; clear MODELR2
model.F              = MODELF;  clear MODELF
model.p              = MODELp;  clear MODELp
model.betas          = BETASB;  clear BETASB
if nb_factors ~= 0
    model.conditions.F   = F_CONDVALUES; clear F_CONDVALUES
    model.conditions.p   = p_CONDVALUES; clear p_CONDVALUES
end
if ~isempty(nb_interactions)
    model.interactions.F = F_INTERVALUES; clear F_INTERVALUES
    model.interactions.p = p_INTERVALUES; clear p_INTERVALUES
end
if nb_continuous ~=0
    model.continuous.F   = F_CONTVALUES; clear F_CONTVALUES
    model.continuous.p   = p_CONTVALUES; clear p_CONTVALUES
end



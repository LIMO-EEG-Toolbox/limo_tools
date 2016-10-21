function model = limo_glm1(varargin)

% General Linear model for 1st level analysis of EEG data
% The model consider trials as independent observations
% i.e. this is similar as running a N-way ANOVA or ANCOVA
% Analyses are performed at one electrode only, but for all
% trials and all time points.
%
% FORMAT:
% model = limo_glm1(Y,LIMO)
% model = limo_glm1(Y,X,nb_conditions,nb_interactions,nb_continuous,method,analysis type,n_freqs,n_times)
%
% INPUTS:
%   Y             = 2D matrix of EEG data with format trials x frames
%   LIMO          = structure that contains the above information (except Y)
%   or input info that are in LIMO.mat
%   X             = 2 dimensional design matrix
%   nb_conditions = a vector indicating the number of conditions per factor
%   nb_interactions = a vector indicating number of columns per interactions
%   nb_continuous = number of covariates
%   method        = 'OLS', 'WLS', 'IRLS' (bisquare)
%   analysis type =  'Time', 'Frequency' or 'Time-Frequency'
%   n_freqs       = the nb of frequency bins
%   n_times       = the nb of time bins
%
% OUTPUTS:
%   model.R2_univariate = the R2 of the model
%   model.F             = the F value of the model
%   model.df            = the df associated to the model F
%   model.p             = the p value of the model
%   model.betas dim     = the beta parameters (dimensions nb of paramters x frames)
%   model.conditions    = categorical effects
%          --> F/p in rows are the factors, in columns time frames
%          --> df row 1 = df, row 2 = dfe, columns are factors
%   model.continuous = continuous effects
%          --> F/p in rows are the variables, in columns time frames
%          --> df column 1 = df, column2 2 = dfe (same for all covariates)
%
% NOTES:
%
% The parameters can be computed using 3 methods: ordinary least squares, 
% weighted least squares and iterative reweighted least squares. 
% Each effect is accounted for given the other effects; this means
% that one can have a different number of trials per conditions/factors 
% provided there is no interactions - for interactions this is not possible
% and no correction is provided - a design created by limo_design_matrix
% would have sampled trials to make sure the number of trials is identical
% across interaction terms. For time*frequency analyses, limo_eeg_tf send a
% vector of length freq*time that is unwrap after running limo_glm1 -
% however, weights cannot be computed that way and thus one calls LIMO.Analysis
% to unwrap and then rewarp to get the right weights and betas -- this
% implies that to the run 'by hand' limo_glm1 for time freqency with WLS,
% you need to run per freq bin.
%
% References
% Christensen, R. 2002. Plane answers to complex questions. 3rd Ed. Springer-Verlag
% Friston et al. 2007. Statitical Parametric Mapping. Academic Press
% Yandell, B.S. 1997. Practical Data Analysis For Designed Experiments. Chapman & Hall
%
% See also
% LIMO_DESIGN_MATRIX, LIMO_WLS, LIMO_IRLS, LIMO_EEG
%
% Cyril Pernet v1 01-01-2011
% Cyril Pernet v2 01-11-2011
% Cyril Pernet v3 07-07-2015 (methods and analysis type)
% ---------------------------------------------------------
%  Copyright (C) LIMO Team 2015

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

% total sum of squares, projection matrix for errors, residuals 
% --------------------------------------------------------------
T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
R     = eye(size(Y,1)) - (X*pinv(X));                                      % Projection on E
E     = (Y'*R*Y);                                                          % SS Error

% compute Beta parameters and weights
if strcmp(method,'OLS')
    
    if strcmp(Analysis,'Time-Frequency')
        W = ones(n_freqs,size(X,1));
    else
        W = ones(size(Y,1),1);
    end
    
    if nb_continuous ~=0 && nb_factors == 0
        Betas = X\Y; % numerically more stable than pinv
    else
        Betas = pinv(X)*Y;
    end
    
elseif strcmp(method,'WLS')
    [Betas,W] = limo_WLS(X,Y);
    
elseif strcmp(method,'WLS-TF')
    % unpack the data
    [n_freq_times, N] = size(Y');
    if n_freq_times ~= n_freqs*n_times
        error('dimensions disagreement to reshape freq*time')
    else
        reshaped = nan(n_freqs, n_times, N);
    end

    for tr = 1:N
        for tm = 1:n_times
            this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 2D tf 
                eft_3d(:,tm) =Y(tr,this_freq_start_index:(this_freq_start_index+n_freqs-1))';
        end
        reshaped(:,:,tr) = eft_3d;
    end
    
     % get estimates per freq band
    Betas = NaN(size(X,2),n_freqs*n_times);
    W = NaN(n_freqs,size(X,1));
    index1 = 1; 
    for f=1:n_freqs
        [Betas(:,index1:6:(n_freqs*n_times)),W(f,:)] = limo_WLS(X,squeeze(reshaped(f,:,:))');
        index1=index1+1; 
    end
    clear reshaped
    
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
H  = (Betas'*X'*M*X*Betas);  % SS Effect
Rsquare   = diag(H)./diag(T); % Variances explained
F_Rsquare = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
p_Rsquare = 1 - fcdf(F_Rsquare, (rank(X)-1), (size(Y,1)-rank(X)));

%% Compute effects
% ------------------

% -------------------------
if nb_factors == 1   %  1-way ANOVA
% -------------------------

    
    % compute F for categorical variables
    % -----------------------------------
    if nb_conditions ~= 0 && nb_continuous == 0
        df_conditions   = rank(C)-1;
        F_conditions    = F_Rsquare;
        pval_conditions = p_Rsquare;
        
    elseif nb_conditions ~= 0 && nb_continuous ~= 0
        C = eye(size(X,2));
        C(:,(nb_conditions+1):size(X,2)) = 0;
        C0 = eye(size(X,2)) - C*pinv(C);
        X0 = X*C0; % Here the reduced model includes the covariates
        R0 = eye(size(Y,1)) - (X0*pinv(X0));
        M  = R0 - R;
        H  = (Betas'*X'*M*X*Betas);
        df_conditions = rank(C)-1;
        F_conditions    = (diag(H)/(rank(C)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
        pval_conditions = 1 - fcdf(F_conditions(:), df_conditions, (size(Y,1)-rank(X)));
    end
    
    model.conditions.F  = F_conditions;
    model.conditions.p  = pval_conditions;
    model.conditions.df = [df_conditions ; (size(Y,1)-rank(X))];
    
% ------------------------------------------------
elseif nb_factors > 1  && isempty(nb_interactions) % N-ways ANOVA without interactions
% ------------------------------------------------
    
    % --------------------------------------
    % compute F and p values of each factor
    % --------------------------------------

    df_conditions = zeros(1,length(nb_conditions));
    F_conditions = zeros(length(nb_conditions),size(Y,2));
    pval_conditions = zeros(length(nb_conditions),size(Y,2));

    eoi = zeros(1,size(X,2));
    eoi(1:nb_conditions(1)) = 1:nb_conditions(1);
    eoni = [1:size(X,2)];
    eoni = find(eoni - eoi);
    
    for f = 1:length(nb_conditions)
        C = eye(size(X,2));
        C(:,eoni) = 0;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0;
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;
        H    = (Betas'*X'*M*X*Betas);
        df_conditions(f) = rank(C)-1;
        F_conditions(f,:)    = (diag(H)/df_conditions(f)) ./ (diag(E)/(size(Y,1)-rank(X)));
        pval_conditions(f,:) = 1 - fcdf(F_conditions(f,:), df_conditions(f), (size(Y,1)-rank(X)));
        
        % update factors
        if f<length(nb_conditions)
            update = max(find(eoi));
            eoi = zeros(1,size(X,2));
            eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
            eoni = [1:size(X,2)];
            eoni = find(eoni - eoi);
        end
    end
    model.conditions.F  = F_conditions;
    model.conditions.p  = pval_conditions;
    model.conditions.df = [df_conditions ; repmat((size(Y,1)-rank(X)),1,numel(df_conditions))]';

       
% ------------------------------------------------
elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways ANOVA with interactions
% ------------------------------------------------
    
    % ---------------------------------------------------
    % start by ANOVA without interaction for main effects
    % ---------------------------------------------------

    df_conditions = zeros(1,length(nb_conditions));
    F_conditions = zeros(length(nb_conditions),size(Y,2));
    pval_conditions = zeros(length(nb_conditions),size(Y,2));
    
    % covariates
    covariate_columns = [(sum(nb_conditions)+sum(nb_interactions)+1):(size(X,2)-1)];
    
    % main effects
    dummy_columns = 1:sum(nb_conditions);
    
    % re-define X
    x = [X(:,dummy_columns) X(:,covariate_columns) ones(size(X,1),1)];
     
    % run same model as above
    R        = eye(size(Y,1)) - (x*pinv(x));   
    if strcmp(method,'IRLS')
        betas = pinv(W*x)*(W*Y);
    else
        betas = pinv(repmat(W,1,size(x,2)).*x)*(repmat(W,1,size(Y,2)).*Y);
    end
             
    eoi = zeros(1,size(x,2));
    eoi(1:nb_conditions(1)) = 1:nb_conditions(1);
    eoni = [1:size(x,2)];
    eoni = find(eoni - eoi);
    
    for f = 1:length(nb_conditions)
        C = eye(size(x,2));
        C(:,eoni) = 0;
        C0   = eye(size(x,2)) - C*pinv(C);
        X0   = x*C0;
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;
        H(f,:) = diag((betas'*x'*M*x*betas));
        df_conditions(f) = rank(C)-1;
        F_conditions(f,:)    = (H(f,:)./df_conditions(f)) ./ (diag(E)./(size(Y,1)-rank(X)))';
        pval_conditions(f,:) = 1 - fcdf(F_conditions(f,:), df_conditions(f), (size(Y,1)-rank(X)));
        
        % update factors
        if f<length(nb_conditions)
            update = max(find(eoi));
            eoi = zeros(1,size(x,2));
            eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
            eoni = [1:size(x,2)];
            eoni = find(eoni - eoi);
        end
    end
    model.conditions.F  = F_conditions;
    model.conditions.p  = pval_conditions;
    model.conditions.df = [df_conditions ; repmat((size(Y,1)-rank(X)),1,numel(df_conditions))]';

    % ---------------------------
    % now deal with interactions
    % ---------------------------
    
    if nb_factors == 2 && nb_continuous == 0 % the quick way with only one interaction
        HI = diag(T)' - H(1,:) - H(2,:) - diag(E)';
        df_interactions = prod(df_conditions);
        F_interactions  = (HI./df_interactions) ./ (diag(E)/(size(Y,1)-rank(X)))';
        pval_interactions  = 1 - fcdf(F_interactions, df_interactions, (size(Y,1)-rank(X)));
    
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
        % run substituting and/or incrementing parts of X
        for f = 1:length(nb_interactions)
            
            % re-define X with interactions
            test = size(interaction{f},2);
            if test == 2
               x = [Main_effects I{f} Cov_and_Mean];
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
                betas = pinv(repmat(W,1,size(x,2)).*x)*(repmat(W,1,size(Y,2)).*Y);
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
        end
        
        % get appropriate df and F/p values
        df_interactions = zeros(1,length(nb_interactions));
        F_interactions = zeros(length(nb_interactions),size(Y,2));
        pval_interactions = zeros(length(nb_interactions),size(Y,2));

        for f = 1:length(nb_interactions)
            dfs = df_conditions(interaction{f});
            df_interactions(f) = prod(dfs);
            F_interactions(f,:) = (HI(f,:)./df_interactions(f)) ./ (diag(E)/(size(Y,1)-rank(X)))';
            pval_interactions(f,:) = 1 - fcdf(F_interactions(f,:), df_interactions(f), (size(Y,1)-rank(X)));
        end
    end
    model.interactions.F  = F_interactions;
    model.interactions.p  = pval_interactions;
    model.interactions.df = [df_interactions ; repmat((size(Y,1)-rank(X)),1,numel(df_interactions))]';
end


% -----------------------------------
%% compute F for continuous variables
% -----------------------------------

if nb_continuous ~=0
    
    if nb_factors == 0 && nb_continuous == 1 % simple regression
        model.continuous.F  = F_Rsquare;
        model.continuous.p  = p_Rsquare;
        model.continuous.df = [1 (size(Y,1)-rank(X))];
    
    else % ANCOVA type of deisgns
        
        % pre-allocate space
        F_continuous = zeros(nb_continuous,size(Y,2));
        pval_continuous = zeros(nb_continuous,size(Y,2));
        
        % compute
        N_conditions = sum(nb_conditions) + sum(nb_interactions);
        for n = 1:nb_continuous
            C    = zeros(size(X,2));
            C(N_conditions+n,N_conditions+n) = 1;
            C0   = eye(size(X,2)) - C*pinv(C);
            X0   = X*C0;
            R0   = eye(size(Y,1)) - (X0*pinv(X0));
            M    = R0 - R;
            H    = Betas'*X'*M*X*Betas;
            F_continuous(n,:) = (diag(H)./(rank(C))) ./ (diag(E)/(size(Y,1)-rank(X)));
            pval_continuous(n,:) = 1 - fcdf(F_continuous(n,:), 1, (size(Y,1)-rank(X)));
        end
        model.continuous.F  = F_continuous';
        model.continuous.p  = pval_continuous';
        model.continuous.df = [1 (size(Y,1)-rank(X))];
    end
end

% ----------------------------
%% update the model structure
% ----------------------------

model.R2_univariate   = Rsquare;
model.F               = F_Rsquare;
model.p               = p_Rsquare;
model.betas           = Betas;
model.df              = [rank(X)-1 (size(Y,1)-rank(X))];
model.W               = W';


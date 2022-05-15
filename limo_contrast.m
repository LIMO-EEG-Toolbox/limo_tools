function result = limo_contrast(varargin)

% limo_contrast computes contrasts (i.e. differences between regressors)
% using outputs from main statitistical tests (limo_glm.m,
% limo_hotelling.m). The function uses the parameters computed, reads the
% design matrix and compute the contrast and statistical test associated to
% it.
%
% FORMATS:
% result = limo_contrast(Y, Betas, LIMO, contrast type, analysis type ,contrast)
%          --> applies to 1st level, 2nd level regressions and 2nd level ANOVA/ANCOVA
% result = limo_contrast(Yr,LIMO, analysis type ,contrast);
%          --> applies to 2nd level repeated measures ANOVA
%
% INPUT:
% Y              = the data as matrix or file name
% Betas          = the model parameters as matrix or file name 
%                  Betas.mat for analysis_type = 1 or H0_Betas.mat for analysis type = 2
% LIMO           = the LIMO structure or LIMO file name
% contrast type  = 0 or 'T' for T test, 1 or 'F' for F test
% analysis type  = 1 Contrast for 1st level analyses and 2nd level regression/ANOVA/ANCOVA
%                  2 for 1st level analyses and 2nd level bootstrap regression/ANOVA/ANCOVA
% analysis type  = 3 for 2nd level repeated measures ANOVA
%                  4 for 2nd level bootrapped repeated measures ANOVA
% contrast       = optional a contrast to test ; if not specified the
%                  contrast should be in LIMO.contrast and last one is
%                  evaluated
%
% OUTPUT
% con/ess maps saved on disk
% these files are of dimension [nb of channels, time/freq, C*Beta/se/df/t/p]
%
% *****************************************************
% See also limo_contrast_checking, limo_glm, limo_results, limo_contrast_manager
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2021


%% nargin stuff
if nargin == 5 || nargin == 6
    type = varargin{5};
elseif nargin == 4
    type = varargin{3};
end

%% default
result = [];
warning on

%% Analyses

if type == 1 || type == 2
    % ---------------------------------------------------------------------
    %         1st level / 2nd level regressions and ANOVA/ANCOVA
    % ---------------------------------------------------------------------
    Y     = varargin{1};
    if ischar(Y)
        Y = load(varargin{1});
        Y = Y.(cell2mat(fieldnames(Y)));
    end
    
    Betas     = varargin{2};
    if ischar(Betas)
        Betas = load(varargin{2});
        Betas = Betas.(cell2mat(fieldnames(Betas)));
        if type == 2 && size(Betas,numel(size(Betas))) < 101
            warning('input Betas file is not a H0 one, checking for a H0 boostraps file')
            if exist(fullfile(fileparts(varargin{2}),['H0' filesep 'H0_Betas.mat']),'file')
                Betas = load(fullfile(fileparts(varargin{2}),['H0' filesep 'H0_Betas.mat']));
                Betas = Betas.(cell2mat(fieldnames(Betas)));    
                if size(Betas,numel(size(Betas))) < 101
                    error('loading H0_Betas.mat but this seems to have less than 101 bootstraps?')
                end
            else
                warning('contrast with boostrap aborded no H0 file found')
                return
            end
        end
    end
    
    LIMO     = varargin{3};
    if ischar(LIMO)
        LIMO = load(varargin{3});
        LIMO = LIMO.LIMO;
    end
    
    if contains(LIMO.design.name,'Repeated','IgnoreCase',true)
        error('2nd level Repeated measure Analysis detected ; switch analysis type');
    end
    
    X       = LIMO.design.X;
    nb_beta = size(LIMO.design.X,2);
    if isfield(LIMO.model,'model_df')
        dfe = LIMO.model.model_df(:,2:end);
    else
        dfe = size(Y,1)-rank(X); %% happens for 2nd level N-way ANOVA or ANCOVA
    end
    if length(dfe) == 1
        dfe = repmat(dfe,1,size(Y,1)); % dfe per channel
    end
    
    if nargin == 6 && type == 1 % <------ nargin = 6, user input a contrast
        if isfield(LIMO,'contrast')
            contrast_nb = size(LIMO.contrast,2)+1;
        else
            contrast_nb = 1;
        end
        out = limo_contrast_checking(LIMO.dir,LIMO.design.X,varargin{6}); % add zeros if needed
        
        if limo_contrast_checking(out,LIMO.design.X) % if contrast is valid
            if varargin{4} == 1 || strcmpi(varargin{4},'T')
                if size(out,1) == 1
                    LIMO.contrast{contrast_nb}.V = 'T';
                else
                    warning('the specificed contrast is on multiples rows, using F constrast')
                    LIMO.contrast{contrast_nb}.V = 'F';
                end
            else
                if size(out,1) == 1
                    warning('the specificed contrast is on one row, using T constrast')
                    LIMO.contrast{contrast_nb}.V = 'T';
                else
                    LIMO.contrast{contrast_nb}.V = 'F';
                end
            end
            LIMO.contrast{contrast_nb}.C = out;
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        else
           error('invalid contrast ass input') 
        end
    elseif nargin == 6 && type == 2 % <---- find the index of the contrast to bootstrap
        if ~isfield(LIMO,'contrast')
            warning('analysis type = 2; no constrast to boostrap found like the one as input')
            return
        else
            allC  = cellfun(@(x) x.C,LIMO.contrast,'UniformOutput',false);
            out   = limo_contrast_checking(LIMO.dir,LIMO.design.X,varargin{6});
            contrast_nb = max(cellfun(@(x) all(x==out), allC));
        end
        
        if contrast_nb == 0
            warning('analysis type = 2; no constrast to boostrap found like the one as input')
            return
        end
        
    elseif nargin == 5 %<--- nothing specifed = bootstrap the last one
        contrast_nb = size(LIMO.contrast,2);
        go = limo_contrast_checking(LIMO.contrast{contrast_nb}.C,LIMO.design.X);
        if go ==0
            error('the analysis of the %g contrast in LIMO.contrast failed, invalid contrast',contrast_nb)
        end
    end
    C       = LIMO.contrast{contrast_nb}.C;
    Method  = LIMO.design.type_of_analysis;
    
    % legacy naming convention
    if strcmpi(varargin{4},'T')
        Test = 0;
    elseif strcmpi(varargin{4},'F')
        Test = 1;
    else
        Test = varargin{4};
    end
    
elseif type == 3 || type == 4
    % ---------------------------------------------------------------------
    %                  2nd level repeated measures ANOVA
    % ---------------------------------------------------------------------
    Yr         = varargin{1};
    if ischar(Yr)
        Yr = load(varargin{1});
        Yr = Yr.(cell2mat(fieldnames(Yr)));
    end
    
    LIMO       = varargin{2};
    if ischar(LIMO)
        LIMO = load(LIMO);
        LIMO = LIMO.LIMO;
    end
    if LIMO.Level == 1
        error('1st level Analysis detected - limo_contrast wrong case ; switch analysis type');
    elseif ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
        error('2nd level Analysis but not a Repeated measure Analysis ; switch analysis type');
    end
    gp_values  = LIMO.design.nb_conditions;
    
    if nargin == 4 && type == 3
        if isfield(LIMO,'contrast')
            LIMO.contrast{end+1}.V = 'F';
            LIMO.contrast{end}.C = varargin{4};
        else
            LIMO.contrast{1}.V = 'F';
            LIMO.contrast{1}.C = varargin{4};
        end
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    end
    
    if ~isfield(LIMO,'contrast')
        error('no contrast found to evaluate')
    else
        if nargin == 4 && type == 4
            allC  = cellfun(@(x) x.C,LIMO.contrast,'UniformOutput',false);
            index = max(cellfun(@(x) all(x==varargin{4}), allC));
            if index == 0
               warning('analysis type = 2; no constrast to boostrap found like the one as input')
               return
            end
        else
            index      = size(LIMO.contrast,2);
        end
        C          = LIMO.contrast{index}.C;
        Test       = 2; % always a F-test
    end
end
clear varargin


%% start the analysis
switch type
    
    case{1}
        % -----------------------------------------------------------------
        % Contrast for 1st level analyses and 2nd level regression/ANOVA/ANCOVA
        % -----------------------------------------------------------------
        
        % get residuals
        Res   = load([LIMO.dir filesep 'Res.mat']);   
        Res   = Res.(cell2mat(fieldnames(Res)));
        
        % string time-frequency for OLS and IRLS
        if strcmp(LIMO.Analysis ,'Time-Frequency') && strcmpi(LIMO.design.method,'OLS') || ...
                strcmp(LIMO.Analysis ,'Time-Frequency') && strcmpi(LIMO.design.method,'IRLS')
            Y     = limo_tf_4d_reshape(Y);
            Betas = limo_tf_4d_reshape(Betas);
            Res   = limo_tf_4d_reshape(Res);
        end
        
        if strcmp(Method,'Mass-univariate')
            if strcmp(LIMO.Analysis ,'Time-Frequency') && strcmpi(LIMO.design.method,'WLS')
                % create con or ess file
                if Test == 0
                    con      = NaN(size(Y,1),size(Y,2),size(Y,3),5); % dim 3 = C*Beta/se/df/t/p
                    filename = sprintf('con_%g.mat',size(LIMO.contrast,2));
                else
                    ess      = NaN(size(Y,1),size(Y,2),size(Y,3),size(C,1)+4); % dim 3 = C*Beta/se/df/F/p
                    filename = sprintf('ess_%g.mat',size(LIMO.contrast,2));
                end
                
                warning off;
                array = find(~isnan(Y(:,1,1))); % skip empty channels
                for e = 1:length(array)
                    channel = array(e); 
                    if strcmp(LIMO.Type,'Channels')
                        fprintf('applying contrast on channel %g/%g \n',e,size(array,1));
                    else
                        fprintf('applying contrast on component %g/%g \n',e,size(array,1));
                    end
                    
                    % contrasts
                    % -----------
                    for freq = 1:size(Y,2)
                        if Test == 0 % T contrast
                            
                            % Update con file [mean value, se, df, t, p]
                            var                    = (squeeze(Res(channel,freq,:,:))*squeeze(Res(channel,freq,:,:))') / dfe(channel,freq);
                            con(channel,freq,:,1)  = C*squeeze(Betas(channel,freq,:,:))';
                            con(channel,freq,:,3)  = dfe(channel,freq);
                            WX                     = X.*repmat(squeeze(LIMO.design.weights(channel,freq,:)),1,size(X,2));
                            con(channel,freq,:,2)  = sqrt(diag(var)'.*(C*pinv(WX'*WX)*C')); % var is weighted already
                            con(channel,freq,:,4)  = (C*squeeze(Betas(channel,freq,:,:))') ./ sqrt(diag(var)'.*(C*pinv(WX'*WX)*C'));
                            con(channel,freq,:,5)  = (1-tcdf(squeeze(abs(con(channel,freq,:,4))), dfe(channel,freq))).*2;
                        else % F contrast
                            % Update ess file [mean values, se, df, F, p]
                            E = diag(squeeze(Res(channel,freq,:,:))*squeeze(Res(channel,freq,:,:))');
                            ess(channel,freq,:,1:size(C,1)) = (C*squeeze(Betas(channel,freq,:,:))')' ;
                            ess(channel,freq,:,end-3)       = E/dfe(channel,freq);
                            if rank(diag(C)) == 1
                                df = 1;
                            else
                                df = rank(diag(C)) - 1;
                            end
                            ess(channel,freq,:,end-2) = df;
                            
                            c  = zeros(length(C));
                            C0 = eye(size(c,1)) - diag(C)*pinv(diag(C));
                            WX = X.*repmat(squeeze(LIMO.design.weights(channel,freq,:)),1,size(X,2));
                            R  = eye(size(Y,4)) - (WX*pinv(WX));
                            X0 = X*C0;
                            R0 = eye(size(Y,4)) - (X0*pinv(X0));
                            M  = R0 - R;
                            H  = (squeeze(Betas(channel,freq,:,:))*X'*M*X*squeeze(Betas(channel,freq,:,:))');
                            ess(channel,freq,:,end-1) = (diag(H)/df)./(E/dfe(channel,freq));  % F value
                            ess(channel,freq,:,end)   = 1 - fcdf(ess(channel,freq,:,end-1), df, dfe(channel,freq)); % p value
                        end
                    end
                end
                warning on;
                
            else % all other data/methods
                
                % create con or ess file
                if Test == 0
                    con      = NaN(size(Y,1),size(Y,2),5); % dim 3 = C*Beta/se/df/t/p
                    filename = sprintf('con_%g.mat',size(LIMO.contrast,2));
                else
                    ess      = NaN(size(Y,1),size(Y,2),size(C,1)+4); % dim 3 = C*Beta/se/df/F/p
                    filename = sprintf('ess_%g.mat',size(LIMO.contrast,2));
                end
                
                % update con/ess file
                warning off;
                array = find(~isnan(Y(:,1,1))); % skip empty channels
                for e = 1:length(array)
                    channel = array(e);
                    if strcmp(LIMO.Type,'Channels')
                        fprintf('applying contrast on channel %g/%g \n',e,size(array,1));
                    else
                        fprintf('applying contrast on component %g/%g \n',e,size(array,1));
                    end
                    
                    % contrasts
                    % -----------
                    if Test == 0 % T contrast
                        
                        % Update con file [mean value, se, df, t, p]
                        if strcmpi(LIMO.design.method,'OLS') || strcmpi(LIMO.design.method,'WLS')
                            var                      = (squeeze(Res(channel,:,:))*squeeze(Res(channel,:,:))') / dfe(channel); % sum of (xi-mean)^2 since res are xi-mean take res^2, dived by dfe ie n-dimensions of the mean
                            con(channel,:,1)         = C*squeeze(Betas(channel,:,:))'; % how do we scale axes of WX
                            WX                       = X.*repmat(LIMO.design.weights(channel,:)',1,size(X,2)); 
                            con(channel,:,2)         = sqrt(diag(var)'.*(C*pinv(WX'*WX)*C')); % var = avg distance to model projected into the contrast space
                            con(channel,:,3)         = dfe(channel);
                            con(channel,:,4)         = (C*squeeze(Betas(channel,:,:))') ./ sqrt(diag(var)'.*(C*pinv(WX'*WX)*C'));
                            con(channel,:,5)         = (1-tcdf(squeeze(abs(con(channel,:,4))), dfe(channel))).*2; 
                        elseif strcmpi(LIMO.design.method,'IRLS')
                            var                      = diag((squeeze(Res(channel,:,:))*squeeze(Res(channel,:,:))') ./ dfe(channel,:));
                            for frame = 1:size(Betas,2)
                                con(channel,frame,1) = C*squeeze(Betas(channel,frame,:));
                                WX                   = X.*repmat(squeeze(LIMO.design.weights(channel,frame,:)),1,size(X,2));
                                con(channel,frame,2) = sqrt(var(frame).*(C*pinv(WX'*WX)*C'));
                                con(channel,frame,3) = dfe(channel,frame);
                                con(channel,frame,4) = (C*squeeze(Betas(channel,frame,:))) ./ sqrt(var(frame).*(C*pinv(WX'*WX)*C'));
                                con(channel,frame,5) = (1-tcdf(squeeze(abs(con(channel,frame,4))), dfe(channel,frame))).*2;
                            end
                        end
                    else % F contrast
                        % Update ess file [mean values, se, df, F, p]
                        E = diag(squeeze(Res(channel,:,:))*squeeze(Res(channel,:,:))');
                        ess(channel,:,1:size(C,1)) = (C*squeeze(Betas(channel,:,:))')' ;
                        if rank(diag(C)) == 1
                            df = 1;
                        else
                            df = rank(diag(C)) - 1;
                        end
                        ess(channel,:,end-2) = df;
                        
                        c  = zeros(length(C));
                        C0 = eye(size(c,1)) - diag(C)*pinv(diag(C));
                        if strcmpi(LIMO.design.method,'OLS') || strcmpi(LIMO.design.method,'WLS')
                            if isfield(LIMO.design,'weights')
                                WX = X.*repmat(squeeze(LIMO.design.weights(channel,:)'),1,size(X,2));
                            else
                                WX = X;
                            end
                            R  = eye(size(Y,3)) - (WX*pinv(WX));
                            X0 = X*C0;
                            R0 = eye(size(Y,3)) - (X0*pinv(X0));
                            M  = R0 - R;
                            H  = (squeeze(Betas(channel,:,:))*X'*M*X*squeeze(Betas(channel,:,:))');
                            ess(channel,:,end-3) = E/dfe(channel);
                            ess(channel,:,end-1) = (diag(H)/df)./(E/dfe(channel));  % F value
                            ess(channel,:,end)   = 1 - fcdf(ess(channel,:,end-1), df, dfe(channel)); % p value
                        else % IRLS
                            for frame = 1:size(Betas,2)
                                WX = X.*repmat(squeeze(LIMO.design.weights(channel,frame,:)),1,size(X,2));
                                R  = eye(size(Y,3)) - (WX*pinv(WX));
                                X0 = X*C0;
                                R0 = eye(size(Y,3)) - (X0*pinv(X0));
                                M  = R0 - R;
                                H  = (squeeze(Betas(channel,frame,:))'*X'*M*X*squeeze(Betas(channel,frame,:)));
                                ess(channel,:,end-3) = E(frame)/dfe(channel);
                                ess(channel,frame,end-1) = (H/df)./(E(frame)/dfe(channel));  % F value
                                ess(channel,frame,end)   = 1 - fcdf(ess(channel,frame,end-1), df, dfe(channel)); % p value
                            end
                        end
                    end
                end
                warning on;
               
                % reshape Time-Frequency files
                if strcmp(LIMO.Analysis ,'Time-Frequency')
                    if Test == 0
                        con = limo_tf_4d_reshape(con);
                    else
                        ess = limo_tf_4d_reshape(ess);
                    end
                end
            end
            
            % save files
            if nargout == 1 && Test == 0
                result = con;
            elseif nargout == 1 && Test == 1
                result = ess;
            else
                if Test == 0
                    save(fullfile(LIMO.dir,filename),'con'); clear con
                else
                    save (fullfile(LIMO.dir,filename),'ess'); clear ess
                end
            end
            
            if LIMO.design.tfce == 1
                if ~exist(fullfile(LIMO.dir,['tfce' filesep 'tfce_' filename]),'file')
                    limo_tfce_handling(fullfile(LIMO.dir,filename));
                end
            end
            
        elseif strcmp(Method,'Multivariate')
            % ------------------------------
            
            con = NaN(size(Y,2),2); %  F /p values (always the same no matter RoY or Pillai)
            for time = 1:size(Y,2)
                fprintf('time frame %g \n',time);
                
                E = (squeeze(Y(:,time,:))*R*squeeze(Y(:,time,:))');
                c = zeros(length(C));
                for n=1:length(C)
                    c(n,n) = C(n);
                end
                
                try
                    C0 = eye(rank(X)+1) - c*pinv(c);
                catch ME
                    C0 = eye(rank(X)) - c*pinv(c);
                end
                X0 = X*C0;
                R0 = eye(size(Y,2)) - (X0*pinv(X0));
                M = R0 - R;
                H = (squeeze(Betas(:,time,:))*X'*M*X*squeeze(Betas(:,time,:))');
                
                multivariate.EV    = limo_decomp(E,H);
                multivariate.theta = max(multivariate.EV) / (1+max(multivariate.EV));
                multivariate.V     = sum(multivariate.EV ./ (1+multivariate.EV));
                multivariate.df    = size(Y,2);
                multivariate.dfe   = abs(size(Y,1) - (nb_beta-1) - (multivariate.df-1));
                multivariate.T_contrast    = sqrt((mean(dfe)*max(multivariate.EV))/multivariate.df);
                multivariate.pval_contrast = 1-fcdf(multivariate.T_contrast, multivariate.df, mean(abs(dfe)));
                % to do save into LIMO + con file
            end
        end
        
    case{2}
        % -----------------------------------------------------------------
        % bootstraps
        % -----------------------------------------------------------------
        nboot = LIMO.design.bootstrap;
        if nboot == 1
            nboot = 800;
        end
        
        if strcmp(LIMO.Analysis ,'Time-Frequency') && strcmpi(LIMO.design.method,'OLS') || ...
                strcmp(LIMO.Analysis ,'Time-Frequency') && strcmpi(LIMO.design.method,'IRLS')
            Y  = limo_tf_4d_reshape(Y);
        end
         
         % make data files
        % ----------------
        if Test == 0
            H0_con   = NaN(size(Y,1),size(Y,2),2,nboot); % dim 3 = t/p
            filename = sprintf('H0_con_%g.mat',size(LIMO.contrast,2));
        else
            H0_ess   = NaN(size(Y,1),size(Y,2),2,nboot); % dim 3 = F/p
            filename = sprintf('H0_ess_%g.mat',size(LIMO.contrast,2));
        end
        
        
        % prepare data for bootstrap as in limo_glm_boot
        % ---------------------------------------------
        for channel=size(Y,1):-1:1
            centered_data(channel,:,:) = limo_glm_null(squeeze(Y(channel,:,:))',...
                X,LIMO.design.nb_conditions,LIMO.design.nb_interactions)';
        end
        
        % start the analysis
        % -------------------
        boot_table = load(fullfile(LIMO.dir,['H0' filesep 'boot_table.mat']));
        boot_table = boot_table.(cell2mat(fieldnames(boot_table)));
        array = find(~isnan(Y(:,1,1))); % skip empty channels
        design = X;
        
        if strcmp(Method,'Mass-univariate')
            % ---------------------------------
            warning off;
            for e = 1:length(array)
                channel = array(e);
                for B = 1:nboot
                    if ~iscell(boot_table)
                        resampling_index = boot_table(:,B); % 1st level boot_table all the same
                    else
                        resampling_index = boot_table{channel}(:,B);
                    end
                    
                    % create data under H0
                    Y = squeeze(centered_data(channel,:,resampling_index))';
                    
                    if strcmp(LIMO.design.method,'OLS') || strcmp(LIMO.design.method,'WLS')
                        fprintf('compute bootstrap channel %g ... \n',channel)
                        trials_to_keep = ~isnan(Y(:,1));
                        Y              = Y(trials_to_keep,:);
                        X              = design(trials_to_keep,:); % do not resample X
                        W              = LIMO.design.weights(channel,trials_to_keep)';
                        if any(isnan(Y(:,1))) && ...
                                LIMO.design.nb_continuous ~= 0 && ...
                                LIMO.design.zscore == 1 % rezscore the covariates
                            N = LIMO.design.nb_conditions + LIMO.design.nb_interactions;
                            if any(mean(X(:,N+1:end-1),1) > 10e-15)
                                X(:,N+1:end-1) = zscore(X(:,N+1:end-1));
                            end
                        end
                        
                        % compute Projection onto the error
                        WX = X .* repmat(W,1,size(X,2));
                        R  = eye(size(Y,1)) - WX*pinv(WX);
                        
                        % T contrast
                        % -----------
                        if Test == 0
                            var   = ((R*Y)'*(R*Y)) / dfe(channel); % error of H0 data
                            H0_con(channel,:,1,B) = (C*squeeze(Betas(channel,:,:,B))') ./ sqrt(diag(var)'.*(C*pinv(WX'*WX)*C')); % T value
                            H0_con(channel,:,2,B) = 1-tcdf(squeeze(H0_con(channel,:,1,B)), dfe(channel)); % p value
                            
                            % F contrast
                            % ----------
                        else
                            E = (Y'*R*Y);
                            c = zeros(length(C));
                            for n=1:length(C)
                                c(n,n) = C(n);
                            end
                            C0 = eye(size(c,2)) - c*pinv(c);
                            X0 = WX*C0;
                            R0 = eye(size(Y,1)) - (X0*pinv(X0));
                            M = R0 - R;
                            H = (squeeze(Betas(channel,:,:,B))*X'*M*X*squeeze(Betas(channel,:,:,B))');
                            df = rank(c) - 1;
                            if df == 0
                                df = 1;
                            end
                            H0_ess(channel,:,1,B) = (diag(H)/df)./(diag(E)/dfe(channel));  % F value
                            H0_ess(channel,:,2,B) = 1 - fcdf(H0_ess(channel,:,1,B), rank(c)-1, dfe(channel));   % p value
                        end
                        
                        
                    else % -------- IRLS ------------
                        for frame = 1:size(Y,2)
                           fprintf('compute bootstrap channel %g frame %g/%g ... \n',channel, frame,size(Y,2))
                            X = design; % do not resample X
                            W = squeeze(LIMO.design.weights(channel,frame,~isnan(Y(:,1))));
                            if any(isnan(Y(:,1)))
                                Y = Y(~isnan(Y(:,1)),:);
                                X = X(~isnan(Y(:,1)),:);
                                W = W(~isnan(Y(:,1)),:);
                                if LIMO.design.nb_continuous ~= 0 && LIMO.design.zscore == 1 % rezscore the covariates
                                    N = LIMO.design.nb_conditions + LIMO.design.nb_interactions;
                                    if N==0
                                        if sum(mean(X(:,1:end-1),1)) > 10e-15
                                            X(:,1:end-1) = zscore(X(:,1:end-1));
                                        end
                                    else
                                        if sum(mean(X(:,N+1:end-1),1)) > 10e-15
                                            X(:,N+1:end-1) = zscore(X(:,N+1:end-1));
                                        end
                                    end
                                end
                            end
                            
                            WX  = X .* repmat(W,1,size(X,2));
                            HM  = WX*pinv(WX);
                            R   = eye(size(Y,1)) - HM;
                            dfe = trace((eye(size(HM))-HM)'*(eye(size(HM))-HM));
                            
                            % T contrast
                            % -----------
                            if Test == 0
                                var   = ((R*Y(:,frame))'*(R*Y(:,frame))) / dfe; % error of H0 data
                                H0_con(channel,frame,1,B) = (C*squeeze(Betas(channel,frame,:,B))) ./ sqrt(var.*(C*pinv(X'*X)*C')); % T value
                                H0_con(channel,frame,2,B) = 1-tcdf(squeeze(H0_con(channel,frame,2,B)), dfe); % p value
                                
                                % F contrast
                                % ----------
                            else
                                E = (Y(:,frame)'*R*Y(:,frame));
                                c = zeros(length(C));
                                for n=1:length(C)
                                    c(n,n) = C(n);
                                end
                                C0 = eye(size(c,2)) - c*pinv(c);
                                X0 = WX*C0;
                                R0 = eye(size(Y,1)) - (X0*pinv(X0));
                                M = R0 - R;
                                H = squeeze(Betas(channel,frame,:,B))'*X'*M*X*squeeze(Betas(channel,frame,:,B));
                                df = rank(c) - 1;
                                if df == 0
                                    df = 1;
                                end
                                H0_ess(channel,frame,1,B) = (diag(H)/df)./(diag(E)/dfe);  % F value
                                H0_ess(channel,frame,2,B) = 1 - fcdf(H0_ess(channel,frame,end-1,B), rank(c)-1, dfe);   % p value
                            end
                        end
                    end
                end
            end
            warning on;
            
            if Test == 0
                save (fullfile(LIMO.dir,['H0' filesep filename]), 'H0_con'); clear H0_con; 
            else
                save (fullfile(LIMO.dir,['H0' filesep filename]), 'H0_ess'); clear H0_ess; 
            end
        end
        
        if LIMO.design.tfce == 1
            if ~exist(fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' filename]),'file')
                limo_tfce_handling(fullfile(LIMO.dir,filename(4:end)),'checkfile','no');
            end
        end
            
        % ----------------------------------------
        if strcmp(Method,'Multivariate')
            % ----------------------------------------
            
            warning off;
            for e = 1:size(Y,1)
                channel = array(e);
                fprintf('compute bootstrap channel %g ... \n',channel)
                for B = 1:nboot
                    % create data under H0
                    if LIMO.design.nb_continuous == 0
                        % sample from the centered data in categorical designs
                        Y = centered_data(boot_table(:,B));
                        X = design(boot_table(:,B)); % resample X as Y
                    else
                        % sample and break the link between Y and (regression and AnCOVA designs)
                        Y = Y(boot_table(:,B));
                        if LIMO.design.zscore == 1 % rezscore the covariates
                            N = LIMO.design.nb_conditions + LIMO.design.nb_interactions;
                            if sum(mean(X(:,N+1:end-1),1)) ~= 0
                                X(:,N+1:end-1) = zscore(X(:,N+1:end-1));
                            end
                        end
                    end
                    
                    % compute Projection onto the error
                    R = eye(size(Y,1)) - (X*pinv(X));
                    
                    E = (Y'*R*Y);
                    c = zeros(length(C));
                    for n=1:length(C)
                        c(n,n) = C(n);
                    end
                    
                    try
                        C0 = eye(rank(X)+1) - c*pinv(c);
                    catch ME
                        C0 = eye(rank(X)) - c*pinv(c);
                    end
                    X0 = X*C0;
                    R0 = eye(size(Y,1)) - (X0*pinv(X0));
                    M = R0 - R;
                    H = (Betas'*X'*M*X*Betas);
                    
                    multivariate.EV    = limo_decomp(E,H);
                    multivariate.theta = max(multivariate.EV) / (1+max(multivariate.EV));
                    multivariate.V     = sum(multivariate.EV ./ (1+multivariate.EV));
                    multivariate.df    = size(Y,2);
                    multivariate.dfe   = abs(size(Y,1) - (nb_beta-1) - (multivariate.df-1));
                    multivariate.T_contrast    = sqrt((mean(dfe)*max(multivariate.EV))/multivariate.df);
                    multivariate.pval_contrast = 1-fcdf(multivariate.T_contrast, multivariate.df, mean(abs(dfe)));
                    result = multivariate;
                end
            end
            warning on;
       end
        
    case(3)
        % --------------------------------------------
        %              Repeated Measure ANOVA
        % ---------------------------------------------
        
        cd(LIMO.dir);
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tmp = NaN(size(Yr,1), size(Yr,2)*size(Yr,3),size(Yr,4),size(Yr,5));
            for measure = 1:size(Yr,5)
                if size(Yr,1) == 1
                    tmp(:,:,:,measure) = limo_tf_4d_reshape(Yr(1,:,:,:,measure));
                else
                    tmp(:,:,:,measure) = limo_tf_4d_reshape(squeeze(Yr(:,:,:,:,measure)));
                end
            end
            clear Yr; Yr = tmp; clear tmp
        end
        
        % [mean value, se, df, F, p])
        if gp_values == 1
            ess = zeros(size(Yr,1),size(Yr,2),5);
            array = find(nansum(squeeze((Yr(:,1,:,1))),2));    
            for c = 1:length(array)
                channel = array(c);
                fprintf('channel %g \n',channel);
                % Inputs
                tmp = squeeze(Yr(channel,:,:,:));
                Y   = tmp(:,find(~isnan(tmp(1,:,1))),:);
                gp  = LIMO.data.Cat(find(~isnan(tmp(1,:,1))),:);
                % mean, se, df
                n = size(Y,2);
                if strcmpi(LIMO.design.method,'Mean')
                    for time=1:size(Y,1)
                        ess(channel,time,1) = nanmean(C(1:size(Y,3))*squeeze(Y(time,:,:))',2);
                        ess(channel,time,2) = sqrt(C(1:size(Y,3))*cov(squeeze(Y(time,:,:)))*C(1:size(Y,3))');
                    end
                else
                    ess(channel,:,1) = limo_trimmed_mean([1 Y]);
                    for time=1:size(Y,1)
                        ess(channel,time,2) = sqrt(C(1:size(Y,3))*cov(squeeze(ess(channel,time,1)))*C(1:size(Y,3))');
                    end
                end
                df  = rank(C); dfe = n-df;
                ess(channel,:,3) = dfe;
                % F and p
                if strcmpi(LIMO.design.method,'Mean')
                    result = limo_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(Y,3)));
                else
                    result = limo_robust_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(Y,3)));
                end
                ess(channel,:,4) = result.F;
                ess(channel,:,5) = result.p;
            end
        else
            ess  = zeros(size(Yr,1),size(Yr,2),5); % dim rep measures, F,p
            ess2 = zeros(size(Yr,1),size(Yr,2),5); % dim gp*interaction F,p
            
            % design matrix for gp effects
            gp_vector = LIMO.data.Cat;
            gp_values = unique(gp_vector); 
            k         = length(gp_values); 
            X         = NaN(size(gp_vector,1),k+1);
            for g =1:k
                X(:,g) = gp_vector == gp_values(g); 
            end
            X(:,end) = 1; 
            
            % call rep anova
            for channel = 1:size(Yr,1)
                fprintf('channel %g \n',channel);
                % Inputs
                tmp = squeeze(Yr(channel,:,:,:));
                Y   = tmp(:,find(~isnan(tmp(1,:,1))),:);
                gp  = LIMO.data.Cat(find(~isnan(tmp(1,:,1))),:);
                XB  = X(find(~isnan(tmp(1,:,1))),:);
                % mean, se, df
                n = size(Y,2);
                if strcmpi(LIMO.design.method,'Mean')
                    g = 0; % < ----------- bellow the code is trimmed mean and winsorized variance
                           %               but with g = 0 this is regular mean and varance
                else
                    g = floor((20/100)*n);
                end
                
                for frame=1:size(Y,1)
                    [v,indices]          = sort(squeeze(Y(frame,:,:))); % sorted data
                    TD(frame,:,:)         = v((g+1):(n-g),:);           % trimmed data
                    ess(channel,frame,1)  = nanmean(C(1:size(TD,3))*squeeze(TD(frame,:,:))',2);
                    I                    = zeros(1,1,n); 
                    I(1,1,:)             = (C(1:size(TD,3))*squeeze(Y(frame,:,:))')'; % interaction
                    ess2(channel,frame,1) = limo_trimmed_mean(I);
                    v(1:g+1,:)           = repmat(v(g+1,:),g+1,1);
                    v(n-g:end,:)         = repmat(v(n-g,:),g+1,1);      % winsorized data
                    [~,reorder]          = sort(indices);
                    for j = 1:size(Y,3)
                        SD(:,j) = v(reorder(:,j),j); % restore the order of original data
                    end 
                    S(frame,:,:)          = cov(SD);  % winsorized covariance
                    ess(channel,frame,2)  = sqrt(C(1:size(TD,3))*squeeze(S(frame,:,:))*C(1:size(TD,3))');
                    ess2(channel,frame,2) = NaN;
                end
                df  = rank(C); dfe = n-df;
                ess(channel,:,3) = dfe;
                
                % F and p values
                if strcmpi(LIMO.design.method,'Mean')
                    result = limo_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(TD,3)),XB);
                else
                    result = limo_robust_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(TD,3)),XB);
                end
                ess(channel,:,4)  = result.repeated_measure.F;
                ess(channel,:,5)  = result.repeated_measure.p;
                ess2(channel,:,4) = result.interaction.F;
                ess2(channel,:,5) = result.interaction.p;
            end
        end
        
        filename = sprintf('ess_%g.mat',index);
        if strcmp(LIMO.Analysis,'Time-Frequency') 
            ess = limo_tf_4d_reshape(ess);
        end
        save(filename, 'ess', '-v7.3');
        
        if exist('ess2','var')
            if strcmp(LIMO.Analysis,'Time-Frequency')
                ess = limo_tf_4d_reshape(ess2);
            else
                ess = ess2;
            end
            filename2 = sprintf('ess_gp_interaction_%g.mat',index);
            save(filename2, 'ess', '-v7.3');
        end
        
        % tfce if needed
        if LIMO.design.tfce ~= 0
            if ~exist(fullfile(LIMO.dir,['tfce' filesep 'tfce_' filename]),'file')
                limo_tfce_handling(fullfile(LIMO.dir,filename));               
            end
            
            if exist('ess2','var') && ~exist(fullfile(LIMO.dir,['tfce' filesep 'tfce_' filename2]),'file')
                limo_tfce_handling(fullfile(LIMO.dir,filename2));               
            end
        end
        
    case(4)
        % --------------------------------------------
        %              bootstrap
        % ---------------------------------------------
        
        filename = fullfile(LIMO.dir,['H0' filesep 'H0_ess_' num2str(index) '.mat']);
        % prepare the boostrap with centering the data
        cd([LIMO.dir filesep 'H0']);
        if ~exist('centered_data.mat','file') || ~exist('boot_table.mat','file')
            error('H0 data and/or resampling table missing')
        else
            centered_data = load('centered_data'); centered_data = centered_data.(cell2mat(fieldnames(centered_data)));
            boot_table    = load('boot_table');    boot_table = boot_table.(cell2mat(fieldnames(boot_table)));
        end
        
        if gp_values == 1
            if strcmp(LIMO.Analysis,'Time-Frequency')
                H0_ess = NaN(size(Yr,1),size(Yr,2)*size(Yr,3),2,LIMO.design.bootstrap);
            else
                H0_ess = NaN(size(Yr,1),size(Yr,2),2,LIMO.design.bootstrap);
            end
            clear Yr
            
            %  compute
            array = find(nansum(squeeze((centered_data(:,1,:,1))),2));          
            fprintf('bootstrapping contrast ...\n');
            parfor b = 1:LIMO.design.bootstrap
                H0_ess_sub = NaN(size(centered_data,1),size(centered_data,2),2);
                for c = 1:length(array)
                    channel = array(c);
                    if c == 1
                        fprintf('parallel boot %g channel %g',b,channel);
                    elseif c==length(array)
                        fprintf(' %g\n',channel);
                    else
                        fprintf(' %g',channel);
                    end
                    % Inputs
                    tmp = squeeze(centered_data(channel,:,boot_table{channel}(:,b),:));
                    Y   = tmp(:,find(~isnan(tmp(1,:,1))),:); % resampling should not have NaN, JIC
                    gp  = LIMO.data.Cat(find(~isnan(tmp(1,:,1))));
                    % F and p
                    if strcmpi(LIMO.design.method,'Mean')
                        result = limo_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(Y,3)));
                    else
                        result = limo_robust_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(Y,3)));
                    end
                    H0_ess_sub(channel,:,1) = result.F;
                    H0_ess_sub(channel,:,2) = result.p;
                end
                H0_ess(:,:,:,b) = H0_ess_sub;
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                H0_ess = limo_tf_5d_reshape(H0_ess,LIMO);
            end
            save(filename, 'H0_ess', '-v7.3');

        else %% group*repeated measures
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                H0_ess = NaN(size(Yr,1),size(Yr,2)*size(Yr,3),2,LIMO.design.bootstrap);  
                H0_ess2 = NaN(size(Yr,1),size(Yr,2)*size(Yr,3),2,LIMO.design.bootstrap); 
            else
                H0_ess  = NaN(size(Yr,1),size(Yr,2),2,LIMO.design.bootstrap); 
                H0_ess2 = NaN(size(Yr,1),size(Yr,2),2,LIMO.design.bootstrap); 
            end
            clear Yr
            
            % design matrix for gp effects
            gp_vector = LIMO.data.Cat; 
            gp_values = unique(gp_vector); 
            k         = length(gp_values); 
            X         = NaN(size(gp_vector,1),k+1);
            for g =1:k
                X(:,g) = gp_vector == gp_values(g); 
            end 
            X(:,end) = 1; 
            
            % call rep anova
            array = find(nansum(squeeze((centered_data(:,1,:,1))),2));
            fprintf('bootstrapping contrast ...\n');
            parfor b = 1:LIMO.design.bootstrap
                H0_ess_sub  = NaN(size(centered_data,1), size(centered_data,2),2);
                H0_ess2_sub = NaN(size(centered_data,1), size(centered_data,2),2);
                for c = 1:length(array)
                    channel = array(c);
                    if c == 1
                        fprintf('parallel boot %g channel %g',b,channel);
                    elseif c==length(array)
                        fprintf(' %g\n',channel);
                    else
                        fprintf(' %g',channel);
                    end
                    % Inputs
                    tmp = squeeze(centered_data(channel,:,boot_table{channel}(:,b),:));
                    Y   = tmp(:,find(~isnan(tmp(1,:,1))),:);
                    gp  = LIMO.data.Cat(find(~isnan(tmp(1,:,1))),:); % adjust gp as Y
                    XB  = X(find(~isnan(tmp(1,:,1))),:); % adjust X as well
                    % F and p values
                    if strcmpi(LIMO.design.method,'Mean')
                        result = limo_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(Y,3)));
                    else
                        result = limo_robust_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(Y,3)));
                    end
                    H0_ess_sub(channel,:,1)  = result.repeated_measure.F;
                    H0_ess_sub(channel,:,2)  = result.repeated_measure.p;
                    H0_ess2_sub(channel,:,1) = result.interaction.F;
                    H0_ess2_sub(channel,:,2) = result.interaction.p;
                end
                H0_ess(:,:,:,b)  = H0_ess_sub;
                H0_ess2(:,:,:,b) = H0_ess2_sub;
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                H0_ess = limo_tf_5d_reshape(H0_ess,LIMO);
            end
            save(filename, 'H0_ess', '-v7.3');

            if exist('H0_ess2','var')
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    H0_ess = limo_tf_5d_reshape(H0_ess2,LIMO);
                else
                    H0_ess = H0_ess2;
                end
                filename2 = fullfile(LIMO.dir,['H0' filesep 'H0_ess_gp_interaction_' num2str(index) '.mat']);
                save(filename2, 'H0_ess', '-v7.3');
            end
        end
        cd(LIMO.dir); 
        
        % tfce if needed
        if LIMO.design.tfce ~= 0
            if ~exist(fullfile(LIMO.dir,['H0' filesep 'tfce_H0_ess_' num2str(index) '.mat']),'file')
                limo_tfce_handling(fullfile(LIMO.dir,['ess_' num2str(index) '.mat']),'checkfile','no');
            end
            
            if exist('ess2','var') && ~exist(fullfile(LIMO.dir,['H0' filesep 'tfce_H0_ess_gp_interaction_' num2str(index) '.mat']),'file')
                limo_tfce_handling(fullfile(LIMO.dir,['ess_gp_interaction_' num2str(index) '.mat']),'checkfile','no');
            end
        end
        disp('contrast bootstrap done')
end
save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');


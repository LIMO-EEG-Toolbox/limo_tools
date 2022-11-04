function limo_glm_handling(LIMO)

% data handling function for GLM
% this function calls limo_glm, limo_glm_boot to get the analysis done, and
% organize all the files around that - externalized from limo_eeg(4)
%
% FORMAT limo_glm_handling(LIMO)
%
% Cyril Pernet
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2020

cd(LIMO.dir);
warning on;

%% Compute GLM and save stats files

if strcmp(LIMO.design.status,'to do')
    Yr    = load('Yr');    Yr    = Yr.(cell2mat(fieldnames(Yr)));  
    Yhat  = load('Yhat');  Yhat  = Yhat.(cell2mat(fieldnames(Yhat)));
    Res   = load('Res');   Res   = Res.(cell2mat(fieldnames(Res)));
    R2    = load('R2');    R2    = R2.(cell2mat(fieldnames(R2)));
    Betas = load('Betas'); Betas = Betas.(cell2mat(fieldnames(Betas)));
  
    % check method and change parametyers accordingly
    % -----------------------------------------------
    if size(Yr,1) == 1 % in any cases, just one channel/component
        array = 1;
    else
        if LIMO.Level == 2        % second level we can have missing data because of 
            array = [1:size(Yr,1)]'; %#ok<NBRAK> % bad channels for some subjects - adjust X
        else % level 1 = skip empty channels
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                array = find(~isnan(Yr(:,1,1,1))); 
            else
                array = find(~isnan(Yr(:,1,1)));
            end
        end
    end
    
    if strcmpi(LIMO.design.method,'IRLS') % 1st or 2nd level
        N = size(Yr,numel(size(Yr)));
        if N < 50
            LIMO.design.method = 'OLS';
            warning('with %g observations detected, IRLS won''t converge, switching to OLS',N)
        end        
    end
    
    % check dimensions (3D vs 4D)
    % --------------------------------------
    if strcmpi(LIMO.Analysis,'Time-Frequency') 
         [~,n_freqs,n_times,~] = size(Yr); 
         Yr    = limo_tf_4d_reshape(Yr); 
         Yhat  = limo_tf_4d_reshape(Yhat); % reshape to 3D
         Res   = limo_tf_4d_reshape(Res);
         R2    = limo_tf_4d_reshape(R2);
         Betas = limo_tf_4d_reshape(Betas);
    end

    % ------------ prepare condition/covariates -------------------
    if LIMO.design.nb_conditions ~=0
        tmp_Condition_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_conditions),2);
    end
    
    if LIMO.design.nb_interactions ~=0
        tmp_Interaction_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_interactions),2);
    end
    
    if LIMO.design.nb_continuous ~=0
        tmp_Covariate_effect = NaN(size(Yr,1),size(Yr,2),LIMO.design.nb_continuous,2);
    end
    
    % ------------- prepare weight matrix  -------------------------------------
    if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            W = ones(size(Yr,1),n_freqs,size(Yr,4));
        else
            W = ones(size(Yr,1),size(Yr,3));
        end
    elseif strcmp(LIMO.design.method,'IRLS')
        W = ones(size(Yr));
    end
    
    % ------------ run limo_glm per channels ---------------------------
    update = 1;
    X      = LIMO.design.X;
    if isfield(LIMO,'model')
        LIMO = rmfield(LIMO,'model');
    end
    
    warning off;
    for e = 1:length(array)
        channel = array(e); 
        if LIMO.Level == 2
            fprintf('analyzing channel %g/%g \n',e,size(array,1));
            Y             = squeeze(Yr(channel,:,:));
            index         = find(~isnan(Y(1,:))); % which subjects to keep
            if isempty(index)
                index     = 1:size(Y,2);
            end
            Y             = Y(:,index);
            LIMO.design.X = X(index,:);
            model         = limo_glm(Y',LIMO); warning on;
        else % level 1 we should not have any NaNs because we use 'array'
            if strcmp(LIMO.Type,'Channels')
                fprintf('analyzing channel %g/%g \n',e,size(array,1));
            else
                fprintf('analyzing component %g/%g \n',e,size(array,1));
            end
            index = 1:size(Yr,3);
            model = limo_glm(squeeze(Yr(channel,:,:))',LIMO);
        end
        
        % update the LIMO.mat 
        if update == 1 && strcmpi(LIMO.design.method,'OLS')
            LIMO.model.model_df = model.df;
            if LIMO.design.nb_conditions ~=0
                LIMO.model.conditions_df  = model.conditions.df;
            end
            if LIMO.design.nb_interactions ~=0
                LIMO.model.interactions_df  = model.interactions.df;
            end
            if LIMO.design.nb_continuous ~=0
                LIMO.model.continuous_df  = model.continuous.df;
            end
            update = 0;
        
        elseif update == 1 && ~strcmpi(LIMO.design.method,'OLS') % each channel can have different weighting and thus different df 
            
            % store temporarily as cell everything
            LIMO.model.model_df{channel} = model.df;          
            if LIMO.design.nb_conditions ~=0
                LIMO.model.conditions_df{channel}  = squeeze(model.conditions.df);
            end
            if LIMO.design.nb_interactions ~=0
                LIMO.model.interactions_df{channel}  = squeeze(model.interactions.df);
            end
            if LIMO.design.nb_continuous ~=0
                LIMO.model.continuous_df{channel}  = squeeze(model.continuous.df);
            end

            % 1 cell per channel 
            if e == size(array,1)
                tmp = cell2mat(LIMO.model.model_df)'; % dim (elec*[df dfe]) x 1 or time
                if size(tmp,2) == size(Yr,1)*2        % when dim 2 is shorter matlab can switch dim around in limo_glm :-(
                   tmp = tmp'; 
                end
                df  = tmp(1:2:end,1); % a single value over time
                dfe = tmp(2:2:end,:); % could be different over time
                LIMO.model = rmfield(LIMO.model,'model_df');
                LIMO.model.model_df = [df dfe]; clear tmp
                
                if LIMO.design.nb_conditions ~=0
                    tmp = cell2mat(LIMO.model.conditions_df)'; % dim (elec*[df dfe]) * 1
                    if size(tmp,1) == size(Yr,1)*2 
                        df  = tmp(1:2:end,:); dfe = tmp(2:2:end,:);
                    elseif size(tmp,1) == size(Yr,1)
                        df  = tmp(1,:); dfe = tmp(2,:);
                    end
                    LIMO.model = rmfield(LIMO.model,'conditions_df');
                    LIMO.model.conditions_df = [df dfe]; clear tmp
                end
                
                if LIMO.design.nb_interactions ~=0
                    tmp =  cell2mat(LIMO.model.interactions_df)'; % dim (elec*[df dfe]) * 1
                    if size(tmp,1) == size(Yr,1)*2 
                        df  = tmp(1:2:end,:); dfe = tmp(2:2:end,:);
                    elseif size(tmp,1) == size(Yr,1)
                        df  = tmp(1,:); dfe = tmp(2,:);
                    end
                    LIMO.model = rmfield(LIMO.model,'interactions_df');
                    LIMO.model.interactions_df = [df dfe]; clear tmp
                end
                
                if LIMO.design.nb_continuous ~=0
                    tmp = cell2mat(LIMO.model.continuous_df)'; % dim (elec*[df dfe]) * n
                    if size(tmp,1) == size(Yr,1)*2 
                        df  = tmp(1:2:end,1); dfe = tmp(2:2:end,:);
                    elseif size(tmp,1) == size(Yr,1)
                        df  = tmp(:,1); dfe = tmp(:,2:end);
                    end
                    LIMO.model = rmfield(LIMO.model,'continuous_df');
                    LIMO.model.continuous_df = [df dfe]; clear tmp
                end
            end
        end
        
        % update the files to be stored on the disk
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            if strcmp(LIMO.design.method,'IRLS')
                W(channel,:,index) = model.W';
                for ft=size(W,2):-1:1 % each freq*time has different weighting
                    WX = LIMO.design.X .* repmat(squeeze(W(channel,ft,:)),1,size(X,2));
                    fitted_data(:,ft)  = (WX*squeeze(model.betas(:,ft,:)));
                end
            elseif strcmp(LIMO.design.method,'WLS')
                W(channel,:,index) = model.W'; 
                for f=n_freqs:-1:1 % each freq has different weighting
                    WX = LIMO.design.X .* repmat(squeeze(W(channel,f,:)),1,size(X,2));
                    fitted_data(1,f,:,:)  = (WX*squeeze(model.betas(:,f,:)))';
                end
                fitted_data = squeeze(limo_tf_4d_reshape(fitted_data))';
                % reshape beta freq to ft
                for c=size(model.betas,1):-1:1
                    tmp(c,:) = reshape(model.betas(c,:,:), [n_freqs*n_times,1]);
                end
                model.betas = tmp; clear tmp

            else % OLS, W is already ones
                fitted_data = LIMO.design.X*model.betas;
            end
        else
            if strcmp(LIMO.design.method,'IRLS')
                W(channel,:,index) = model.W';
            elseif strcmp(LIMO.design.method,'WLS')
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    W(channel,:,index) = model.W;
                else
                    W(channel,index) = model.W;
                end
            end
            fitted_data = LIMO.design.X*model.betas;
        end
        
        % all these always 3D - reshape before saving
        Yhat(channel,:,index) = fitted_data';
        Res(channel,:,index)  = squeeze(Yr(channel,:,index)) - fitted_data';
        clear fitted_data
        R2(channel,:,1)       = model.R2_univariate;
        R2(channel,:,2)       = model.F;
        R2(channel,:,3)       = model.p;
        Betas(channel,:,:)    = model.betas';

        if prod(LIMO.design.nb_conditions) ~=0
            if length(LIMO.design.nb_conditions) == 1
                tmp_Condition_effect(channel,:,1,1) = model.conditions.F;
                tmp_Condition_effect(channel,:,1,2) = model.conditions.p;
            else
                for i=1:length(LIMO.design.nb_conditions)
                    tmp_Condition_effect(channel,:,i,1) = model.conditions.F(i,:);
                    tmp_Condition_effect(channel,:,i,2) = model.conditions.p(i,:);
                end
            end
        end
        
        if LIMO.design.fullfactorial == 1
            for i=1:length(LIMO.design.nb_interactions)
                tmp_Interaction_effect(channel,:,i,1) = model.interactions.F(i,:);
                tmp_Interaction_effect(channel,:,i,2) = model.interactions.p(i,:);
            end
        end
        
        if LIMO.design.nb_continuous ~=0
            if LIMO.design.nb_continuous == 1
                tmp_Covariate_effect(channel,:,1,1) = model.continuous.F;
                tmp_Covariate_effect(channel,:,1,2) = model.continuous.p;
            else
                for i=1:LIMO.design.nb_continuous
                    tmp_Covariate_effect(channel,:,i,1) = model.continuous.F(:,i);
                    tmp_Covariate_effect(channel,:,i,2) = model.continuous.p(:,i);
                end
            end
        end
        clear model
    end
    warning on;
    
    % save data on the disk and clean out
    disp('saving data to disk')
    LIMO.design.X       = X;
    LIMO.design.weights = W;
    LIMO.design.status  = 'done';
    if ~isfield(LIMO.design,'name')
        LIMO.design.name    = 'GLM';
    end
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO'); 
    
    if strcmpi(LIMO.Analysis,'Time-Frequency')
         Yhat  = limo_tf_4d_reshape(Yhat);
         Res   = limo_tf_4d_reshape(Res);
         R2    = limo_tf_4d_reshape(R2);
         Betas = limo_tf_4d_reshape(Betas);
    end
    save(fullfile(LIMO.dir,'Yhat.mat'),  'Yhat',  '-v7.3');
    save(fullfile(LIMO.dir,'Res.mat'),   'Res',   '-v7.3'); 
    save(fullfile(LIMO.dir,'Betas.mat'), 'Betas', '-v7.3');
    save(fullfile(LIMO.dir,'R2.mat'),    'R2',    '-v7.3');
    clear Yhat Res Betas R2
    
    if prod(LIMO.design.nb_conditions) ~=0
        for i=1:length(LIMO.design.nb_conditions)
            name = sprintf('Condition_effect_%g.mat',i);
            if size(tmp_Condition_effect,1) == 1
                tmp                     = squeeze(tmp_Condition_effect(1,:,i,:));
                Condition_effect        = NaN(1,size(tmp_Condition_effect,2),2);
                Condition_effect(1,:,:) = tmp;
            else
                Condition_effect        = squeeze(tmp_Condition_effect(:,:,i,:));
            end
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                Condition_effect = limo_tf_4d_reshape(Condition_effect);
            end
            save(fullfile(LIMO.dir,name),'Condition_effect','-v7.3')
        end
        clear Condition_effect tmp_Condition_effect
    end
    
    if LIMO.design.fullfactorial == 1
        for i=1:length(LIMO.design.nb_interactions)
            name = sprintf('Interaction_effect_%g.mat',i);
            if size(tmp_Interaction_effect,1) == 1
                tmp                       = squeeze(tmp_Interaction_effect(1,:,i,:));
                Interaction_effect        = NaN(1,size(tmp_Interaction_effect,2),2);
                Interaction_effect(1,:,:) = tmp;
            else
                Interaction_effect        = squeeze(tmp_Interaction_effect(:,:,i,:));
            end
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                Interaction_effect = limo_tf_4d_reshape(Interaction_effect);
            end
            save(fullfile(LIMO.dir,name),'Interaction_effect','-v7.3')
        end
        clear Interaction_effect tmp_Interaction_effect
    end
    
    if LIMO.design.nb_continuous ~=0
        for i=1:LIMO.design.nb_continuous
            name = sprintf('Covariate_effect_%g.mat',i);
            if size(tmp_Covariate_effect,1) == 1
                tmp                     = squeeze(tmp_Covariate_effect(1,:,i,:));
                Covariate_effect        = NaN(1,size(tmp_Covariate_effect,2),2);
                Covariate_effect(1,:,:) = tmp;
            else
                Covariate_effect        = squeeze(tmp_Covariate_effect(:,:,i,:));
            end
            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                Covariate_effect = limo_tf_4d_reshape(Covariate_effect);
            end
            save(fullfile(LIMO.dir,name),'Covariate_effect','-v7.3')
        end
        clear Covariate_effect tmp_Covariate_effect
    end
    clear file channel filename model reg dir i W
end


%% Bootstrap under H0
% ----------------------------------------------------------
% ----------------------------------------------------------

if LIMO.design.bootstrap ~=0
    
    % being pretty time consuming, let users have one core free
    if isempty(gcp('nocreate'))
        parpool(feature('numCores')-1);
    end
    
    % avoid overwriting / recomputing H0 if done
    % (limo_eeg(4) called via the results interface)
    if exist('H0','dir')
        if ~exist('TFCE','dir') && LIMO.design.tfce == 1
            overwrite_H0boot = questdlg('H0 present for tfce, overwrite?','limo check','yes','no','no');
        else
            overwrite_H0boot = questdlg('overwrite H0?','limo check','yes','no','yes');
            if strcmp(overwrite_H0boot,'no') || isempty(overwrite_H0boot)
                if exist('warndlg2','file')
                    warndlg2('Analysis stopped - not overwriting H0')
                else
                    warndlg('Analysis stopped - not overwriting H0')
                end
                return
            end
        end
    else
        overwrite_H0boot = 'yes';
    end
    
    if strcmp(overwrite_H0boot,'yes')
        try
            mkdir H0;
            fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Bootstrapping GLM, ... \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
            
            Yr = load('Yr');
            Yr = Yr.Yr; % reload in any cases - ensuring right dimensions
            if size(Yr,1) == 1 % in any cases, just one channel/component
                array = 1;
            else
                if LIMO.Level == 2        % second level we can have missing subjects because of
                    array = 1:size(Yr,1); % bad channels for some subjects - just adjust X
                else % level 1 = skip empty channels - can't be missing trials
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        array = find(~isnan(Yr(:,1,1,1)));
                    else
                        array = find(~isnan(Yr(:,1,1)));
                    end
                end
            end
            
            if LIMO.design.bootstrap < 800
                if LIMO.design.bootstrap == 101
                    fprintf('bootstrap set to 101, this is a testing hack, otherwise the minimum required would be 800\n')
                else
                    fprintf('setting bootstrap to the minimum required, i.e. 800 instead of %g\n',LIMO.design.bootstrap)
                    LIMO.design.bootstrap = 800;
                end
            end
            nboot = LIMO.design.bootstrap;
            
            if LIMO.Level == 2
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    boot_table = limo_create_boot_table(squeeze(Yr(:,1,:,:)),nboot);
                else
                    boot_table = limo_create_boot_table(Yr,nboot);
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    boot_table = randi(size(Yr,4),size(Yr,4),nboot);
                else
                    boot_table = randi(size(Yr,3),size(Yr,3),nboot);
                end
            end
            
            % make file of the right size to avoid reshaping 5D files
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                H0_R2    = NaN(size(Yr,1), size(Yr,2), size(Yr,3), 3, nboot); % stores R, F and p values for each boot
                H0_Betas = NaN(size(Yr,1), size(Yr,2), size(Yr,3), size(LIMO.design.X,2), nboot);
                if LIMO.design.nb_conditions ~= 0
                    tmp_H0_Conditions = NaN(size(Yr,1), size(Yr,2), size(Yr,3), length(LIMO.design.nb_continuous), 2, nboot);
                end
                if LIMO.design.nb_interactions ~=0
                    tmp_H0_Interaction_effect = NaN(size(Yr,1), size(Yr,2),  size(Yr,3),length(LIMO.design.nb_interactions), 2, nboot);
                end
                if LIMO.design.nb_continuous ~= 0
                    tmp_H0_Covariates = NaN(size(Yr,1), size(Yr,2), size(Yr,3), LIMO.design.nb_continuous, 2, nboot);
                end
            else
                H0_R2    = NaN(size(Yr,1), size(Yr,2), 3, nboot); % stores R, F and p values for each boot
                H0_Betas = NaN(size(Yr,1), size(Yr,2), size(LIMO.design.X,2), nboot);
                if LIMO.design.nb_conditions ~= 0
                    tmp_H0_Conditions = NaN(size(Yr,1), size(Yr,2), length(LIMO.design.nb_continuous), 2, nboot);
                end
                if LIMO.design.nb_interactions ~=0
                    tmp_H0_Interaction_effect = NaN(size(Yr,1), size(Yr,2), length(LIMO.design.nb_interactions), 2, nboot);
                end
                if LIMO.design.nb_continuous ~= 0
                    tmp_H0_Covariates = NaN(size(Yr,1), size(Yr,2), LIMO.design.nb_continuous, 2, nboot);
                end
            end
            
            % run the analysis, loop per channel
            % limo_glm_boot then uses parfor to bootstrap the data under the null
            warning off;
            X = LIMO.design.X;
            h = waitbar(0,'bootstraping data','name','% done');
            
            for e = 1:length(array)
                channel = array(e);
                waitbar(e/size(array,2))
                fprintf('bootstrapping channel %g \n',channel);
                if LIMO.Level == 2
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
                            Y = squeeze(Yr(channel,:,:,:));
                            index = find(~isnan(Y(1,1,:))); % because across subjects, we can have missing data
                            for f=1:size(Yr,2)
                                Weights = squeeze(LIMO.design.weights(channel,f,index));
                                model{f} = limo_glm_boot(squeeze(Y(f,:,index))',X(index,:), Weights,...
                                    LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,...
                                    LIMO.design.method,LIMO.Analysis,boot_table{channel});
                            end
                        elseif strcmp(LIMO.design.method,'IRLS')
                            Y = squeeze(Yr(channel,:,:,:));
                            index = find(~isnan(Y(1,1,:))); % because across subjects, we can have missing data
                            if numel(size(LIMO.design.weights)) == 3
                                LIMO.design.weights = limo_tf_4d_reshape(LIMO.design.weights,LIMO.data.size4D);
                            end
                            for f=1:size(Yr,2)
                                Weights = squeeze(LIMO.design.weights(channel,f,:,index));
                                model{f} = limo_glm_boot(squeeze(Y(f,:,index))',X(index,:), Weights,...
                                    LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,...
                                    LIMO.design.method,LIMO.Analysis,boot_table{channel});
                            end
                        end
                    else
                        if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
                            Y = squeeze(Yr(channel,:,:));
                            index = find(~isnan(Y(1,:)));
                            Weights = squeeze(LIMO.design.weights(channel,index))';
                            model = limo_glm_boot(squeeze(Y(:,index))',X(index,:),Weights,...
                                LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,...
                                LIMO.design.method,LIMO.Analysis,boot_table{channel});
                        elseif strcmp(LIMO.design.method,'IRLS')
                            Y = squeeze(Yr(channel,:,:));
                            index = find(~isnan(Y(1,:)));
                            Weights = squeeze(LIMO.design.weights(channel,:,index));
                            model = limo_glm_boot(squeeze(Y(:,index))',X(index,:),Weights,...
                                LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,...
                                LIMO.design.method,LIMO.Analysis,boot_table{channel});
                        end
                    end
                else % LIMO.Level == 1
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
                            for f=1:size(Yr,2)
                                LIMO.Weights = squeeze(LIMO.design.weights(channel,f,:));
                                model{f} = limo_glm_boot(squeeze(Yr(channel,f,:,:))',LIMO,boot_table);
                            end
                        elseif strcmp(LIMO.design.method,'IRLS')
                            if numel(size(LIMO.design.weights)) == 3
                                LIMO.design.weights = limo_tf_4d_reshape(LIMO.design.weights,LIMO.data.size4D);
                            end
                            for f=1:size(Yr,2)
                                LIMO.Weights = squeeze(LIMO.design.weights(channel,f,:,:));
                                model{f} = limo_glm_boot(squeeze(Yr(channel,f,:,:))',LIMO,boot_table);
                            end
                        end
                    else
                        if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
                            LIMO.Weights = squeeze(LIMO.design.weights(channel,:))';
                            model        = limo_glm_boot(squeeze(Yr(channel,:,:))',LIMO,boot_table);
                        elseif strcmp(LIMO.design.method,'IRLS')
                            LIMO.Weights = squeeze(LIMO.design.weights(channel,:,:));
                            model        = limo_glm_boot(squeeze(Yr(channel,:,:))',LIMO,boot_table);
                        end
                    end
                end
                
                % update the files to be stored on the disk
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    for f=1:length(model)
                        for B = 1:nboot % now loop because we use cells
                            H0_Betas(channel,f,:,:,B) = model{f}.betas{B};
                            H0_R2(channel,f,:,1,B)  = model{f}.R2_univariate{B};
                            H0_R2(channel,f,:,2,B)  = model{f}.F{B};
                            H0_R2(channel,f,:,3,B)  = model{f}.p{B};
                            
                            if prod(LIMO.design.nb_conditions) ~=0
                                if length(LIMO.design.nb_conditions) == 1
                                    tmp_H0_Conditions(channel,f,:,1,1,B) = model{f}.conditions.F{B};
                                    tmp_H0_Conditions(channel,f,:,1,2,B) = model{f}.conditions.p{B};
                                else
                                    for i=1:length(LIMO.design.nb_conditions)
                                        tmp_H0_Conditions(channel,f,:,i,1,B) = model{f}.conditions.F{B}(i,:);
                                        tmp_H0_Conditions(channel,f,:,i,2,B) = model{f}.conditions.p{B}(i,:);
                                    end
                                end
                            end
                            
                            if LIMO.design.fullfactorial == 1
                                if length(LIMO.design.nb_interactions) == 1
                                    tmp_H0_Interaction_effect(channel,f,:,1,1,B) = model{f}.interactions.F{B};
                                    tmp_H0_Interaction_effect(channel,f,:,1,2,B) = model{f}.interactions.p{B};
                                else
                                    for i=1:length(LIMO.design.nb_interactions)
                                        tmp_H0_Interaction_effect(channel,f,:,i,1,B) = model{f}.interactions.F{B}(i,:);
                                        tmp_H0_Interaction_effect(channel,f,:,i,2,B) = model{f}.interactions.p{B}(i,:);
                                    end
                                end
                            end
                            
                            if LIMO.design.nb_continuous ~=0
                                if LIMO.design.nb_continuous == 1
                                    tmp_H0_Covariates(channel,f,:,1,1,B) = model{f}.continuous.F{B};
                                    tmp_H0_Covariates(channel,f,:,1,2,B) = model{f}.continuous.p{B};
                                else
                                    for i=1:LIMO.design.nb_continuous
                                        tmp_H0_Covariates(channel,f,:,i,1,B) = model{f}.continuous.F{B}(:,i);
                                        tmp_H0_Covariates(channel,f,:,i,2,B) = model{f}.continuous.p{B}(:,i);
                                    end
                                end
                            end
                        end
                    end
                else % erp or spec
                    for B = 1:nboot
                        H0_Betas(channel,:,:,B) = model.betas{B};
                        H0_R2(channel,:,1,B)    = model.R2_univariate{B};
                        H0_R2(channel,:,2,B)    = model.F{B};
                        H0_R2(channel,:,3,B)    = model.p{B};
                        
                        if prod(LIMO.design.nb_conditions) ~=0
                            if length(LIMO.design.nb_conditions) == 1
                                tmp_H0_Conditions(channel,:,1,1,B) = model.conditions.F{B};
                                tmp_H0_Conditions(channel,:,1,2,B) = model.conditions.p{B};
                            else
                                for i=1:length(LIMO.design.nb_conditions)
                                    tmp_H0_Conditions(channel,:,i,1,B) = model.conditions.F{B}(i,:);
                                    tmp_H0_Conditions(channel,:,i,2,B) = model.conditions.p{B}(i,:);
                                end
                            end
                        end
                        
                        if LIMO.design.fullfactorial == 1
                            if length(LIMO.design.nb_interactions) == 1
                                tmp_H0_Interaction_effect(channel,:,1,1,B) = model.interactions.F{B};
                                tmp_H0_Interaction_effect(channel,:,1,2,B) = model.interactions.p{B};
                            else
                                for i=1:length(LIMO.design.nb_interactions)
                                    tmp_H0_Interaction_effect(channel,:,i,1,B) = model.interactions.F{B}(i,:);
                                    tmp_H0_Interaction_effect(channel,:,i,2,B) = model.interactions.p{B}(i,:);
                                end
                            end
                        end
                        
                        if LIMO.design.nb_continuous ~=0
                            if LIMO.design.nb_continuous == 1
                                tmp_H0_Covariates(channel,:,1,1,B) = model.continuous.F{B};
                                tmp_H0_Covariates(channel,:,1,2,B) = model.continuous.p{B};
                            else
                                for i=1:LIMO.design.nb_continuous
                                    if all(size(squeeze(tmp_H0_Covariates(channel,:,i,1,B))) == size(squeeze(model.continuous.F{B}(:,i)))) || ...
                                            all(size(squeeze(tmp_H0_Covariates(channel,:,i,1,B))) == size(squeeze(model.continuous.F{B}(:,i))'))
                                        tmp_H0_Covariates(channel,:,i,1,B) = model.continuous.F{B}(:,i);
                                        tmp_H0_Covariates(channel,:,i,2,B) = model.continuous.p{B}(:,i);
                                    else
                                        tmp_H0_Covariates(channel,:,i,1,B) = model.continuous.F{B}(i,:);
                                        tmp_H0_Covariates(channel,:,i,2,B) = model.continuous.p{B}(i,:);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            close(h); warning on;
            clear Yr
            
            % save data on the disk and clear out
            save([LIMO.dir filesep 'H0' filesep 'H0_R2.mat'],'H0_R2','-v7.3');
            save([LIMO.dir filesep 'H0' filesep 'boot_table.mat'],'boot_table');
            save([LIMO.dir filesep 'H0' filesep 'H0_Betas.mat'],'H0_Betas','-v7.3');
            clear H0_R2 boot_table H0_Betas
            
            if prod(LIMO.design.nb_conditions) ~=0
                for i=1:length(LIMO.design.nb_conditions)
                    name = sprintf('H0_Condition_effect_%g',i);
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tmp = squeeze(tmp_H0_Conditions(:,:,:,i,:,:));
                    else
                        tmp = squeeze(tmp_H0_Conditions(:,:,i,:,:));
                    end
                    
                    if isfield(LIMO.design,'electrode')
                        if ~isempty(LIMO.design.electrode)
                            H0_Condition_effect = NaN([1 size(tmp)]);
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                H0_Condition_effect(1,:,:,:,:) = tmp;
                            else
                                H0_Condition_effect(1,:,:,:) = tmp;
                            end
                        else
                            H0_Condition_effect = tmp;
                        end
                    else
                        H0_Condition_effect = tmp;
                    end
                    save(fullfile(LIMO.dir,['H0' filesep name]),'H0_Condition_effect','-v7.3');
                    clear tmp H0_Condition_effect
                end
                clear tmp_H0_Conditions
            end
            
            if LIMO.design.fullfactorial == 1
                for i=1:length(LIMO.design.nb_interactions)
                    name = sprintf('H0_Interaction_effect_%g',i);
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tmp = squeeze(tmp_H0_Interaction_effect(:,:,:,i,:,:));
                    else
                        tmp = squeeze(tmp_H0_Interaction_effect(:,:,i,:,:));
                    end
                    
                    if isfield(LIMO.design,'electrode')
                        if ~isempty(LIMO.design.electrode)
                            H0_Interaction_effect = NaN([1 size(tmp)]);
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                H0_Interaction_effect(1,:,:,:,:) = tmp;
                            else
                                H0_Interaction_effect(1,:,:,:) = tmp;
                            end
                        else
                            H0_Interaction_effect = tmp;
                        end
                    else
                        H0_Interaction_effect = tmp;
                    end
                    save(fullfile(LIMO.dir,['H0' filesep name]),'H0_Interaction_effect','-v7.3');
                    clear H0_Interaction_effect
                end
                clear tmp_H0_Interaction_effect
            end
            
            if LIMO.design.nb_continuous ~=0
                for i=1:LIMO.design.nb_continuous
                    name = sprintf('H0_Covariate_effect_%g',i);
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tmp = squeeze(tmp_H0_Covariates(:,:,:,i,:,:));
                    else
                        tmp = squeeze(tmp_H0_Covariates(:,:,i,:,:));
                    end
                    
                    if isfield(LIMO.design,'electrode')
                        if ~isempty(LIMO.design.electrode)
                            H0_Covariate_effect = NaN([1 size(tmp)]);
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                H0_Covariate_effect(1,:,:,:,:) = tmp;
                            else
                                H0_Covariate_effect(1,:,:,:) = tmp;
                            end
                        else
                            H0_Covariate_effect = tmp;
                        end
                    else
                        H0_Covariate_effect = tmp;
                    end
                    save(fullfile(LIMO.dir,['H0' filesep name]),'H0_Covariate_effect','-v7.3');
                    clear tmp H0_Covariate_effect
                end
                clear tmp_H0_Covariates
            end
            
            clear channel model H0_R2;
            cd(LIMO.dir); disp(' ');
            
        catch boot_error
            disp('an error occured while attempting to bootstrap the data')
            error('%s \n',boot_error.message);
        end
    end
end

%% TFCE 
% --------------
if LIMO.design.tfce == 1
    
    if exist('TFCE','dir')
        if strcmp(questdlg('TFCE directory detected, overwrite?','data check','Yes','No','No'),'No')
            warndlg2('Analysis stopped - not overwriting TFCE')
            return
        end
    end
    
    % check if there is a neighbouring matrix 
    % (since TFCE integrates over clusters)
    if ~isfield(LIMO.data,'neighbouring_matrix')
       warning('no neighbouring matrix found, this is required for TFCE')
       answer = questdlg('load or compute neighbouring matrix?','channel neighbouring definition','Load','Compute','Compute');
       if strcmp(answer,'Load')
           [file,newpath,whatsup] = uigetfile('*.mat','select neighbourghing matrix (or expected chanloc file)');
           if whatsup == 0
               disp('selection aborded');
               return
           else
               tmp   = load(fullfile(newpath,file));
               fn    = fieldnames(tmp);
               index = find(ismember(fn,'channeighbstructmat'));
               if isempty(index)
                   error('no neighbouring matrix ''channeighbstructmat'' found')
               else
                   LIMO.data.neighbouring_matrix = getfield(tmp,fn{index});
                   save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
               end
           end
       else
           [~, LIMO.data.neighbouring_matrix] = limo_expected_chanlocs(LIMO.data.data, LIMO.data.data_dir);
           if isempty(LIMO.data.neighbouring_matrix)
               error('no neighbouring matrix returned, try creating with limo tools')
           else
               save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
           end
       end
    end
       
    fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Computing TFCE for GLM takes a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
    mkdir tfce; 
    
    % R2
    limo_tfce_handling(fullfile(LIMO.dir,'R2.mat'),'checkfile','no')
    
    % conditions
    if prod(LIMO.design.nb_conditions) ~=0
        for i=1:length(LIMO.design.nb_conditions)
            name = sprintf('Condition_effect_%g.mat',i);
            limo_tfce_handling(fullfile(LIMO.dir,name),'checkfile','no')
        end
    end
    
    % interactions
    if LIMO.design.fullfactorial == 1
        for i=1:length(LIMO.design.fullfactorial)
            name = sprintf('Interaction_effect_%g.mat',i);
            limo_tfce_handling(fullfile(LIMO.dir,name),'checkfile','no')
        end
    end
    
    % covariates / continuous regressors
    if LIMO.design.nb_continuous ~=0
        for i=1:LIMO.design.nb_continuous
            name = sprintf('Covariate_effect_%g.mat',i); 
            limo_tfce_handling(fullfile(LIMO.dir,name),'checkfile','no')
        end
    end
end




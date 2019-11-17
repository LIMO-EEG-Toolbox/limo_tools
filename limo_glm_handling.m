function limo_glm_handling(LIMO)

% data handling function for GLM
% this function calls limo_glm, limo_glm_boot to get the analysis done, and
% organize all the files around that - externalized from limo_eeg(4)
%
% FORMAT limo_glm_handling(LIMO)
%
% Cyril Pernet
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2019

cd(LIMO.dir);

%% Compute GLM and save stats files

if strcmp(LIMO.design.status,'to do')
    Yr    = load('Yr');    Yr    = Yr.(cell2mat(fieldnames(Yr)));  
    Yhat  = load('Yhat');  Yhat  = Yhat.(cell2mat(fieldnames(Yhat)));
    Res   = load('Res');   Res   = Res.(cell2mat(fieldnames(Res)));
    R2    = load('R2');    R2    = R2.(cell2mat(fieldnames(R2)));
    Betas = load('Betas'); Betas = Betas.(cell2mat(fieldnames(Betas)));
  
    % check dimensions (3D vs 4D) 
    % --------------------------------------
    if strcmpi(LIMO.Analysis,'Time-Frequency') 
         [~,n_freqs,~,~] = size(Yr);
         Yr    = limo_tf_4d_reshape(Yr); % reshape to 3D
         Yhat  = limo_tf_4d_reshape(Yhat);
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
    
    % -------------- loop the analysis electrode per electrode
    if size(Yr,1) == 1
        array = 1;
    else
        array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
    end
    
    % ------------- prepare weight matrix  -------------------------------------
    try
        limo_pcout(squeeze(Yr(array(1),:,:))');
    catch pcout_error
        if strcmp(pcout_error,'Principal Component Projection cannot be computed, more observations than variables are needed')
            LIMO.design.method = 'OLS'; disp('Cannot use WLS, not enough observations - switching to OLS')
        end
    end
    
    
    if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            W = ones(size(Yr,1),n_freqs,size(Yr,3));
        else
            W = ones(size(Yr,1),size(Yr,3));
        end
    elseif strcmp(LIMO.design.method,'IRLS')
        W = ones(size(Yr));
    end
    
    % ------------ run limo_glm per electrodes ---------------------------
    update = 1;
    X = LIMO.design.X;
    for e = 1:size(array,1)
        electrode = array(e); warning off;
        if LIMO.Level == 2
            fprintf('analyzing channel %g/%g \n',e,size(array,1));
            Y             = squeeze(Yr(electrode,:,:));
            index         = find(~isnan(Y(1,:)));
            if isempty(index)
                index     = 1:size(Y,2);
            end
            Y             = Y(:,index);
            LIMO.design.X = X(index,:);
            model         = limo_glm(Y',LIMO); warning on;
        else % level 1 we should not have any NaNs
            if strcmp(LIMO.Type,'Channels')
                fprintf('analyzing channel %g/%g \n',e,size(array,1));
            else
                fprintf('analyzing component %g/%g \n',e,size(array,1));
            end
            index = 1:size(Yr,3);
            model = limo_glm(squeeze(Yr(electrode,:,:))',LIMO);
        end
        
        % update the LIMO.mat (do it only once)
        if update == 1 && ~strcmpi(LIMO.design.method,'IRLS')
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
        elseif update == 1 &&  strcmpi(LIMO.design.method,'IRLS')
            LIMO.model.model_df(electrode,:) = model.df;
            if LIMO.design.nb_conditions ~=0
                LIMO.model.conditions_df(electrode,:,:)  = model.conditions.df;
            end
            if LIMO.design.nb_interactions ~=0
                LIMO.model.interactions_df(electrode,:,:)  = model.interactions.df;
            end
            if LIMO.design.nb_continuous ~=0
                LIMO.model.continuous_df(electrode,:,:)  = model.continuous.df;
            end
        end
        
        % update the files to be stored on the disk
        if strcmp(LIMO.design.method,'IRLS')
            W(electrode,:,index) = model.W';
        elseif strcmp(LIMO.design.method,'WLS')
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                W(electrode,:,index) = model.W;
            else
                W(electrode,index) = model.W;
            end
        end
        fitted_data             = LIMO.design.X*model.betas;
        Yhat(electrode,:,index) = fitted_data';
        Res(electrode,:,index)  = squeeze(Yr(electrode,:,index)) - fitted_data';
        clear fitted_data
        R2(electrode,:,1)       = model.R2_univariate;
        R2(electrode,:,2)       = model.F;
        R2(electrode,:,3)       = model.p;
        Betas(electrode,:,:)    = model.betas';
        
        if prod(LIMO.design.nb_conditions) ~=0
            if length(LIMO.design.nb_conditions) == 1
                tmp_Condition_effect(electrode,:,1,1) = model.conditions.F;
                tmp_Condition_effect(electrode,:,1,2) = model.conditions.p;
            else
                for i=1:length(LIMO.design.nb_conditions)
                    tmp_Condition_effect(electrode,:,i,1) = model.conditions.F(i,:);
                    tmp_Condition_effect(electrode,:,i,2) = model.conditions.p(i,:);
                end
            end
        end
        
        if LIMO.design.fullfactorial == 1
            for i=1:length(LIMO.design.nb_interactions)
                tmp_Interaction_effect(electrode,:,i,1) = model.interactions.F(i,:);
                tmp_Interaction_effect(electrode,:,i,2) = model.interactions.p(i,:);
            end
        end
        
        if LIMO.design.nb_continuous ~=0
            if LIMO.design.nb_continuous == 1
                tmp_Covariate_effect(electrode,:,1,1) = model.continuous.F;
                tmp_Covariate_effect(electrode,:,1,2) = model.continuous.p;
            else
                for i=1:LIMO.design.nb_continuous
                    tmp_Covariate_effect(electrode,:,i,1) = model.continuous.F(:,i);
                    tmp_Covariate_effect(electrode,:,i,2) = model.continuous.p(:,i);
                end
            end
        end
    end
    
    % save data on the disk and clean out
    disp('saving data to disk')
    LIMO.design.X       = X;
    LIMO.design.weights = W;
    LIMO.design.status  = 'done';
    LIMO.design.name    = 'GLM';
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
    clear file electrode filename model reg dir i W
end


%% Bootstrap under H0
% -------------------------------
if LIMO.design.bootstrap ~=0
    % avoid overwriting / recomputing H0 if done
    % (limo_eeg(4) called via the results interface)
    if exist('H0','dir')
        ow = questdlg('overwrite H0?','limo check','yes','no','yes');
        if strcmp(ow,'no') || isempty(ow)
            warndlg2('Analysis stopped - not overwriting H0')
            return
        end
    else
        try
            mkdir H0;
            fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Bootstrapping GLM, ... \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
            
            Yr = load('Yr'); 
            Yr = Yr.Yr; % reload in any cases - ensuring right dimensions
            if size(Yr,1) == 1
                array = 1;
            else
                array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
            end
            
            if LIMO.design.bootstrap <= 800
                if LIMO.design.bootstrap == 101
                    fprintf('bootstrap set to 101, this is a testing hack, otherwise the minimum required would be 800\n')
                else
                    fprintf('setting bootstrap to the minimum required, i.e. 800 instead of %g\n',LIMO.design.bootstrap)
                    LIMO.design.bootstrap = 800;
                end
            end
            nboot = LIMO.design.bootstrap;
            
            if LIMO.Level == 2
                boot_table = limo_create_boot_table(Yr,nboot);
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
            
            % run the analysis
            warning off;
            X = LIMO.design.X;
            h = waitbar(0,'bootstraping data','name','% done');
            for e = 1:size(array,1)
                electrode = array(e);
                waitbar(e/size(array,1))
                fprintf('bootstrapping electrode %g \n',electrode);
                if LIMO.Level == 2
                    Y = squeeze(Yr(electrode,:,:));
                    index = find(~isnan(Y(1,:))); % because across subjects, we can have missing data
                    model = limo_glm_boot(Y(:,index)',X(index,:),LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,LIMO.design.method,LIMO.Analysis,[],[],boot_table{electrode});
                else
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        for f=1:size(Yr,2)
                            model{f} = limo_glm_boot(squeeze(Yr(electrode,f,:,:))',LIMO,boot_table);
                        end
                    else
                        model = limo_glm_boot(squeeze(Yr(electrode,:,:))',LIMO,boot_table);
                    end
                end
                
                % update the files to be stored on the disk
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    for f=1:length(model)
                        for B = 1:nboot % now loop because we use cells
                            H0_Betas(electrode,f,:,:,B) = model{f}.betas{B};
                            H0_R2(electrode,f,:,1,B)  = model{f}.R2_univariate{B};
                            H0_R2(electrode,f,:,2,B)  = model{f}.F{B};
                            H0_R2(electrode,f,:,3,B)  = model{f}.p{B};
                            
                            if prod(LIMO.design.nb_conditions) ~=0
                                if length(LIMO.design.nb_conditions) == 1
                                    tmp_H0_Conditions(electrode,f,:,1,1,B) = model{f}.conditions.F{B};
                                    tmp_H0_Conditions(electrode,f,:,1,2,B) = model{f}.conditions.p{B};
                                else
                                    for i=1:length(LIMO.design.nb_conditions)
                                        tmp_H0_Conditions(electrode,f,:,i,1,B) = model{f}.conditions.F{B}(i,:);
                                        tmp_H0_Conditions(electrode,f,:,i,2,B) = model{f}.conditions.p{B}(i,:);
                                    end
                                end
                            end
                            
                            if LIMO.design.fullfactorial == 1
                                if length(LIMO.design.nb_interactions) == 1
                                    tmp_H0_Interaction_effect(electrode,f,:,1,1,B) = model{f}.interactions.F{B};
                                    tmp_H0_Interaction_effect(electrode,f,:,1,2,B) = model{f}.interactions.p{B};
                                else
                                    for i=1:length(LIMO.design.nb_interactions)
                                        tmp_H0_Interaction_effect(electrode,f,:,i,1,B) = model{f}.interactions.F{B}(i,:);
                                        tmp_H0_Interaction_effect(electrode,f,:,i,2,B) = model{f}.interactions.p{B}(i,:);
                                    end
                                end
                            end
                            
                            if LIMO.design.nb_continuous ~=0
                                if LIMO.design.nb_continuous == 1
                                    tmp_H0_Covariates(electrode,f,:,1,1,B) = model{f}.continuous.F{B};
                                    tmp_H0_Covariates(electrode,f,:,1,2,B) = model{f}.continuous.p{B};
                                else
                                    for i=1:LIMO.design.nb_continuous
                                        tmp_H0_Covariates(electrode,f,:,i,1,B) = model{f}.continuous.F{B}(:,i);
                                        tmp_H0_Covariates(electrode,f,:,i,2,B) = model{f}.continuous.p{B}(:,i);
                                    end
                                end
                            end
                        end
                    end
                else % erp or spec
                    for B = 1:nboot 
                        H0_Betas(electrode,:,:,B) = model.betas{B};
                        H0_R2(electrode,:,1,B)    = model.R2_univariate{B};
                        H0_R2(electrode,:,2,B)    = model.F{B};
                        H0_R2(electrode,:,3,B)    = model.p{B};
                        
                        if prod(LIMO.design.nb_conditions) ~=0
                            if length(LIMO.design.nb_conditions) == 1
                                tmp_H0_Conditions(electrode,:,1,1,B) = model.conditions.F{B};
                                tmp_H0_Conditions(electrode,:,1,2,B) = model.conditions.p{B};
                            else
                                for i=1:length(LIMO.design.nb_conditions)
                                    tmp_H0_Conditions(electrode,:,i,1,B) = model.conditions.F{B}(i,:);
                                    tmp_H0_Conditions(electrode,:,i,2,B) = model.conditions.p{B}(i,:);
                                end
                            end
                        end
                        
                        if LIMO.design.fullfactorial == 1
                            if length(LIMO.design.nb_interactions) == 1
                                tmp_H0_Interaction_effect(electrode,:,1,1,B) = model.interactions.F{B};
                                tmp_H0_Interaction_effect(electrode,:,1,2,B) = model.interactions.p{B};
                            else
                                for i=1:length(LIMO.design.nb_interactions)
                                    tmp_H0_Interaction_effect(electrode,:,i,1,B) = model.interactions.F{B}(i,:);
                                    tmp_H0_Interaction_effect(electrode,:,i,2,B) = model.interactions.p{B}(i,:);
                                end
                            end
                        end
                        
                        if LIMO.design.nb_continuous ~=0
                            if LIMO.design.nb_continuous == 1
                                tmp_H0_Covariates(electrode,:,1,1,B) = model.continuous.F{B};
                                tmp_H0_Covariates(electrode,:,1,2,B) = model.continuous.p{B};
                            else
                                for i=1:LIMO.design.nb_continuous
                                    tmp_H0_Covariates(electrode,:,i,1,B) = model.continuous.F{B}(:,i);
                                    tmp_H0_Covariates(electrode,:,i,2,B) = model.continuous.p{B}(:,i);
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
                        H0_Condition_effect = squeeze(tmp_H0_Conditions(:,:,:,i,:,:));
                    else
                        H0_Condition_effect = squeeze(tmp_H0_Conditions(:,:,i,:,:));
                    end
                    save(fullfile(LIMO.dir,['H0' filesep name]),'H0_Condition_effect','-v7.3');
                    clear H0_Condition_effect
                end
                clear tmp_H0_Conditions
            end
            
            if LIMO.design.fullfactorial == 1
                for i=1:length(LIMO.design.nb_interactions)
                    name = sprintf('H0_Interaction_effect_%g',i);
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        H0_Interaction_effect = squeeze(tmp_H0_Interaction_effect(:,:,:,i,:,:));
                    else
                        H0_Interaction_effect = squeeze(tmp_H0_Interaction_effect(:,:,i,:,:));
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
                        H0_Covariate_effect = squeeze(tmp_H0_Covariates(:,:,:,i,:,:));
                    else
                        H0_Covariate_effect = squeeze(tmp_H0_Covariates(:,:,i,:,:));
                    end
                    save(fullfile(LIMO.dir,['H0' filesep name]),'H0_Covariate_effect','-v7.3');
                    clear H0_Covariate_effect
                end
                clear tmp_H0_Covariates
            end
            
            clear electrode model H0_R2; 
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
       [~, LIMO.data.neighbouring_matrix] = limo_expected_chanlocs;
       if isempty(LIMO.data.neighbouring_matrix)
           return
       else
           save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
       end
    end
       
    fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Computing TFCE for GLM takes a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
    mkdir TFCE; neighbouring_matrix = LIMO.data.neighbouring_matrix;
    
    % R2
    R2 = load('R2'); 
    R2 = R2.(cell2mat(fieldnames(R2)));
    fprintf('Creating R2 TFCE scores \n');
    if size(R2,1) == 1
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score(1,:,:) = limo_tfce(2, squeeze(R2(:,:,:,2)),[]); % no neighbouring, time-freq cluster
        else
            tfce_score(1,:)   = limo_tfce(1, squeeze(R2(:,:,2)),neighbouring_matrix);
        end
    else
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score        = limo_tfce(3, squeeze(R2(:,:,:,2)),neighbouring_matrix);
        else
            tfce_score        = limo_tfce(2, squeeze(R2(:,:,2)),neighbouring_matrix);
        end
    end
    save(fullfile(LIMO.dir,['TFCE' filesep 'tfce_R2']),'tfce_score','-v7.3');
    clear R2;
    
    if exist(['H0' filesep 'H0_R2.mat'],'file')
        fprintf('Thresholding H0_R2 using TFCE \n');
        H0_R2 = load(['H0' filesep 'H0_R2.mat']);
        H0_R2 = H0_R2.(cell2mat(fieldnames(H0_R2)));
        if size(H0_R2,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(1,size(H0_R2,2),size(H0_R2,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,:,b) = limo_tfce(2,squeeze(H0_R2(:,:,:,2,b)),[],0);
                end
            else
                tfce_H0_score = NaN(1,size(H0_R2,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_R2(:,:,2,b)),neighbouring_matrix,0);
                end
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(size(H0_R2,1),size(H0_R2,2),size(H0_R2,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_R2(:,:,:,2,b)),neighbouring_matrix,0);
                end
            else
                tfce_H0_score = NaN(size(H0_R2,1),size(H0_R2,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_R2(:,:,2,b)),neighbouring_matrix,0);
                end
            end
        end
        save(fullfile(LIMO.dir,['H0' filesep 'tfce_H0_R2']),'tfce_H0_score','-v7.3'); 
        clear H0_R2 tfce_H0_score; 
    end
    
    % conditions
    if prod(LIMO.design.nb_conditions) ~=0
        for i=1:length(LIMO.design.nb_conditions)
            name = sprintf('Condition_effect_%g.mat',i); 
            fprintf('Creating %s TFCE score \n',name)
            Condition_effect = load(name);
            Condition_effect = Condition_effect.(cell2mat(fieldnames(Condition_effect)));
            if size(Condition_effect,1) == 1
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_score(1,:,:) = limo_tfce(2, squeeze(Condition_effect(:,:,:,1)),[]);                    
                else
                    tfce_score(1,:) = limo_tfce(1, squeeze(Condition_effect(:,:,1)),neighbouring_matrix);
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_score      = limo_tfce(3, squeeze(Condition_effect(:,:,:,1)),neighbouring_matrix);
                else
                    tfce_score      = limo_tfce(2, squeeze(Condition_effect(:,:,1)),neighbouring_matrix);
                end
            end
            save(fullfile(LIMO.dir,['TFCE' filesep 'tfce_' name]),'tfce_score','-v7.3');
            clear Condition_effect tfce_score; 
        end
        
        for i=1:length(LIMO.design.nb_conditions)
            name = sprintf('H0_Condition_effect_%g.mat',i);
            if exist(['H0' filesep name],'file')
                fprintf('Creating %s TFCE score \n',name);
                H0_Condition_effect = load(['H0' filesep name]);
                H0_Condition_effect = H0_Condition_effect.(cell2mat(fieldnames(H0_Condition_effect)));
                if size(H0_Condition_effect,1) == 1
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tfce_H0_score = NaN(1,size(H0_Condition_effect,2),size(H0_Condition_effect,3),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(1,:,:,b) = limo_tfce(2,squeeze(H0_Condition_effect(:,:,:,1,b)),[],0);
                        end
                    else
                        tfce_H0_score = NaN(1,size(H0_Condition_effect,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Condition_effect(:,:,1,b)),neighbouring_matrix,0);
                        end
                    end
                else
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_Condition_effect,2),size(H0_Condition_effect,3),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_Condition_effect(:,:,:,1,b)),neighbouring_matrix,0);
                        end
                    else
                        tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_Condition_effect,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Condition_effect(:,:,1,b)),neighbouring_matrix,0);
                        end
                    end
                end
                save(fullfile(LIMO.dir,['H0' filesep 'tfce_' name]),'tfce_H0_score','-v7.3');
                clear H0_Condition_effect tfce_H0_score;
            end
        end
    end
    
    % interactions
    if LIMO.design.fullfactorial == 1
        for i=1:length(LIMO.design.fullfactorial)
            name = sprintf('Interaction_effect_%g.mat',i);
            fprintf('Creating %s TFCE score \n',name)
            Interaction_effect = load(name);
            Interaction_effect = Interaction_effect.(cell2mat(fieldnames(Interaction_effect)));
            if size(Interaction_effect,1) == 1
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_score(1,:,:) = limo_tfce(2,squeeze(Interaction_effect(:,:,:,1)),[]);
                else
                    tfce_score(1,:) = limo_tfce(1,squeeze(Interaction_effect(:,:,1)),neighbouring_matrix);
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_score = limo_tfce(3,squeeze(Interaction_effect(:,:,:,1)),neighbouring_matrix);
                else
                    tfce_score = limo_tfce(2,squeeze(Interaction_effect(:,:,1)),neighbouring_matrix);
                end
            end
            save(fullfile(LIMO.dir,['TFCE' filesep 'tfce_' name]),'tfce_score','-v7.3');
            clear Interaction_effect tfce_score; 
        end
        
        for i=1:length(LIMO.design.fullfactorial)
            name = sprintf('H0_Interaction_effect_%g.mat',i);
            if exist(['H0' filesep name],'file')
                fprintf('Creating %s TFCE score \n',name);
                H0_Interaction_effect = load(['H0' filesep name]);
                H0_Interaction_effect = H0_Interaction_effect.(cell2mat(fieldnames(H0_Interaction_effect)));
                if size(H0_Interaction_effect,1) == 1
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tfce_H0_score = NaN(1,size(H0_Interaction_effect,2),size(H0_Interaction_effect,3),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(1,:,:,b) = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,:,1,b)),[],0);
                        end
                    else
                        tfce_H0_score = NaN(1,size(H0_Interaction_effect,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Interaction_effect(:,:,1,b)),neighbouring_matrix,0);
                        end
                    end
                else
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tfce_H0_score = NaN(size(H0_Interaction_effect,1),size(H0_Interaction_effect,2),size(H0_Interaction_effect,3),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_Interaction_effect(:,:,:,1,b)),neighbouring_matrix,0);
                        end
                    else
                        tfce_H0_score = NaN(size(H0_Interaction_effect,1),size(H0_Interaction_effect,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,1,b)),neighbouring_matrix,0);
                        end
                    end
                end
                save(fullfile(LIMO.dir,['H0' filesep 'tfce_' name]),'tfce_H0_score','-v7.3');
                clear H0_Interaction_effect tfce_H0_score;
            end
        end
    end
    
    % covariates / continuous regressors
    if LIMO.design.nb_continuous ~=0
        for i=1:LIMO.design.nb_continuous
            name = sprintf('Covariate_effect_%g.mat',i); 
            fprintf('Creating %s TFCE score \n',name)
            Covariate_effect = load(name);
            Covariate_effect = Covariate_effect.(cell2mat(fieldnames(Covariate_effect)));
            if size(Covariate_effect,1) == 1
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_score(1,:,:) = limo_tfce(2,squeeze(Covariate_effect(:,:,:,1)),[]);
                else
                    tfce_score(1,:) = limo_tfce(1,squeeze(Covariate_effect(:,:,1)),neighbouring_matrix);
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_score = limo_tfce(3,squeeze(Covariate_effect(:,:,:,1)),neighbouring_matrix);
                else
                    tfce_score = limo_tfce(2,squeeze(Covariate_effect(:,:,1)),neighbouring_matrix);
                end
            end
            save(fullfile(LIMO.dir,['TFCE' filesep 'tfce_' name]),'tfce_score','-v7.3');
            clear Covariate_effect tfce_score; 
        end
        
        fprintf('Creating H0 Covariate(s) TFCE scores \n');
        for i=1:LIMO.design.nb_continuous
            name = sprintf('H0_Covariate_effect_%g.mat',i); 
            if exist(['H0' filesep name],'file')
                fprintf('Creating %s TFCE score \n',name);
                H0_Covariate_effect = load(['H0' filesep name]);
                H0_Covariate_effect = H0_Covariate_effect.(cell2mat(fieldnames(H0_Covariate_effect)));
                if size(H0_Covariate_effect,1) == 1
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),size(H0_Covariate_effect,3),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(1,:,:,b) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,:,1,b)),neighbouring_matrix,0);
                        end
                    else
                        tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,1,b)),neighbouring_matrix,0);
                        end
                    end
                else
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),size(H0_Covariate_effect,3),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(:,:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,:,1,b)),neighbouring_matrix,0);
                        end
                    else
                        tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,b)),neighbouring_matrix,0);
                        end
                    end
                end
                save(fullfile(LIMO.dir,['H0' filesep 'tfce_' name]),'tfce_H0_score','-v7.3');
                clear H0_Covariate_effect tfce_H0_score
            end
        end
    end
end




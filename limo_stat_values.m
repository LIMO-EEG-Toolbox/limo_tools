function [M, mask, mytitle] = limo_stat_values(varargin)

% find corrected p values and mask from data under H0
%
% FORMAT [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO)
%
% INPUTS
%         Type = the type of plots (matters for titles only)
%         FileName = Name of the file selected
%         p = p value for thresholding
%         MCC = multiple comparisons option
%               1 none (useful for WLS)
%               2 clustering
%               3 TFCE
%               4 Max
%         LIMO = LIMO.mat
%
% OUTPUTS
%         M the (corrected) p values
%         mask the significant data
%         mytitle is the right title to use to make the figure for type 1 and 2
%
% see limo_display_results
%
% Cyril Pernet, Andrew Stewart, Marianne Latinus, Guilaume Rousselet 
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2019


Type      = varargin{1}; % type of plot
FileName  = varargin{2}; % Name of the file selected
p         = varargin{3}; % p value
MCC       = varargin{4}; % multiple comparison option
LIMO      = varargin{5}; % LIMO.mat

if isfield(LIMO,'Type')
    if strcmp(LIMO.Type,'Components') && MCC ~= 1, MCC = 4; end
else
    LIMO.Type = 'channels';
end

load (FileName);
M = []; mask =[]; mytitle=[];
c = clock; disp(' ');

if MCC ~= 1
    fprintf('limo_display_results %gh %gmin %gsec: computing statistical correction...\n',c(4),c(5),c(6));
    if MCC == 2
        disp('Ref for Clustering:')
        disp('Maris, E. & Oostenveld, R. 2007')
        disp('Nonparametric statistical testing of EEG- and MEG-data.')
        disp('Journal of Neuroscience Methods, 164, 177-190')
        disp(' ');
        disp('Pernet, C., Latinus, M., Nichols, T. & Rousselet, G.A. (2015).')
        disp('Cluster-based computational methods for mass univariate analyses')
        disp('of event-related brain potentials/fields: A simulation study.')
        disp('Journal of Neuroscience methods,250,83-95')
        
    elseif MCC == 3
        disp('Ref for TFCE:')
        disp('Pernet, C., Latinus, M., Nichols, T. & Rousselet, G.A. (2015).')
        disp('Cluster-based computational methods for mass univariate analyses')
        disp('of event-related brain potentials/fields: A simulation study.')
        disp('Journal of Neuroscience methods,250,83-95')
    end
else
    fprintf('limo_display_results %gh %gmin %gsec: making figure...\n',c(4),c(5),c(6));
end


%% Deal with each case of FileName


% -------------------------------
%% R2.mat (from 1st or 2nd level)
% -------------------------------

if strcmp(FileName,'R2.mat')
    
    M        = squeeze(R2(:,:,2)); % F values
    MCC_data = 'H0_R2.mat';
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if strcmp(LIMO.design.method,'WLS') 
            if all(isnan(squeeze(R2(:,size(R2,2)/2,3)))) && exist(['H0' filesep MCC_data],'file')
                load(['H0' filesep MCC_data]);
                H0_values = squeeze(H0_R2(:,:,2,:)); clear H0_R2;
                [mask,M]  = boot_threshold(H0_values,p,M);
                mytitle   = sprintf('R^2 : uncorrected threshold \n using bootstraped F values');
                R2(:,:,3) = M; save(FileName,'R2','-v7.3');
            elseif all(isnan(squeeze(R2(:,size(R2,2)/2,3)))) && ~exist(['H0' filesep MCC_data],'file')
                M         = squeeze(R2(:,:,3));
                mask      = ones(size(R2,1),size(R2,2));
                mytitle   = sprintf('unthresholded R^2 values \n no p-values available without bootstrap');
            else
                M         = squeeze(R2(:,:,3));
                mask      = squeeze(R2(:,:,3))< p;
                mytitle   = sprintf('R^2 : uncorrected threshold \n using bootstraped F values');
            end
        else
            mask    = squeeze(R2(:,:,3)) < p;
            M       = squeeze(R2(:,:,3)); % p values
            mytitle = sprintf('R^2 values: \n uncorrected threshold');
        end
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        if exist(['H0' filesep MCC_data],'file')
            try 
                H0_R2 = load(['H0' filesep MCC_data]); 
                bootM = squeeze(H0_R2.H0_R2(:,:,2,:)); % get all F values under H0
                bootP = squeeze(H0_R2.H0_R2(:,:,3,:)); % get all P values under H0
                if all(isnan(bootP(:)))
                    disp('WLS 1st pass - getting p values for null data')
                    Pboot = cell(1,size(bootP,3));
                    parfor boot = 1:size(bootP,3)
                        index = 1:size(bootP,3); index(boot) = [];
                        [~,Pboot{boot}]  = boot_threshold(bootM(:,:,index),p,bootM(:,:,boot));
                    end
                    bootP = reshape(cell2mat(Pboot),size(bootM)); clear Pboot;
                    H0_R2 = H0_R2.H0_R2; H0_R2(:,:,3,:) = bootP;
                    save(['H0' filesep MCC_data],'H0_R2','-v7.3'); 
                    disp('null p-values saved')
                end
                clear H0_R2
                
                if all(isnan(squeeze(R2(:,size(R2,2)/2,3))))
                    [~,M]  = boot_threshold(bootM,p,M);
                    R2(:,:,3) = M; save(FileName,'R2','-v7.3');
                end
                [mask,M] = limo_clustering(M,squeeze(R2(:,:,3)),bootM,bootP,LIMO,MCC,p); % mask and cluster p values
                Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
                if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
                mytitle = sprintf('R^2 values \n cluster correction (%g %s)', Nclust, Mclust);
            catch ME
                fprintf('someting went wrong load H0 data %s \n',ME.message)
                return
            end
        else
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        if exist(['H0' filsesep MCC_data],'file')
            try 
                H0_R2    = load(['H0' filesep MCC_data]); 
                bootM    = squeeze(H0_R2.H0_R2(:,:,2,:)); % get all F values under H0
                clear H0_R2; 
                [mask,M] = limo_max_correction(M,bootM,p);
                mytitle  = sprintf('R^2 values: \n correction by F max');
            catch ME
                fprintf('someting went wrong load H0 data %s \n',ME.message)
                return
            end
        else
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        if exist(['TFCE' filesep 'tfce_R2.mat'],'file')
            try 
                score    = load(['TFCE' filesep 'tfce_R2.mat']); 
                H0_score = load(['H0' filesep 'tfce_H0_R2.mat']); 
                [mask,M] = limo_max_correction(score.tfce_score,H0_score.tfce_H0_score,p);
                mytitle  = sprintf('R^2 values: \n correction using TFCE');
            catch ME
                fprintf('someting went wrong load H0 data %s \n',ME.message)
                return
            end
        else
            errordlg('no tfce bootstrap or tfce file was found to compute the tfce distribution','missing data')
        end
    end
    
    % ---------------------------------
    %% Condition_effect (from 1st level)
    % ---------------------------------
    
elseif strncmp(FileName,'Condition_effect',16)
    
    effect_nb = eval(FileName(18:end-4));
    M         = squeeze(Condition_effect(:,:,1)); % F values
    MCC_data  = sprintf('H0_Condition_effect_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if strcmp(LIMO.design.method,'WLS') 
            if exist(['H0' filesep MCC_data],'file')
                load(['H0' filesep MCC_data]);
                H0_values = squeeze(H0_Condition_effect(:,:,1,:)); clear H0_Condition;
                [mask,M]= boot_threshold(H0_values,p,M);
                mytitle = sprintf('Condition %g: uncorrected threshold \n using bootstraped F values',effect_nb);
            else
                mask = squeeze(Condition_effect(:,:,2)) ;
                M = squeeze(Condition_effect(:,:,2)); 
                mytitle = sprintf('unthresholded effect \n no p-values available without bootstrap');
            end
        else
            mask    = squeeze(Condition_effect(:,:,2)) < p;
            M       = squeeze(Condition_effect(:,:,2));
            mytitle = sprintf('Condition %g: F-values \n uncorrected threshold',effect_nb);
        end

        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Condition_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Condition_effect(:,:,2,:)); % get all P values under H0
            clear H0_Conditions;
            [mask,M] = limo_clustering(M,squeeze(Condition_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Condition %g: F-values \n cluster correction (%g %s)', effect_nb, Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Condition_effect(:,:,1,:)); % get all F values under H0
            clear H0_Condition_effect; [mask,M] = limo_max_correction(M,bootM,p);
            mytitle = sprintf('Condition %g: F-values \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            clear H0_Condition_effect; [mask,M] = limo_max_correction(tfce_score,tfce_H0_score,p);
            mytitle = sprintf('Condition %g: F-values \n correction using TFCE',effect_nb);
        catch ME
            errordlg('no tfce bootstrap or tfce file was found to compute the tfce distribution','missing data')
            return
        end
    end
    
    % ---------------------------------
    %% Covariate_effect (from 1st level)
    % ---------------------------------
    
elseif strncmp(FileName,'Covariate_effect',16)
    
    effect_nb = eval( FileName(18:end-4));
    M         = squeeze(Covariate_effect(:,:,1)); % F values
    MCC_data  = sprintf('H0_Covariate_effect_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if strcmp(LIMO.design.method,'WLS') 
            if exist(['H0' filesep MCC_data],'file')
                load(['H0' filesep MCC_data]);
                H0_values = squeeze(H0_Covariate_effect(:,:,1,:)); clear H0_Covariates;
                [mask,M]= boot_threshold(H0_values,p,M);
                mytitle = sprintf('Covariate %g: uncorrected threshold \n using bootstraped F values',effect_nb);
            else
                mask = squeeze(Covariate_effect(:,:,2));
                M = squeeze(Covariate_effect(:,:,2)); 
                mytitle = sprintf('unthresholded covariate effect \n no p-values available without bootstrap');
            end
        else
            mask    = squeeze(Covariate_effect(:,:,2)) < p;
            M       = squeeze(Covariate_effect(:,:,2)); % p values
            mytitle = sprintf('Covariate %g: F-values \n uncorrected threshold ',effect_nb);
        end
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Covariate_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Covariate_effect(:,:,2,:)); % get all P values under H0
            clear H0_Covariate_effect;
            [mask,M] = limo_clustering(M,squeeze(Covariate_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Covariate %g: F-values \n cluster correction (%g %s)', effect_nb, Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Covariate_effect(:,:,1,:)); % get all F values under H0
            [mask,M] = limo_max_correction(M,bootM,p);
            mytitle = sprintf('Covariate %g: F-values \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            [mask,M] = limo_max_correction(tfce_score,tfce_H0_score,p);
            mytitle = sprintf('Covariate %g: F-values \n correction using TFCE',effect_nb);
        catch ME
            errordlg('no tfce bootstrap file was found to compute the tfce distribution','missing data')
            return
        end
    end
    
    % ---------------------------------
    %% Interaction effect (from 1st level)
    % ---------------------------------
    
elseif strncmp(FileName,'Interaction_effect',18)
    
    effect_nb = eval( FileName(20:end-4));
    M         = squeeze(Interaction_effect(:,:,1)); % F values
    MCC_data  = sprintf('H0_Interaction_effect_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if strcmp(LIMO.design.method,'WLS') 
            if exist(['H0' filesep MCC_data],'file')
                load(['H0' filesep MCC_data]);
                H0_values = squeeze(H0_Interaction(:,:,1,:)); clear H0_Interaction;
                [mask,M]= boot_threshold(H0_values,p,M);
                mytitle = sprintf('Interaction %g: uncorrected threshold \n using bootstraped F values',effect_nb);
            else
                mask = squeeze(Interaction_effect(:,:,2));
                M = squeeze(Interaction_effect(:,:,2));
                mytitle = sprintf('unthresholded interaction \n no p-values available without bootstrap');
            end
        else
            mask    = squeeze(Interaction_effect(:,:,2)) < p;
            M       = squeeze(Interaction_effect(:,:,2));
            mytitle = sprintf('Interaction %g: F-values \n uncorrected threshold',effect_nb);
        end
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Interaction_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Interaction_effect(:,:,2,:)); % get all P values under H0
            clear H0_Interaction;
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Interaction %g: F-values \n cluster correction (%g %s)', effect_nb, Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Interaction_effect(:,:,1,:)); % get all F values under H0
            clear H0_Interaction_effect; [mask,M] = limo_max_correction(M,bootM,p);
            mytitle = sprintf('Interaction %g: F-values \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            clear H0_Interaction_effect; [mask,M] = limo_max_correction(tfce_score,tfce_H0_score,p);
            mytitle = sprintf('Interaction %g: F-values \n correction using TFCE',effect_nb);
        catch ME
            errordlg('no tfce bootstrap or tfce file was found to compute the tfce distribution','missing data')
            return
        end
    end
    
    % ---------------------------------
    %% Semi-partial coef (from 1st level)
    % ---------------------------------
    
elseif strncmp(FileName,'semi_partial_coef',17)
    
    effect_nb = eval(FileName(19:end-4));
    M         = squeeze(semi_partial_coef(:,:,2)); % F values
    MCC_data  = sprintf('H0_semi_partial_coef_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if strcmp(LIMO.design.method,'WLS') 
            if exist(['H0' filesep MCC_data],'file')
                load(['H0' filesep MCC_data]);
                H0_values = squeeze(H0_semi_partial_coef(:,:,2,:)); % get F values
                clear H0_semi_partial_coef;
                [mask,M]= boot_threshold(H0_values,p,M);
                mytitle = sprintf('Semi partial coef %g: uncorrected threshold \n using on bootstraped F values',effect_nb);
            else
                mask = squeeze(semi_partial_coef(:,:,3));
                M = squeeze(semi_partial_coef(:,:,3));
                mytitle = sprintf('unthresholded effect \n no p-values available without bootstrap');
            end
        else
            mask = squeeze(semi_partial_coef(:,:,3)) < p;
            M = squeeze(semi_partial_coef(:,:,3));
            mytitle = sprintf('Semi partial coef %g: F-values \n uncorrected threshold ',effect_nb);
        end
         
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_semi_partial_coef(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_semi_partial_coef(:,:,2,:)); % get all P values under H0
            clear H0_semi_partial_coef;
            [mask,M] = limo_clustering(M,squeeze(semi_partial_coef(:,:,3)),bootM,bootP,LIMO,MCC,p);
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Semi-parital coef %g: F-values \n cluster correction (%g %s)', effect_nb, Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_semi_partial_coef(:,:,1,:)); % get all F values under H0
            [mask,M] = limo_max_correction(M,bootM,p);
            mytitle = sprintf('Semi partial coef %g: F-values \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            clear H0_semi partial coef; [mask,M] = limo_max_correction(tfce_score,tfce_H0_score,p);
            mytitle = sprintf('Semi partial coef %g: F-values \n correction using TFCE',effect_nb);
        catch ME
            errordlg('no tfce bootstrap file was found to compute the tfce distribution','missing data')
            return
        end
    end
    
    
    % ---------------------------------
    %% Contrast T (from 1st or 2nd level)
    % --------------------
    
elseif strncmp(FileName,'con_',4)
    
    effect_nb = eval(FileName(5:end-4));
    M         = squeeze(con(:,:,4));
    MCC_data  = sprintf('H0_con_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if strcmp(LIMO.design.method,'WLS') 
            if exist(['H0' filesep MCC_data],'file')
                load(['H0' filesep MCC_data]);
                H0_values  = squeeze(boot_H0_con(:,:,2,:)); % T values under H0
                clear boot_H0_con
                [mask,M]= boot_threshold(H0_values,p,M);
            else
                mask = con(:,:,5);
                M = con(:,:,5);
                mytitle = sprintf('unthresholded contrast \n no p-values available without bootstrap');
            end
        else
            mask    = con(:,:,5) <= p;
            M       = con(:,:,5);
            mytitle = sprintf('Contrast %g: T-values \n uncorrected threshold',effect_nb);
        end
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        try cd H0; load(MCC_data); % dim electrode, frames, param/t/p, nboot
            bootT = H0_con(:,:,2,:); % get all t values under H0
            bootP = H0_con(:,:,3,:); % get all p values under H0
            clear H0_con
            [mask,M] = limo_clustering(M.^2,squeeze(con(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Contrast %g: T-values \n cluster correction (%g %s)', effect_nb, Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        
        % T max correction for multiple testing
        % --------------------------------------
    elseif MCC == 4 % Stat max
        try load(MCC_data);
            bootT = squeeze(H0_con(:,:,2,:)); % take all T values under H0
            [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % absolute T values
            mytitle = sprintf('Contrast %g: T-values \n correction by T max', effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            [mask,M] = limo_max_correction(tfce_score,tfce_H0_score,p);
            mytitle = sprintf('Contrast %g: T-values \n correction using TFCE', effect_nb);
        catch ME
            errordlg('no tfce bootstrap file was found to compute the tfce distribution','missing data')
            return
        end
    end
    
    
    % ---------------------------------
    %% Contrast F (1st or 2nd level)
    % ------------------------------------------
    
elseif strncmp(FileName,'ess_',4)
    
    start_at = max(strfind(FileName,'_'))+1;
    effect_nb = eval(FileName(start_at:end-4));
    if ~exist('ess','var')
        try
            ess = eval(['ess' num2str(effect_nb)]);
            clear(['ess' num2str(effect_nb)])
        catch
            ess = ess1; clear ess1
        end
    end
    M        = squeeze(ess(:,:,end-1));
    MCC_data = sprintf('H0_ess_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if strcmp(LIMO.design.method,'WLS') 
            if exist(['H0' filesep MCC_data],'file')
                load(['H0' filesep MCC_data]);
                % sort all F values
                H0_values = squeeze(boot_H0_ess(:,:,end-1,:));
                [mask,M]= boot_threshold(H0_values,p,M);
                mytitle = sprintf('Contrast F %g: threshold using bootstrapped F values',effect_nb);
            else
                mask = squeeze(ess(:,:,end));
                M = squeeze(ess(:,:,end));
                mytitle = sprintf('unthresholded contrast \n no p-values available without bootstrap');
            end
        else
            mask    = squeeze(ess(:,:,end)) < p;
            M       = squeeze(ess(:,:,end));
            mytitle = sprintf('Contrast %g: F-values \n uncorrected threshold', effect_nb);
        end
        
        % clustering correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd H0; load(MCC_data);
            bootF = squeeze(H0_ess(:,:,end-1,:));
            bootP = squeeze(H0_ess(:,:,end,:));
            clear H0_ess
            [mask,M] = limo_clustering(M,squeeze(ess(:,:,end)),bootF,bootP,LIMO,MCC,p);
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Contrast %g: F-values \n cluster correction (%g %s)', effect_nb, Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster correction','missing data')
            return
        end
        
        % T max correction for multiple testing
        % --------------------------------------
    elseif MCC == 4
        try cd H0; load(MCC_data);
            bootF  = squeeze(H0_ess(:,:,end-1,:)); clear H0_ess;
            [mask,M] = limo_max_correction(M,bootF,p); % absolute T values
            mytitle = sprintf('Contrast %g: F-values \n correction by F max', effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            [mask,M] = limo_max_correction(tfce_score,tfce_H0_score,p);
            mytitle = sprintf('Contrast %g: F-values \n correction using TFCE', effect_nb);
        catch ME
            errordlg('no tfce bootstrap file was found to compute the tfce distribution','missing data')
            return
        end
    end
    
    
    % ------------------------------------------
    %% One sample t-test
    % ------------------------------------------
    
elseif strncmp(FileName,'one_sample',10)
    
    try effect_nb = eval(FileName(28:end-4)); end

    if size(one_sample,1)>1
        M = squeeze(one_sample(:,:,4)); % T values
    else
        M = one_sample(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask    = one_sample(:,:,5) <= p;
        M       = squeeze(one_sample(:,:,5));
        mytitle = sprintf('One sample t-test: uncorrected threshold');
        
        
        % 2D cluster and 1D correction for multiple testing
        % ------------------------------------------
    elseif MCC == 2
        try load(MCC_data);
            bootT = squeeze(H0_one_sample(:,:,1,:)); % get all T values under H0
            bootP = squeeze(H0_one_sample(:,:,2,:)); % get all P values under H0
            if size(one_sample,1) == 1
                tmp = NaN(1,size(one_sample,2),size(H0_one_sample,4));
                tmp(1,:,:) = bootT; bootT = tmp;
                tmp(1,:,:) = bootP; bootP = tmp;
                clear tmp
            end
            clear H0_one_sample
          
            [mask,M] = limo_clustering(M.^2,squeeze(one_sample(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('One Sample t-test: t-values \n cluster correction (%g %s)', Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        
        % T max correction for multiple testing
        % -------------------------------------
    elseif MCC == 4 % Stat max
        try load(MCC_data)
            bootT = squeeze(H0_one_sample(:,:,1,:)); % get all T values under H0
            if size(one_sample,1) == 1
                tmp = NaN(1,size(one_sample,2),size(H0_one_sample,4));
                tmp(1,:,:) = bootT; bootT = tmp; clear tmp
            end
            [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('One Sample t-test: t values \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % Correction using TFCE
        % -------------------------------------
    elseif MCC == 3 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
            load(MCC_tfce_data)
            [mask,M] = limo_max_correction(tfce_one_sample, tfce_H0_one_sample,p);
            mytitle = sprintf('One Sample t-test: t values  \n correction using TFCE');
        catch ME
            errordlg('no tfce bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end
    
    % -----------------------------------------
    %% two samples t-test
    % -----------------------------------------
    
elseif strncmp(FileName,'two_samples',11)
    
    %effect_nb = eval(FileName(29:end-4));
    if size(two_samples,1)>1
        M = squeeze(two_samples(:,:,4)); % T values
    else
        M = two_samples(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask    = two_samples(:,:,5) <= p;
        M       = squeeze(two_samples(:,:,5));
        mytitle = sprintf('Two samples t-test: t values \n uncorrected threshold');
        
        
        % 2D cluster and 1D correction for multiple testing
        % ------------------------------------------
    elseif MCC == 2
        try load(MCC_data);
            bootT = squeeze(H0_two_samples(:,:,1,:)); % get all T values under H0
            bootP = squeeze(H0_two_samples(:,:,2,:)); % get all P values under H0
            if size(two_samples,1) == 1
                tmp = NaN(1,size(two_samples,2),size(H0_two_samples,4));
                tmp(1,:,:) = bootT; bootT = tmp;
                tmp(1,:,:) = bootP; bootP = tmp;
                clear tmp
            end
            [mask,M] = limo_clustering(M.^2,squeeze(two_samples(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Two Samples t-test: t-values \n cluster correction (%g %s)', Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        
        % T max correction for multiple testing
        % -------------------------------------
    elseif MCC == 4 % Stat max
        try load(MCC_data)
            bootT = squeeze(H0_two_samples(:,:,1,:)); % get all T values under H0
            if size(two_samples,1) == 1
                tmp = NaN(1,size(two_samples,2),size(H0_two_samples,4));
                tmp(1,:,:) = bootT; bootT = tmp; clear tmp
            end
            [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('Two Samples t-test: t values \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % Correction using TFCE
        % -------------------------------------
    elseif MCC == 3 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
            load(MCC_tfce_data)
            [mask,M] = limo_max_correction(tfce_two_samples, tfce_H0_two_samples,p);
            mytitle = sprintf('Two Samples t-test: t values \n correction using TFCE');
        catch ME
            errordlg('no tfce bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end
    
    % ---------------------
    %% paired t-test
    % --------------------
    
elseif strncmp(FileName,'paired_samples',14)
    
    %effect_nb = eval(FileName(32:end-4));
    if size(paired_samples,1)>1
        M = squeeze(paired_samples(:,:,4)); % T values
    else
        M = paired_samples(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = paired_samples(:,:,5) <= p;
        M = squeeze(paired_samples(:,:,5));
        mytitle = sprintf('Paired samples t-test: t values \n uncorrected threshold');
        
        
        % 2D cluster and 1D correction for multiple testing
        % ------------------------------------------
    elseif MCC == 2
        if size(M,1) == 1; MCC =3; end
        try load(MCC_data);
            bootT = squeeze(H0_paired_samples(:,:,1,:)); % get all T values under H0
            bootP = squeeze(H0_paired_samples(:,:,2,:)); % get all P values under H0
            if size(paired_samples,1) == 1
                tmp = NaN(1,size(paired_samples,2),size(H0_paired_samples,4));
                tmp(1,:,:) = bootT; bootT = tmp;
                tmp(1,:,:) = bootP; bootP = tmp;
                clear tmp
            end
            [mask,M] = limo_clustering(M.^2,squeeze(paired_samples(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('Paired t-test: t-values \n cluster correction (%g %s)', Nclust, Mclust);
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        
        % T max correction for multiple testing
        % -------------------------------------
    elseif MCC == 4 % Stat max
        try load(MCC_data)
            bootT = squeeze(H0_paired_samples(:,:,1,:)); % get all T values under H0
            if size(paired_samples,1) == 1
                tmp = NaN(1,size(paired_samples,2),size(H0_paired_samples,4));
                tmp(1,:,:) = bootT; bootT = tmp; clear tmp
            end
            [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('Paired Samples t-test: t values \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % Correction using TFCE
        % -------------------------------------
    elseif MCC == 3 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
            load(MCC_tfce_data)
            [mask,M] = limo_max_correction(tfce_paired_samples, tfce_H0_paired_samples,p);
            mytitle = sprintf('Paired Samples t-test: t values \n correction using TFCE');
        catch ME
            errordlg('no tfce bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end
    
    % ------------------------
    %% Repeated measure ANOVA
    % ------------------------
    
elseif strncmp(FileName,'Rep_ANOVA',9)
    
    % all files have dim electrode x time frames x F/p
    if strncmp(FileName,'Rep_ANOVA_Interaction',21)
        M    = Rep_ANOVA_Interaction_with_gp(:,:,1); % get the F values
        PVAL = Rep_ANOVA_Interaction_with_gp(:,:,2);
    elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
        M    = Rep_ANOVA_Gp_effect(:,:,1); % get the F values
        PVAL = Rep_ANOVA_Gp_effect(:,:,2);
    elseif strncmp(FileName,'Rep_ANOVA',9)
        M    = Rep_ANOVA(:,:,1); % get the F values
        PVAL = Rep_ANOVA(:,:,2);
    end
    
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = PVAL <= p;
        M = PVAL;
        if strncmp(FileName,'Rep_ANOVA_Interaction',21)
            mytitle = sprintf('Rep ANOVA Interaction: F-values \n uncorrected threshold');
        elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
            mytitle = sprintf('Rep ANOVA Gp effect: F-values \n uncorrected threshold');
        elseif strncmp(FileName,'Rep_ANOVA',9)
            mytitle = sprintf('Rep ANOVA: F-values \n uncorrected threshold');
        end
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try load(MCC_data);
            if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                bootT = squeeze(H0_Rep_ANOVA_Interaction_with_gp(:,:,1,:));
                bootP = squeeze(H0_Rep_ANOVA_Interaction_with_gp(:,:,2,:));
                if size(Rep_ANOVA_Interaction_with_gp,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA_Interaction_with_gp,2),size(H0_Rep_ANOVA_Interaction_with_gp,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
                clear H0_Rep_ANOVA_Interaction_with_gp
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                bootT = squeeze(H0_Rep_ANOVA_Gp_effect(:,:,1,:));
                bootP = squeeze(H0_Rep_ANOVA_Gp_effect(:,:,2,:));
                if size(Rep_ANOVA_Gp_effect,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA_Gp_effect,2),size(H0_Rep_ANOVA_Gp_effect,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
                clear H0_Rep_ANOVA_Gp_effect
            elseif strncmp(FileName,'Rep_ANOVA',9)
                bootT = squeeze(H0_Rep_ANOVA(:,:,1,:)); % get all F values under H0
                bootP = squeeze(H0_Rep_ANOVA(:,:,2,:)); % get all P values under H0
                if size(Rep_ANOVA,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA,2),size(H0_Rep_ANOVA,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
                clear H0_Rep_ANOVA
            end
            
            [mask,M] = limo_clustering(M,PVAL,bootT,bootP,LIMO,MCC,p);
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
             if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                 mytitle = sprintf('Rep ANOVA Interaction: F-values \n cluster correction (%g %s)', Nclust, Mclust);
             elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                 mytitle = sprintf('Rep ANOVA Gp effect: F-values \n cluster correction (%g %s)', Nclust, Mclust);
             elseif strncmp(FileName,'Rep_ANOVA',9)
                 mytitle = sprintf('Rep ANOVA: F-values \n cluster correction (%g %s)', Nclust, Mclust);
             end
            
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        
        % T max correction for multiple testing
        % -------------------------------------
    elseif MCC == 4 % Stat max
        try load(MCC_data);
            if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                bootT = H0_Rep_ANOVA_Interaction_with_gp(:,:,1,:);
                if size(Rep_ANOVA_Interaction_with_gp,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA_Interaction_with_gp,2),size(H0_Rep_ANOVA_Interaction_with_gp,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                bootT = H0_Rep_ANOVA_Gp_effect(:,:,1,:);
                if size(Rep_ANOVA_Gp_effect,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA_Gp_effect,2),size(H0_Rep_ANOVA_Gp_effect,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
            elseif strncmp(FileName,'Rep_ANOVA',9)
                bootT = H0_Rep_ANOVA(:,:,1,:); % get all F values under H0
                if size(Rep_ANOVA,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA,2),size(H0_Rep_ANOVA,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
            end
            
            [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                mytitle = sprintf('Rep ANOVA Interaction: \n correction by T max');
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                mytitle = sprintf('Rep ANOVA Gp effect: \n correction by T max');
            elseif strncmp(FileName,'Rep_ANOVA',9)
                mytitle = sprintf('Rep ANOVA: \n correction by T max');
            end
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % Correction using TFCE
        % -------------------------------------
    elseif MCC == 3 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
            load(MCC_tfce_data)
            
            if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                [mask,M] = limo_max_correction(tfce_Rep_ANOVA_Interaction_with_gp, tfce_H0_Rep_ANOVA_Interaction_with_gp,p);
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                [mask,M] = limo_max_correction(tfce_Rep_ANOVA_Gp_effect, tfce_H0_Rep_ANOVA_Gp_effect,p);
            elseif strncmp(FileName,'Rep_ANOVA',9)
                [mask,M] = limo_max_correction(tfce_Rep_ANOVA, tfce_H0_Rep_ANOVA,p);
            end
            if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                mytitle = sprintf('Rep ANOVA Interaction: \n correction using TFCE');
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                mytitle = sprintf('Rep ANOVA Gp effect: \n correction using TFCE');
            elseif strncmp(FileName,'Rep_ANOVA',9)
                mytitle = sprintf('Rep ANOVA: \n correction using TFCE');
            end
        catch ME
            errordlg('no tfce bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end
    
    
    % -----------------------
    %% Lateralization maps
    % -----------------------
    
elseif strncmp(FileName,'LI_Map',6)
    
    M = squeeze(LI_Map(:,:,4)); % T values
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = LI_Map(:,:,5) <= p;
        mytitle = sprintf('LI Map: one sample T values \n threshold based on theoretical p values');
        
        
        % 2D cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        MCC_data = sprintf('boot_%s',FileName);
        try load(MCC_data);
            bootT = squeeze(boot_LI_Map(:,:,2,:)); % get all T values under H0
            bootP = squeeze(boot_LI_Map(:,:,3,:)); % get all P values under H0
            [mask,M] = limo_clustering(M.^2,squeeze(LI(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
             Nclust = unique(M(~isnan(M))); Nclust = length(Nclust) ;
             if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('LI Map: T-values \n cluster correction (%g %s)', Nclust, Mclust);
            
        catch ME
            errordlg('no bootstrap file was found to compute the cluster correction','missing data')
            return
        end
        
        
        % T max correction for multiple testing
        % --------------------------------------
    elseif MCC == 4
        MCC_data = sprintf('boot_%s',FileName);
        try load(MCC_data);
            bootT  = squeeze(boot_LI_Map(:,:,2,:)); % take all T values under H0
            [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('LI Map: One sample T values \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % Correction using TFCE
        % -------------------------------------
    elseif MCC == 3 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
            load(MCC_tfce_data)
            [mask,M] = limo_max_correction(tfce_LI, tfce_H0_LI,p);
            mytitle = sprintf('LI Map: One Sample t-test \n correction using TFCE');
        catch ME
            errordlg('no tfce bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end
    
    %% any other files
else
    errordlg2('unidentified FileName - no thresholding done');
end


function [mask,M]= boot_threshold(H0_values,p,M)
% computes the cell-wise p value using bootrap

sorted_values = sort(H0_values,3); clear H0_values;
if all(sorted_values(:))>=0 % i.e. F values
    U = round((1-p)*size(sorted_values,3));
    mask = (M >= sorted_values(:,:,U));
    for row = 1:size(M,1)
        for column = 1:size(M,2)
            tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
        end
    end
    M = 1- (tmp ./ size(sorted_values,3)) ; % p values
else % i.e. T values
    low = round(p*size(sorted_values,3)/2);
    high = size(sorted_values,3) - low;
    mask = (M <= sorted_values(:,:,low))+(M >= sorted_values(:,:,high));
    for row = 1:size(M,1)
        for column = 1:size(M,2)
            tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
        end
    end
    M = min((tmp ./ size(sorted_values,3)), 1- (tmp ./ size(sorted_values,3))) ; % p values
end

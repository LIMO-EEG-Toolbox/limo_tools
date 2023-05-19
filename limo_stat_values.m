function [M, mask, mytitle] = limo_stat_values(varargin)

% find corrected p values and mask from data under H0
%
% FORMAT [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO)
%
% INPUTS
%         FileName = Name of the file selected
%         p        = p value for thresholding
%         MCC      = multiple comparisons option
%                    1 none
%                    2 clustering
%                    3 TFCE
%                    4 Max
%         LIMO     = LIMO.mat structure
%
% OUTPUTS
%         M        = the (corrected) p values
%         mask     = the significant data
%         mytitle  = title that include the filename, effect and method used
%
% see limo_display_results
%
% Cyril Pernet, Andrew Stewart, Marianne Latinus, Guilaume Rousselet
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2020

root = fileparts(which('limo_eeg'));
pathCell = regexp(path, pathsep, 'split');
onPath = all([sum(strcmp([root filesep 'help'],pathCell))~=0,...
    sum(strcmp([root filesep 'limo_cluster_functions'],pathCell))~=0,...
    sum(strcmp([root filesep 'external' filesep 'psom'],pathCell))~=0,...
    sum(strcmp([root filesep 'deprecated'], pathCell))~=0]);
if onPath == 0
    addpath([root filesep 'limo_cluster_functions'])
    addpath([root filesep 'external'])
    addpath([root filesep 'external' filesep 'psom'])
    addpath([root filesep 'help'])
    addpath([root filesep 'deprecated'])
end

FileName  = varargin{1}; % Name of the file selected
p         = varargin{2}; % p value
MCC       = varargin{3}; % multiple comparison option
LIMO      = varargin{4}; % LIMO.mat

% check the appropriate method is used
% -----------------------------------
if isfield(LIMO,'Type')
    if strcmp(LIMO.Type,'Components') && MCC ~= 1
        MCC = 4; % for ICA only max stat since we can't cluster them based a topography
        disp('only maximum statistics can be used for ICA')
    end
end

% check that a neighbouring matrix is there for clustering
% -------------------------------------------------------
if MCC == 2
    limo_check_neighbourghs(LIMO)
end

% load data and set outputs to empty
% ----------------------------------
if exist(FileName,'file')
    matfile = load(FileName);
else
    matfile = load(fullfile(LIMO.dir,FileName));
end
M       = []; 
mask    = []; 
mytitle = [];
disp(' ');

% disp some references for this
% -----------------------------
if MCC ~= 1
    if MCC == 2
        disp('Ref for Clustering & Bootstrap:')
        disp('Maris, E. & Oostenveld, R. 2007')
        disp('Nonparametric statistical testing of EEG- and MEG-data.')
        disp('Journal of Neuroscience Methods, 164, 177-190')
        disp(' ');
        disp('Pernet, C., Latinus, M., Nichols, T. & Rousselet, G.A. (2015).')
        disp('Cluster-based computational methods for mass univariate analyses')
        disp('of event-related brain potentials/fields: A simulation study.')
        disp('Journal of Neuroscience methods, 250, 83-95')
        disp(' ');
    elseif MCC == 3
        disp('Ref for TFCE:')
        disp('Pernet, C., Latinus, M., Nichols, T. & Rousselet, G.A. (2015).')
        disp('Cluster-based computational methods for mass univariate analyses')
        disp('of event-related brain potentials/fields: A simulation study.')
        disp('Journal of Neuroscience methods, 250, 83-95')
        disp(' ');
    end
end

if MCC ~= 1
    fprintf('computing corrected statistics at %s...\n',datetime('now','Format','hh:mm:ss'));
end


%% Deal with each case of FileName

% -------------------------------
%% GLM (from 1st or 2nd level) also robust regresion, robust ANOVA
% ------------------------------------------------------------------------
if strcmpi(LIMO.Analysis,'Time-Frequency')
    if strcmp(FileName,'R2.mat')
        M         = squeeze(matfile.R2(:,:,:,2)); % F values
        Pval      = squeeze(matfile.R2(:,:,:,3)); % P values
        MCC_data  = 'H0_R2.mat';
        titlename = 'R^2 Coef';
    elseif strncmp(FileName,'Condition_effect',16)
        effect_nb = eval(FileName(18:end-4));
        M         = squeeze(matfile.Condition_effect(:,:,:,1));
        Pval      = squeeze(matfile.Condition_effect(:,:,:,2));
        MCC_data  = sprintf('H0_Condition_effect_%g.mat',effect_nb);
        titlename = sprintf('Condition effect %g  F values',effect_nb);
    elseif strncmp(FileName,'Covariate_effect',16)
        effect_nb = eval(FileName(18:end-4));
        M         = squeeze(matfile.Covariate_effect(:,:,:,1));
        Pval      = squeeze(matfile.Covariate_effect(:,:,:,2));
        MCC_data  = sprintf('H0_Covariate_effect_%g.mat',effect_nb);
        titlename = sprintf('Covariate effect %g  F values',effect_nb);
    elseif strncmp(FileName,'Interaction_effect',18)
        effect_nb = eval(FileName(20:end-4));
        M         = squeeze(matfile.Interaction_effect(:,:,:,1));
        Pval      = squeeze(matfile.Interaction_effect(:,:,:,2));
        MCC_data  = sprintf('H0_Interaction_effect_%g.mat',effect_nb);
        titlename = sprintf('Interaction effect %g  F values',effect_nb);
    elseif strncmp(FileName,'semi_partial_coef',17)
        effect_nb = eval(FileName(19:end-4));
        M         = squeeze(matfile.semi_partial_coef(:,:,:,2));
        Pval      = squeeze(matfile.semi_partial_coef(:,:,:,3));
        MCC_data  = sprintf('H0_semi_partial_coef_%g.mat',effect_nb);
        titlename = sprintf('Semi Partial Coef %g',effect_nb);
    elseif strncmp(FileName,'con_',4)
        effect_nb = eval(FileName(5:end-4));
        M         = squeeze(matfile.con(:,:,:,4));
        Pval      = squeeze(matfile.con(:,:,:,5));
        MCC_data  = sprintf('H0_con_%g.mat',effect_nb);
        titlename = sprintf('Contrast %g T values',effect_nb);
    elseif contains(FileName,'ttest') || contains(FileName,'LI_Map')
        matfile   = matfile.(cell2mat(fieldnames(matfile)));
        M         = matfile(:,:,:,4); % T values
        Pval      = matfile(:,:,:,5);
        MCC_data  = sprintf('H0_%s', FileName);
        name      = FileName(1:strfind(FileName,'ttest')+4);
        name(strfind(name,'_')) = ' ';
        titlename = sprintf('%s t-test T values',name);
    elseif strncmp(FileName,'ess_',4)
        effect_nb = eval(FileName(5:end-4));
        M         = squeeze(matfile.ess(:,:,:,end-1));
        Pval      = squeeze(matfile.ess(:,:,:,end));
        MCC_data  = sprintf('H0_ess_%g.mat',effect_nb);
        titlename = sprintf('Contrast %g F values',effect_nb);
    end
    
else  % same with one less dimention
    if strcmp(FileName,'R2.mat')
        M         = squeeze(matfile.R2(:,:,2)); % F values
        Pval      = squeeze(matfile.R2(:,:,3)); % P values
        MCC_data  = 'H0_R2.mat';
        titlename = 'R^2 Coef';
    elseif strncmp(FileName,'Condition_effect',16)
        effect_nb = eval(FileName(18:end-4));
        M         = squeeze(matfile.Condition_effect(:,:,1));
        Pval      = squeeze(matfile.Condition_effect(:,:,2));
        MCC_data  = sprintf('H0_Condition_effect_%g.mat',effect_nb);
        titlename = sprintf('Condition effect %g  F values',effect_nb);
    elseif strncmp(FileName,'Covariate_effect',16)
        effect_nb = eval(FileName(18:end-4));
        M         = squeeze(matfile.Covariate_effect(:,:,1));
        Pval      = squeeze(matfile.Covariate_effect(:,:,2));
        MCC_data  = sprintf('H0_Covariate_effect_%g.mat',effect_nb);
        titlename = sprintf('Covariate effect %g  F values',effect_nb);
    elseif strncmp(FileName,'Interaction_effect',18)
        effect_nb = eval(FileName(20:end-4));
        M         = squeeze(matfile.Interaction_effect(:,:,1));
        Pval      = squeeze(matfile.Interaction_effect(:,:,2));
        MCC_data  = sprintf('H0_Interaction_effect_%g.mat',effect_nb);
        titlename = sprintf('Interaction effect %g  F values',effect_nb);
    elseif strncmp(FileName,'semi_partial_coef',17)
        effect_nb = eval(FileName(19:end-4));
        M         = squeeze(matfile.semi_partial_coef(:,:,2));
        Pval      = squeeze(matfile.semi_partial_coef(:,:,3));
        MCC_data  = sprintf('H0_semi_partial_coef_%g.mat',effect_nb);
        titlename = sprintf('Semi Partial Coef %g',effect_nb);
    elseif strncmp(FileName,'con_',4)
        effect_nb = eval(FileName(5:end-4));
        M         = squeeze(matfile.con(:,:,4));
        Pval      = squeeze(matfile.con(:,:,5));
        MCC_data  = sprintf('H0_con_%g.mat',effect_nb);
        titlename = sprintf('Contrast %g T values',effect_nb);
    elseif contains(FileName,'ttest') || contains(FileName,'LI_Map')
        matfile   = matfile.(cell2mat(fieldnames(matfile)));
        M         = matfile(:,:,4); % T values
        Pval      = matfile(:,:,5);
        MCC_data  = sprintf('H0_%s', FileName);
        name      = FileName(1:strfind(FileName,'ttest')+4);
        name(strfind(name,'_')) = ' ';
        titlename = sprintf('%s T values',name);
    elseif strncmp(FileName,'ess_',4)
        effect_nb = eval(FileName(max(strfind(FileName,'_'))+1:end-4));
        M         = squeeze(matfile.ess(:,:,end-1));
        Pval      = squeeze(matfile.ess(:,:,end));
        MCC_data  = sprintf('H0_ess_%g.mat',effect_nb);
        titlename = sprintf('Contrast %g F values',effect_nb);
    end
end

% no correction for multiple testing
% -----------------------------------
if ~isempty(M) && MCC == 1
    M       = Pval;
    mask    = Pval <= p;
    mytitle = sprintf('%s: uncorrected threshold',titlename);
    
    % cluster correction for multiple testing
    % ---------------------------------------
elseif ~isempty(M) && MCC == 2
    if exist(['H0' filesep MCC_data],'file')
        try
            H0_data = load(['H0' filesep MCC_data]);
            H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                if contains(FileName,'R2') || contains(FileName,'semi_partial')
                    bootM = squeeze(H0_data(:,:,:,2,:)); % get all F values under H0
                    bootP = squeeze(H0_data(:,:,:,3,:)); % get all P values under H0
               else
                    bootM = squeeze(H0_data(:,:,:,1,:));
                    bootP = squeeze(H0_data(:,:,:,2,:));
                end
                
                if size(M,1) == 1
                    tmp = NaN(1,size(M,2),size(M,3),size(bootM,3));
                    tmp(1,:,:,:) = bootM; bootM = tmp;
                    tmp(1,:,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
            else
                if contains(FileName,'R2') || contains(FileName,'semi_partial')
                    bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
                    bootP = squeeze(H0_data(:,:,3,:)); % get all P values under H0
                else
                    bootM = squeeze(H0_data(:,:,1,:));
                    bootP = squeeze(H0_data(:,:,2,:));
                end
                
                if size(M,1) == 1
                    tmp = NaN(1,size(M,2),size(bootM,2));
                    tmp(1,:,:,:) = bootM; bootM = tmp;
                    tmp(1,:,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
            end
            
            % finally get cluster mask and corrected p-values
            if contains(FileName,'ttest') || contains(FileName,'LI_Map')
                [mask,M] = limo_clustering(M.^2,Pval,bootM.^2,bootP,LIMO,MCC,p); % mask and cluster p values
            else
                [mask,M] = limo_clustering(M,Pval,bootM,bootP,LIMO,MCC,p); % mask and cluster p values
            end
                
            Nclust   = unique(mask); 
            Nclust   = length(Nclust)-1; % mask = mask>0;
            if Nclust <= 1
                Mclust = 'cluster';
            else
                Mclust = 'clusters';
            end
            mytitle = sprintf('%s cluster correction (%g %s)', titlename, Nclust, Mclust);
        catch ME
            errordlg(sprintf('error log: %s \n',ME.message),'cluster correction failure')
            return
        end
    else
        errordlg(['H0' filesep MCC_data ' not found'],'cluster correction failure')
        return
    end
    
    % correction using the max
    % --------------------------
elseif ~isempty(M) && MCC == 4 % Stat max
    if exist(['H0' filesep MCC_data],'file')
        try
            H0_data = load(['H0' filesep MCC_data]);
            H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                if contains(FileName,'R2') || contains(FileName,'semi_partial')
                    bootM = squeeze(H0_data(:,:,:,2,:)); % get all F values under H0
                else
                    bootM = squeeze(H0_data(:,:,:,1,:));
                end
            else
                if contains(FileName,'R2') || contains(FileName,'semi_partial')
                    bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
                else
                    bootM = squeeze(H0_data(:,:,1,:));
                end
            end
            clear H0_data;
            [mask,M] = limo_max_correction(abs(M),abs(bootM),p);
            mytitle  = sprintf('%s: correction by max',titlename);
        catch ME
            errordlg(sprintf('error log: %s \n',ME.message),'max correction failure')
            return
        end
    else
        errordlg(['H0' filesep MCC_data ' not found'],'max correction failure')
    end
    
    % correction using TFCE
    % --------------------------
elseif ~isempty(M) && MCC == 3 % Stat max
    if exist(fullfile(LIMO.dir,['tfce' filesep 'tfce_' FileName]),'file')
        try
            score    = load(fullfile(LIMO.dir,['tfce' filesep 'tfce_' FileName]));
            score    = score.(cell2mat(fieldnames(score)));
            H0_score = load(fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' FileName]));
            H0_score = H0_score.(cell2mat(fieldnames(H0_score)));
            [mask,M] = limo_max_correction(score,H0_score,p);
            mytitle  = sprintf('%s: correction using TFCE',titlename);
        catch ME
            errordlg(sprintf('error log: %s \n',ME.message),'tfce correction failure')
            return
        end
    else
        errordlg('no tfce tfce file was found','missing data')
    end
end

% ------------------------
%% Repeated measures ANOVA
% ------------------------

if contains(FileName,'Rep_ANOVA')
    
    % all files have dim electrode x [freq/time] frames x F/p
    if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
        M    = matfile.(cell2mat(fieldnames(matfile)))(:,:,:,1);
        PVAL = matfile.(cell2mat(fieldnames(matfile)))(:,:,:,2);
    else
        M    = matfile.(cell2mat(fieldnames(matfile)))(:,:,1);
        PVAL = matfile.(cell2mat(fieldnames(matfile)))(:,:,2);
    end
    MCC_data = fullfile(LIMO.dir,['H0' filesep 'H0_' FileName]);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = PVAL <= p;
        M    = PVAL;
        if contains(FileName,'Rep_ANOVA_Interaction')
            mytitle = sprintf('Interaction F-values uncorrected threshold');
        elseif contains(FileName,'Rep_ANOVA_Gp_effect')
            mytitle = sprintf('Gp effect F-values uncorrected threshold');
        elseif contains(FileName,'Rep_ANOVA_Main')
            mytitle = sprintf('Main Effect F-values uncorrected threshold');
        end
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        if exist(MCC_data,'file')
            try
                H0_data = load(MCC_data);
                H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    bootT = squeeze(H0_data(:,:,:,1,:));
                    bootP = squeeze(H0_data(:,:,:,2,:));
                    if size(M,1) == 1
                        tmp = NaN(1,size(M,2),size(M,3),size(bootT,4));
                        tmp(1,:,:,:) = bootT; bootT = tmp;
                        tmp(1,:,:,:) = bootP; bootP = tmp;
                        clear tmp
                    end
                else
                    bootT = squeeze(H0_data(:,:,1,:));
                    bootP = squeeze(H0_data(:,:,2,:));
                    if size(M,1) == 1
                        tmp = NaN(1,size(M,2),size(bootT,2));
                        tmp(1,:,:) = bootT; bootT = tmp;
                        tmp(1,:,:) = bootP; bootP = tmp;
                        clear tmp
                    end
                end
                
                if size(M,1) == 1
                    [mask,M] = limo_clustering(M,PVAL,bootT,bootP,LIMO,3,p); % temporal clustering
                else
                    [mask,M] = limo_clustering(M,PVAL,bootT,bootP,LIMO,2,p); % spatial-temporal clustering
                end
                Nclust = unique(mask); Nclust = length(Nclust)-1; % mask = mask>0;
                if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
                if contains(FileName,'Rep_ANOVA_Interaction')
                    mytitle = sprintf('Interaction F-values cluster correction (%g %s)', Nclust, Mclust);
                elseif contains(FileName,'Rep_ANOVA_Gp_effect')
                    mytitle = sprintf('Gp effect F-values cluster correction (%g %s)', Nclust, Mclust);
                elseif contains(FileName,'Rep_ANOVA_Main')
                    mytitle = sprintf('Main effect F-values cluster correction (%g %s)', Nclust, Mclust);
                end
                
            catch ME
                errordlg(sprintf('error log: %s \n',ME.message),'cluster correction failure')
                return
            end
        else
            errordlg(['H0' filesep MCC_data ' not found'],'cluster correction failure')
            return
        end
        
        
        % T max correction for multiple testing
        % -------------------------------------
    elseif MCC == 4 % Stat max
        if exist(MCC_data,'file')
            try
                H0_data = load(MCC_data);
                H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    bootT = squeeze(H0_data(:,:,:,1,:));
                    if size(M,1) == 1
                        tmp = NaN(1,size(M,2),size(M,3),size(bootT,4));
                        tmp(1,:,:,:) = bootT; bootT = tmp;
                        clear tmp
                    end
                else
                    bootT = squeeze(H0_data(:,:,1,:));
                    if size(M,1) == 1
                        tmp = NaN(1,size(M,2),size(bootT,3));
                        tmp(1,:,:) = bootT; bootT = tmp;
                        clear tmp
                    elseif size(M,2) == 1 %% for Weight bias testing
                        tmp = NaN(size(M,1),1,size(bootT,2));
                        tmp(:,1,:) = bootT; bootT = tmp;
                        clear tmp
                    end
                end
                
                [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    mytitle = sprintf('Interaction correction by T max');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    mytitle = sprintf('Gp effect correction by T max');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    mytitle = sprintf('Main Effect correction by T max');
                end
            catch ME
                errordlg(sprintf('error log: %s \n',ME.message),'max correction failure')
                return
            end
            
        else
            errordlg(['H0' filesep MCC_data ' not found'],'max correction failure')
            return
        end
        
        % Correction using TFCE
        % -------------------------------------
    elseif MCC == 3 % Stat tfce
        tfce_data    = sprintf('tfce%stfce_%s',filesep, FileName);
        H0_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        if exist(tfce_data,'file') && exist(H0_tfce_data,'file')
            try
                tfce_data    = load(tfce_data);     
                tfce_data    = tfce_data.(cell2mat(fieldnames(tfce_data)));
                H0_tfce_data = load(H0_tfce_data);
                H0_tfce_data = H0_tfce_data.(cell2mat(fieldnames(H0_tfce_data)));
                [mask,M]     = limo_max_correction(tfce_data, H0_tfce_data,p);
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    mytitle = sprintf('Interaction correction using TFCE');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    mytitle = sprintf('Gp effect correction using TFCE');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    mytitle = sprintf('Main Effect correction using TFCE');
                end
            catch ME
                errordlg(sprintf('error log: %s \n',ME.message),'tfce correction failure')
                return
            end
        else
            errordlg('no tfce file or bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end
end


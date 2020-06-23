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
matfile = load(fullfile(LIMO.dir,FileName)); 
M       = []; 
mask    = []; 
mytitle = [];
c       = clock; 
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
fprintf('limo_display_results %gh %gmin %gsec: making figure...\n',c(4),c(5),c(6));


%% Deal with each case of FileName

% -------------------------------
%% GLM (from 1st or 2nd level)
% -------------------------------
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
    elseif strncmp(FileName,'ess_',4)
        effect_nb = eval(FileName(5:end-4));
        M         = squeeze(matfile.ess(:,:,:,end-1));
        Pval      = squeeze(matfile.ess(:,:,:,end));
        MCC_data  = sprintf('H0_ess_%g.mat',effect_nb);
        titlename = sprintf('Contrast %g F values',effect_nb);
    end
else
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
    elseif strncmp(FileName,'ess_',4)
        effect_nb = eval(FileName(5:end-4));
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
                if strcmp(FileName,'R2.mat')
                    bootM = squeeze(H0_data(:,:,:,2,:)); % get all F values under H0
                    bootP = squeeze(H0_data(:,:,:,3,:)); % get all P values under H0
                else
                    bootM = squeeze(H0_data(:,:,:,1,:));
                    bootP = squeeze(H0_data(:,:,:,2,:));
                end
            else
                if strcmp(FileName,'R2.mat')
                    bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
                    bootP = squeeze(H0_data(:,:,3,:)); % get all P values under H0
                else
                    bootM = squeeze(H0_data(:,:,1,:));
                    bootP = squeeze(H0_data(:,:,2,:));
                end
            end
            
            % finally get cluster mask and corrected p-values
            [mask,M] = limo_clustering(M,Pval,bootM,bootP,LIMO,MCC,p); % mask and cluster p values
            Nclust = unique(mask); Nclust = length(Nclust)-1; % mask = mask>0;
            if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
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
                if strcmp(FileName,'R2.mat')
                    bootM = squeeze(H0_data(:,:,:,2,:)); % get all F values under H0
                else
                    bootM = squeeze(H0_data(:,:,:,1,:));
                end
            else
                if strcmp(FileName,'R2.mat')
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
            H0_score = load(fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' FileName]));
            [mask,M] = limo_max_correction(score.tfce_score,H0_score.tfce_H0_score,p);
            mytitle  = sprintf('%s: correction using TFCE',titlename);
        catch ME
            errordlg(sprintf('error log: %s \n',ME.message),'tfce correction failure')
            return
        end
    else
        errordlg('no tfce tfce file was found','missing data')
    end
end

% ------------------------------------------
%% One sample t-test
% ------------------------------------------

if strncmp(FileName,'one_sample',10)
       
    if size(matfile.one_sample,1)>1
        M = squeeze(matfile.one_sample(:,:,4)); % T values
    else
        M = matfile.one_sample(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask    = matfile.one_sample(:,:,5) <= p;
        M       = squeeze(matfile.one_sample(:,:,5));
        mytitle = sprintf('One sample t-test: uncorrected threshold');
        
        % 2D cluster and 1D correction for multiple testing
        % ------------------------------------------
    elseif MCC == 2
        if exist(MCC_data,'file')
            try
                H0_one_sample = load(MCC_data);
                H0_one_sample = H0_one_sample.H0_one_sample;
                bootT         = squeeze(H0_one_sample(:,:,1,:)); % get all T values under H0
                bootP         = squeeze(H0_one_sample(:,:,2,:)); % get all P values under H0
                if size(matfile.one_sample,1) == 1
                    tmp = NaN(1,size(matfile.one_sample,2),size(H0_one_sample,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
                clear H0_one_sample
                
                if size(M,1) == 1
                    [mask,M] = limo_clustering(M.^2,squeeze(matfile.one_sample(:,:,5)),bootT.^2,bootP,LIMO,3,p);
                else
                    [mask,M] = limo_clustering(M.^2,squeeze(matfile.one_sample(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
                end
                Nclust = unique(mask); Nclust = length(Nclust)-1; % mask = mask>0;
                if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
                mytitle = sprintf('One Sample t-values cluster correction (%g %s)', Nclust, Mclust);
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
                H0_one_sample = load(MCC_data);
                H0_one_sample = H0_one_sample.H0_one_sample;
                bootT         = squeeze(H0_one_sample(:,:,1,:)); % get all T values under H0
                if size(matfile.one_sample,1) == 1
                    tmp = NaN(1,size(matfile.one_sample,2),size(H0_one_sample,4));
                    tmp(1,:,:) = bootT; bootT = tmp; clear tmp
                end
                clear H0_one_sample
                [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
                mytitle = sprintf('One Sample t-values correction by T max');
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
                H0_tfce_data = load(H0_tfce_data);
                [mask,M]     = limo_max_correction(tfce_data.tfce_one_sample, H0_tfce_data.tfce_H0_one_sample,p);
                mytitle      = sprintf('One Sample t-values correction using TFCE');
            catch ME
                errordlg(sprintf('error log: %s \n',ME.message),'tfce correction failure')
                return
            end
        else
            errordlg('no tfce file or bootstrap file was found to compute the max distribution','tfce correction failure')
            return
        end
    end
end

% -----------------------------------------
%% two samples t-test
% -----------------------------------------

if strncmp(FileName,'two_samples',11)
       
    if size(matfile.two_samples,1)>1
        M = squeeze(matfile.two_samples(:,:,4)); % T values
    else
        M = matfile.two_samples(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask    = matfile.two_samples(:,:,5) <= p;
        M       = squeeze(matfile.two_samples(:,:,5));
        mytitle = sprintf('Two samples t-values uncorrected threshold');
        
        % 2D cluster and 1D correction for multiple testing
        % ------------------------------------------
    elseif MCC == 2
        if exist(MCC_data,'file')
            try
                H0_two_samples = load(MCC_data);
                H0_two_samples = H0_two_samples.H0_two_samples;
                bootT          = squeeze(H0_two_samples(:,:,1,:)); % get all T values under H0
                bootP          = squeeze(H0_two_samples(:,:,2,:)); % get all P values under H0
                if size(matfile.two_samples,1) == 1
                    tmp = NaN(1,size(matfile.two_samples,2),size(H0_two_samples,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
                clear H0_two_samples
                
                if size(M,1) == 1
                    [mask,M] = limo_clustering(M.^2,squeeze(matfile.two_samples(:,:,5)),bootT.^2,bootP,LIMO,3,p);
                else
                    [mask,M] = limo_clustering(M.^2,squeeze(matfile.two_samples(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
                end
                Nclust = unique(mask); Nclust = length(Nclust)-1; % mask = mask>0;
                if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
                mytitle = sprintf('Two Samples t-values cluster correction (%g %s)', Nclust, Mclust);
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
                H0_two_samples = load(MCC_data);
                H0_two_samples = H0_two_samples.H0_two_samples;
                bootT = squeeze(H0_two_samples(:,:,1,:)); % get all T values under H0
                if size(matfile.two_samples,1) == 1
                    tmp = NaN(1,size(matfile.two_samples,2),size(H0_two_samples,4));
                    tmp(1,:,:) = bootT; bootT = tmp; clear tmp
                end
                clear H0_two_samples
                [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
                mytitle = sprintf('Two Samples t-values correction by T max');
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
                H0_tfce_data = load(H0_tfce_data);
                [mask,M]     = limo_max_correction(tfce_data.tfce_two_samples, H0_tfce_data.tfce_H0_two_samples,p);
                mytitle      = sprintf('Two Samples t-values correction using TFCE');
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

% ---------------------
%% paired t-test
% --------------------

if strncmp(FileName,'paired_samples',14)
    
    %effect_nb = eval(FileName(32:end-4));
    if size(matfile.paired_samples,1)>1
        M = squeeze(matfile.paired_samples(:,:,4)); % T values
    else
        M = matfile.paired_samples(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = matfile.paired_samples(:,:,5) <= p;
        M = squeeze(matfile.paired_samples(:,:,5));
        mytitle = sprintf('Paired samples t-values uncorrected threshold');
        
        % 2D cluster and 1D correction for multiple testing
        % ------------------------------------------
    elseif MCC == 2
        if exist(MCC_data,'file')
            try
                H0_paired_samples = load(MCC_data);
                H0_paired_samples = H0_paired_samples.H0_paired_samples;
                bootT             = squeeze(H0_paired_samples(:,:,1,:)); % get all T values under H0
                bootP             = squeeze(H0_paired_samples(:,:,2,:)); % get all P values under H0
                if size(matfile.paired_samples,1) == 1
                    tmp = NaN(1,size(matfile.paired_samples,2),size(H0_paired_samples,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
                clear H0_paired_samples
                
                if size(M,1) == 1
                    [mask,M] = limo_clustering(M.^2,squeeze(matfile.paired_samples(:,:,5)),bootT.^2,bootP,LIMO,3,p);
                else
                    [mask,M] = limo_clustering(M.^2,squeeze(matfile.paired_samples(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
                end
                Nclust = unique(mask); Nclust = length(Nclust)-1; % mask = mask>0;
                if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
                mytitle = sprintf('Paired t-values cluster correction (%g %s)', Nclust, Mclust);
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
                H0_paired_samples = load(MCC_data);
                H0_paired_samples = H0_paired_samples.H0_paired_samples;
                bootT             = squeeze(H0_paired_samples(:,:,1,:)); % get all T values under H0
                if size(matfile.paired_samples,1) == 1
                    tmp = NaN(1,size(matfile.paired_samples,2),size(H0_paired_samples,4));
                    tmp(1,:,:) = bootT; bootT = tmp; clear tmp
                end
                [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
                mytitle = sprintf('Paired Samples t-values correction by T max');
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
                H0_tfce_data = load(H0_tfce_data);
                [mask,M]     = limo_max_correction(tfce_data.tfce_paired_samples, H0_tfce_data.tfce_H0_paired_samples,p);
                mytitle      = sprintf('Paired Samples t-values correction using TFCE');
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

% ------------------------
%% Repeated measure ANOVA
% ------------------------

if contains(FileName,'Rep_ANOVA')
    
    % all files have dim electrode x [freq/time] frames x F/p
    if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
        if contains(FileName,'Rep_ANOVA_Interaction')
            M    = matfile.Rep_ANOVA_Interaction_with_gp(:,:,:,1); % get F values
            PVAL = matfile.Rep_ANOVA_Interaction_with_gp(:,:,:,2); % get P values
        elseif contains(FileName,'Rep_ANOVA_Gp_effect')
            M    = matfile.Rep_ANOVA_Gp_effect(:,:,:,1); 
            PVAL = matfile.Rep_ANOVA_Gp_effect(:,:,:,2);
        elseif contains(FileName,'Rep_ANOVA_Main')
            M    = matfile.Rep_ANOVA(:,:,:,1); 
            PVAL = matfile.Rep_ANOVA(:,:,:,2);
        end
    else
        if contains(FileName,'Rep_ANOVA_Interaction')
            M    = matfile.Rep_ANOVA_Interaction_with_gp(:,:,1); 
            PVAL = matfile.Rep_ANOVA_Interaction_with_gp(:,:,2);
        elseif contains(FileName,'Rep_ANOVA_Gp_effect')
            M    = matfile.Rep_ANOVA_Gp_effect(:,:,1); 
            PVAL = matfile.Rep_ANOVA_Gp_effect(:,:,2);
        elseif contains(FileName,'Rep_ANOVA_Main')
            M    = matfile.Rep_ANOVA(:,:,1); 
            PVAL = matfile.Rep_ANOVA(:,:,2);
        end
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
                        tmp = NaN(1,size(M,2),size(bootT,4));
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
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    H0_Rep_ANOVA_Interaction_with_gp = load(MCC_data);
                    H0_Rep_ANOVA_Interaction_with_gp = H0_Rep_ANOVA_Interaction_with_gp.H0_Rep_ANOVA_Interaction_with_gp;
                    bootT                            = H0_Rep_ANOVA_Interaction_with_gp(:,:,1,:);
                    if size(matfile.Rep_ANOVA_Interaction_with_gp,1) == 1
                        tmp = NaN(1,size(matfile.Rep_ANOVA_Interaction_with_gp,2),size(H0_Rep_ANOVA_Interaction_with_gp,4));
                        tmp(1,:,:) = bootT; bootT = tmp;
                        clear tmp
                    end
                    clear H0_Rep_ANOVA_Interaction_with_gp
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    H0_Rep_ANOVA_Gp_effect = load(MCC_data);
                    H0_Rep_ANOVA_Gp_effect = H0_Rep_ANOVA_Gp_effect.H0_Rep_ANOVA_Gp_effect;
                    bootT = H0_Rep_ANOVA_Gp_effect(:,:,1,:);
                    if size(matfile.Rep_ANOVA_Gp_effect,1) == 1
                        tmp = NaN(1,size(matfile.Rep_ANOVA_Gp_effect,2),size(H0_Rep_ANOVA_Gp_effect,4));
                        tmp(1,:,:) = bootT; bootT = tmp;
                        clear tmp
                    end
                    clear H0_Rep_ANOVA_Gp_effect
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    H0_Rep_ANOVA = load(MCC_data);
                    H0_Rep_ANOVA = H0_Rep_ANOVA.H0_Rep_ANOVA;
                    bootT = H0_Rep_ANOVA(:,:,1,:); % get all F values under H0
                    if size(matfile.Rep_ANOVA,1) == 1
                        tmp = NaN(1,size(matfile.Rep_ANOVA,2),size(H0_Rep_ANOVA,4));
                        tmp(1,:,:) = bootT; bootT = tmp;
                        clear tmp
                    end
                    clear H0_Rep_ANOVA
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
                H0_tfce_data = load(H0_tfce_data);
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    [mask,M] = limo_max_correction(tfce_data.tfce_Rep_ANOVA_Interaction_with_gp, H0_tfce_data.tfce_H0_Rep_ANOVA_Interaction_with_gp,p);
                    mytitle = sprintf('Interaction correction using TFCE');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    [mask,M] = limo_max_correction(tfce_data.tfce_Rep_ANOVA_Gp_effect, H0_tfce_data.tfce_H0_Rep_ANOVA_Gp_effect,p);
                    mytitle = sprintf('Gp effect correction using TFCE');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    [mask,M] = limo_max_correction(tfce_data.tfce_Rep_ANOVA, H0_tfce_data.tfce_H0_Rep_ANOVA,p);
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

% -----------------------
%% Lateralization maps
% -----------------------
if strncmp(FileName,'LI_Map',6)
    
    M = squeeze(matfile.LI_Map(:,:,4)); % T values
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = matfile.LI_Map(:,:,5) <= p;
        mytitle = sprintf('LI Map T values');
        
        % 2D cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        MCC_data = sprintf('boot_%s',FileName);
        if exist(MCC_data,'file')
            try
                H0_LI_Map = load(MCC_data);
                H0_LI_Map = H0_LI_Map.H0_LI_Map;
                bootT     = squeeze(H0_LI_Map(:,:,2,:)); % get all T values under H0
                bootP     = squeeze(H0_LI_Map(:,:,3,:)); % get all P values under H0
                clear H0_LI_Map
                
                [mask,M]  = limo_clustering(M.^2,squeeze(matfile.LI(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
                Nclust = unique(mask); Nclust = length(Nclust)-1; % mask = mask>0;
                if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
                mytitle = sprintf('LI Map T-values cluster correction (%g %s)', Nclust, Mclust);
                
            catch ME
                errordlg(sprintf('error log: %s \n',ME.message),'cluster correction failure')
                return
            end
        else
            errordlg(['H0' filesep MCC_data ' not found'],'cluster correction failure')
            return
        end
        
        
        % T max correction for multiple testing
        % --------------------------------------
    elseif MCC == 4
        MCC_data = sprintf('boot_%s',FileName);
        if exist(MCC_data,'file')
            try
                H0_LI_Map = load(MCC_data);
                H0_LI_Map = H0_LI_Map.H0_LI_Map;
                bootT  = squeeze(H0_LI_Map(:,:,2,:)); % take all T values under H0
                clear H0_LI_Map
                [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
                mytitle = sprintf('LI Map T values correction by T max');
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
                H0_tfce_data = load(H0_tfce_data);
                [mask,M]     = limo_max_correction(tfce_data.tfce_LI, H0_tfce_data.tfce_H0_LI,p);
                mytitle      = sprintf('LI Map T-values correction using TFCE');
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

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
%  Copyright (C) LIMO Team 2019


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

% load data and set outputs to empty
% ----------------------------------
load(FileName);
M = [];
mask =[];
mytitle=[];
c = clock; disp(' ');

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

if strcmp(FileName,'R2.mat')
    M         = squeeze(R2(:,:,2)); % F values
    Pval      = squeeze(R2(:,:,3)); % P values
    MCC_data  = 'H0_R2.mat';
    titlename = 'R^2 Coef';
elseif strncmp(FileName,'Condition_effect',16)
    effect_nb = eval(FileName(18:end-4));
    M         = squeeze(Condition_effect(:,:,1));
    Pval      = squeeze(Condition_effect(:,:,2));
    MCC_data  = sprintf('H0_Condition_effect_%g.mat',effect_nb);
    titlename = sprintf('Condition effect %g  F values',effect_nb);
elseif strncmp(FileName,'Covariate_effect',16)
    effect_nb = eval( FileName(18:end-4));
    M         = squeeze(Covariate_effect(:,:,1));
    Pval      = squeeze(Covariate_effect(:,:,2));
    MCC_data  = sprintf('H0_Covariate_effect_%g.mat',effect_nb);
    titlename = sprintf('Covariate effect %g  F values',effect_nb);
elseif strncmp(FileName,'Interaction_effect',18)
    effect_nb = eval( FileName(20:end-4));
    M         = squeeze(Interaction_effect(:,:,1));
    Pval      = squeeze(Interaction_effect(:,:,2));
    MCC_data  = sprintf('H0_Interaction_effect_%g.mat',effect_nb);
    titlename = sprintf('Interaction effect %g  F values',effect_nb);
elseif strncmp(FileName,'semi_partial_coef',17)
    effect_nb = eval(FileName(19:end-4));
    M         = squeeze(semi_partial_coef(:,:,2));
    Pval      = squeeze(semi_partial_coef(:,:,3));
    MCC_data  = sprintf('H0_semi_partial_coef_%g.mat',effect_nb);
    titlename = sprintf('Semi Partial Coef %g',effect_nb);
elseif strncmp(FileName,'con_',4)
    effect_nb = eval(FileName(5:end-4));
    M         = squeeze(con(:,:,4));
    Pval      = squeeze(con(:,:,5));
    MCC_data  = sprintf('H0_con_%g.mat',effect_nb);
    titlename = sprintf('Contrast %g T values',effect_nb);
elseif strncmp(FileName,'ess_',4)
    effect_nb = eval(FileName(5:end-4));
    M         = squeeze(ess(:,:,end-1));
    Pval      = squeeze(ess(:,:,end));
    MCC_data  = sprintf('H0_ess_%g.mat',effect_nb);
    titlename = sprintf('Contrast %g F values',effect_nb);
end

% no correction for multiple testing
% -----------------------------------
if ~isempty(M) && MCC == 1
    
    if strcmp(LIMO.design.method,'WLS')
        if  all(isnan(Pval(:))) && exist(['H0' filesep MCC_data],'file')
            H0_data   = load(['H0' filesep MCC_data]);
            H0_data   = H0_data.(cell2mat(fieldnames(H0_data)));
            % get T/F values
            if strcmp(FileName,'R2.mat') 
                H0_values = squeeze(H0_data(:,:,2,:)); clear H0_data;
                [mask,M]  = limo_boot_threshold(M,H0_values,p);
                 R2(:,:,3) = M; save(FileName,'R2','-v7.3');
            else
                H0_values = squeeze(H0_data(:,:,end-1,:)); clear H0_data;
                [mask,M]  = limo_boot_threshold(M,H0_values,p);
                if strncmp(FileName,'Condition_effect',16)
                    Condition_effect(:,:,2) = M; save(FileName,'Condition_effect','-v7.3');
                elseif strncmp(FileName,'Covariate_effect',16)
                    Covariate_effect(:,:,2) = M; save(FileName,'Covariate_effect','-v7.3');
                elseif strncmp(FileName,'Interaction_effect',18)
                    Interaction_effect(:,:,2) = M; save(FileName,'Interaction_effect','-v7.3');
                elseif strncmp(FileName,'semi_partial_coef',17)
                    semi_partial_coef(:,:,2) = M; save(FileName,'semi_partial_coef','-v7.3');
                elseif strncmp(FileName,'con_',4)
                    con(:,:,5) = M; save(FileName,'con','-v7.3');
                elseif strncmp(FileName,'ess_',4)
                    ess(:,:,end) = M; save(FileName,'ess','-v7.3');
                end
            end
            mytitle   = sprintf('%s:\n uncorrected threshold using bootstrap',titlename);
        elseif all(isnan(Pval(:))) && ~exist(['H0' filesep MCC_data '.mat'],'file')
            M         = Pval;
            mask      = ones(size(M,1),size(M,2));
            mytitle   = sprintf('unthresholded %s:\n no p-values available without bootstrap',titlename);
        else % because it has been computed as above and stored on disk
            M         = Pval;
            mask      = Pval < p;
            mytitle   = sprintf('%s:\n uncorrected threshold using bootstrap',titlename);
        end
    else % OLS/IRLS
        M       = Pval;
        mask    = Pval < p;
        mytitle = sprintf('%s:\n uncorrected threshold',titlename);
    end
    
    % cluster correction for multiple testing
    % ---------------------------------------
elseif ~isempty(M) && MCC == 2
    if exist(['H0' filesep MCC_data],'file')
        try
            H0_data = load(['H0' filesep MCC_data]);
            H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
            if strcmp(FileName,'R2.mat') 
                bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
                bootP = squeeze(H0_data(:,:,3,:)); % get all P values under H0
            else
                bootM = squeeze(H0_data(:,:,1,:));
                bootP = squeeze(H0_data(:,:,2,:));
            end
            
            if all(isnan(bootP(:)))
                disp('WLS 1st pass - getting p values for null data')
                Pboot = cell(1,size(bootP,3));
                parfor boot = 1:size(bootP,3)
                    index = 1:size(bootP,3); index(boot) = [];
                    [~,Pboot{boot}]  = limo_boot_threshold(bootM(:,:,boot),bootM(:,:,index),p);
                end
                bootP = reshape(cell2mat(Pboot),size(bootM)); clear Pboot;
               
                % save it so we can reuse this anytime
                if strcmp(FileName,'R2.mat')
                    H0_R2 = H0_data; clear H0_data;
                    H0_R2(:,:,3,:) = bootP; 
                    save(['H0' filesep MCC_data],'H0_R2','-v7.3');
                elseif strncmp(FileName,'Condition_effect',16)
                    H0_Condition_effect = H0_data; clear H0_data;
                    H0_Condition_effect(:,:,2,:) = bootP; 
                    save(['H0' filesep MCC_data],'H0_Condition_effect','-v7.3');
                elseif strncmp(FileName,'Covariate_effect',16)
                    H0_Covariate_effect = H0_data; clear H0_data;
                    H0_Covariate_effect(:,:,2,:) = bootP; 
                    save(['H0' filesep MCC_data],'H0_Covariate_effect','-v7.3');
                elseif strncmp(FileName,'Interaction_effect',18)
                    H0_Interaction_effect = H0_data; clear H0_data;
                    H0_Interaction_effect(:,:,2,:) = bootP; 
                    save(['H0' filesep MCC_data],'H0_Interaction_effect','-v7.3');
                elseif strncmp(FileName,'semi_partial_coef',17)
                    H0_semi_partial_coef = H0_data; clear H0_data;
                    H0_semi_partial_coef(:,:,2,:) = bootP; 
                    save(['H0' filesep MCC_data],'H0_semi_partial_coef','-v7.3');
                elseif strncmp(FileName,'con_',4)
                    H0_con = H0_data; clear H0_data;
                    H0_con(:,:,2,:) = bootP; 
                    save(['H0' filesep MCC_data],'H0_con','-v7.3');
                elseif strncmp(FileName,'ess_',4)
                    H0_ess = H0_data; clear H0_data;
                    H0_ess(:,:,2,:) = bootP; 
                    save(['H0' filesep MCC_data],'H0_ess','-v7.3');
                end
                disp('null p-values saved')
            end
            
            % likely we don't have observed p-val
            if all(isnan(Pval(:)))
                [~,Pval] = limo_boot_threshold(M,bootM,p);
                if strcmp(FileName,'R2.mat')
                    R2(:,:,3) = Pval; save(FileName,'R2','-v7.3');
                elseif strncmp(FileName,'Condition_effect',16)
                    Condition_effect(:,:,2) = Pval; save(FileName,'Condition_effect','-v7.3');
                elseif strncmp(FileName,'Covariate_effect',16)
                    Covariate_effect(:,:,2) = Pval; save(FileName,'Covariate_effect','-v7.3');
                elseif strncmp(FileName,'Interaction_effect',18)
                    Interaction_effect(:,:,2) = Pval; save(FileName,'Interaction_effect','-v7.3');
                elseif strncmp(FileName,'semi_partial_coef',17)
                    semi_partial_coef(:,:,3) = Pval; save(FileName,'semi_partial_coef','-v7.3');
                elseif strncmp(FileName,'con_',4)
                    con(:,:,5) = M; save(FileName,'con','-v7.3');
                elseif strncmp(FileName,'ess_',4)
                    ess(:,:,end-1) = Pval; save(FileName,'ess','-v7.3');
                end
            end
            
            % finally get cluster mask and corrected p-values
            [mask,M] = limo_clustering(M,Pval,bootM,bootP,LIMO,MCC,p); % mask and cluster p values
            Nclust = unique(mask); Nclust = length(Nclust)-1; mask = mask>0;
            if Nclust <= 1; Mclust = 'cluster'; else ; Mclust = 'clusters'; end
            mytitle = sprintf('%s\n cluster correction (%g %s)', titlename, Nclust, Mclust);
        catch ME
            fprintf('someting went wrong loading H0 data %s \n',ME.message)
            return
        end
    else
        errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
        return
    end
    
    % correction using the max
    % --------------------------
elseif ~isempty(M) && MCC == 4 % Stat max
    if exist(['H0' filesep MCC_data],'file')
        try
            H0_data = load(['H0' filesep MCC_data]);
            H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
            if strcmp(FileName,'R2.mat') 
                bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
            else
                bootM = squeeze(H0_data(:,:,1,:));
            end
            clear H0_data;
            [mask,M] = limo_max_correction(abs(M),abs(bootM),p);
            mytitle  = sprintf('%s: \n correction by max',titlename);
        catch ME
            fprintf('someting went wrong loading H0 data %s \n',ME.message)
            return
        end
    else
        errordlg('no bootstrap file was found to compute the max distribution','missing data')
    end
    
    % correction using TFCE
    % --------------------------
elseif ~isempty(M) && MCC == 3 % Stat max
    if exist(['TFCE' filesep 'tfce_R2.mat'],'file')
        try
            score    = load(['TFCE' filesep 'tfce_' FileName]);
            H0_score = load(['H0' filesep 'tfce_' MCC_data]);
            [mask,M] = limo_max_correction(score.tfce_score,H0_score.tfce_H0_score,p);
            mytitle  = sprintf('%s:\n correction using TFCE',titlename);
        catch ME
            fprintf('someting went wrong loading H0 data %s \n',ME.message)
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
end

% -----------------------------------------
%% two samples t-test
% -----------------------------------------

if strncmp(FileName,'two_samples',11)
    
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
end

% ---------------------
%% paired t-test
% --------------------

if strncmp(FileName,'paired_samples',14)
    
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
end

% ------------------------
%% Repeated measure ANOVA
% ------------------------

if strncmp(FileName,'Rep_ANOVA',9)
    
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
end

% -----------------------
%% Lateralization maps
% -----------------------
if strncmp(FileName,'LI_Map',6)
    
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
end

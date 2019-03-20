function [M, mask, mytitle] = limo_stat_values(varargin)

% find corrected p values and mask from data under H0
%
% FORMAT [M, mask, mytitle] = limo_stat_values(varargin)
%
% INPUTS
%         Type = the type of plots (matters for title)
%         FileName = Name of the file selected
%         p = p value for thresholding
%         MCC = multiple comparisons option
%               1 none
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
% Cyril Pernet v1 25-05-2011
% Cyril Pernet v2 25-07-2012
% Marianne Latinus v3 2013 updated code for tfce
% Cyril Pernet v3 changed output + clean up
% Cyril Pernet v4 27-05-2015
%              removed depricated arguments, use 4 MCC, optimized H0
% ---------------------------------------------------------------------
% Copyright (C) LIMO Team 2015


Type      = varargin{1}; % type of plot
FileName  = varargin{2}; % Name of the file selected
p         = varargin{3}; % p value
MCC       = varargin{4}; % multiple comparison option
LIMO      = varargin{5}; % LIMO.mat

if isfield(LIMO,'Type')
    if strcmp(LIMO.Type,'Components') && MCC ~= 1, MCC = 4; end;
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
    
    M = squeeze(R2(:,:,2)); % F values
    MCC_data = 'H0_R2.mat';
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask = squeeze(R2(:,:,3)) < p;
        M = squeeze(R2(:,:,3)); % p values
        mytitle = sprintf('R^2 : uncorrected threshold');
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_R2(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_R2(:,:,2,:)); % get all P values under H0
            [mask,M] = limo_clustering(M,squeeze(R2(:,:,3)),bootM,bootP,LIMO,MCC,p); % mask and cluster p values
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('R^2: correction by \n spatial-temporal cluster');
            else
                mytitle = sprintf('R^2: correction by \n spatial-frequency cluster');
            end
        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
        
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_R2(:,:,2,:)); % get all F values under H0
            [mask,M] = limo_max_correction(M,bootM,p);
            mytitle = sprintf('R^2 : correction by F max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 3 % Stat max
        try cd('TFCE'); load('tfce_R2.mat'); cd ..
            cd('H0');load('tfce_H0_R2.mat'); cd ..
            [mask,M] = limo_max_correction(tfce_score,tfce_H0_score,p);
            mytitle = sprintf('R^2 : correction using TFCE');
        catch ME
            errordlg('no tfce bootstrap or tfce file was found to compute the tfce distribution','missing data')
            return
        end
    end
    
    
    
    % ---------------------------------
    %% Condition_effect (from 1st level)
    % ---------------------------------
    
elseif strncmp(FileName,'Condition_effect',16)
    
    effect_nb = eval(FileName(18:end-4));
    M = squeeze(Condition_effect(:,:,1)); % F values
    MCC_data = sprintf('H0_Condition_effect_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask = squeeze(Condition_effect(:,:,2)) < p;
        M = squeeze(Condition_effect(:,:,2));
        mytitle = sprintf('Condition %g: uncorrected threshold',effect_nb);
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Condition_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Condition_effect(:,:,2,:)); % get all P values under H0
            clear H0_Conditions;
            [mask,M] = limo_clustering(M,squeeze(Condition_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Condition %g: \n correction by spatial-temporal cluster',effect_nb);
            else
                mytitle = sprintf('Condition %g: \n correction by spatial-frequency cluster',effect_nb);
            end
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
            mytitle = sprintf('Condition %g: \n correction by F max',effect_nb);
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
            mytitle = sprintf('Condition %g: \n correction using TFCE',effect_nb);
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
    M = squeeze(Covariate_effect(:,:,1)); % F values
    MCC_data = sprintf('H0_Covariate_effect_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask = squeeze(Covariate_effect(:,:,2)) < p;
        M = squeeze(Covariate_effect(:,:,2)); % p values
        mytitle = sprintf('Covariate %g: uncorrected threshold ',effect_nb);
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Covariate_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Covariate_effect(:,:,2,:)); % get all P values under H0
            clear H0_Covariate_effect;
            [mask,M] = limo_clustering(M,squeeze(Covariate_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Covariate %g: \n correction by spatial-temporal cluster',effect_nb);
            else
                mytitle = sprintf('Covariate %g: \n correction by spatial-frequency cluster',effect_nb);
            end
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
            mytitle = sprintf('Covariate %g: \n correction by F max',effect_nb);
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
            mytitle = sprintf('Covariate %g: \n correction using TFCE',effect_nb);
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
    M = squeeze(Interaction_effect(:,:,1)); % F values
    MCC_data = sprintf('H0_Interaction_effect_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask = squeeze(Interaction_effect(:,:,2)) < p;
        M = squeeze(Interaction_effect(:,:,2));
        mytitle = sprintf('Interaction %g: uncorrected threshold',effect_nb);
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Interaction_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Interaction_effect(:,:,2,:)); % get all P values under H0
            clear H0_Interaction;
            [mask,M] = limo_clustering(M,squeeze(Interaction_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Interaction %g: \n correction by spatial-temporal cluster',effect_nb);
            else
                mytitle = sprintf('Interaction %g: \n correction by spatial-frequency cluster',effect_nb);
            end
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
            mytitle = sprintf('Interaction %g: \n correction by F max',effect_nb);
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
            mytitle = sprintf('Interaction %g: \n correction using TFCE',effect_nb);
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
    M = squeeze(semi_partial_coef(:,:,2)); % F values
    MCC_data = sprintf('H0_semi_partial_coef_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask = squeeze(semi_partial_coef(:,:,3)) < p;
        M = squeeze(semi_partial_coef(:,:,3));
        mytitle = sprintf('Semi partial coef %g: uncorrected threshold ',effect_nb);
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_semi_partial_coef(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_semi_partial_coef(:,:,2,:)); % get all P values under H0
            clear H0_semi_partial_coef;
            [mask,M] = limo_clustering(M,squeeze(semi_partial_coef(:,:,3)),bootM,bootP,LIMO,MCC,p);
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Semi partial coef %g: \n correction by spatial-temporal cluster',effect_nb);
            else
                mytitle = sprintf('Semi partial coef %g: \n correction by spatial-frequency cluster',effect_nb);
            end
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
            mytitle = sprintf('Semi partial coef %g: \n correction by F max',effect_nb);
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
            mytitle = sprintf('Semi partial coef %g: \n correction using TFCE',effect_nb);
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
    M = squeeze(con(:,:,4));
    MCC_data = sprintf('H0_con_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask = con(:,:,5) <= p;
        M = con(:,:,5);
        mytitle = sprintf('Contrast T %g: uncorrected threshold',effect_nb);
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        try cd H0; load(MCC_data); % dim electrode, frames, param/t/p, nboot
            bootT = H0_con(:,:,2,:); % get all t values under H0
            bootP = H0_con(:,:,3,:); % get all p values under H0
            clear H0_con
            [mask,M] = limo_clustering(M.^2,squeeze(con(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Contrast T %g: correction by spatial-temporal cluster', effect_nb);
            else
                mytitle = sprintf('Contrast T %g: correction by spatial-frequency cluster', effect_nb);
            end
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
            mytitle = sprintf('Contrast T %g: correction by T max', effect_nb);
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
            mytitle = sprintf('Contrast T %g: correction using TFCE', effect_nb);
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
    M = squeeze(ess(:,:,end-1));
    MCC_data = sprintf('H0_ess_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        mask = squeeze(ess(:,:,end)) < p;
        M = squeeze(ess(:,:,end));
        mytitle = sprintf('Contrast F %g: uncorrected threshold', effect_nb);
        
        
        % 1D 2D cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2
        
        try cd H0; load(MCC_data);
            bootF = squeeze(H0_ess(:,:,end-1,:));
            bootP = squeeze(H0_ess(:,:,end,:));
            clear H0_ess
            [mask,M] = limo_clustering(M,squeeze(ess(:,:,end)),bootF,bootP,LIMO,MCC,p);
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Contrast F %g: correction by spatial-temporal cluster', effect_nb);
            else
                mytitle = sprintf('Contrast F %g: correction by spatial-frequency cluster', effect_nb);
            end
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
            mytitle = sprintf('Contrast F %g: correction by F max', effect_nb);
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
            mytitle = sprintf('Contrast F %g: correction using TFCE', effect_nb);
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
        mask = one_sample(:,:,5) <= p;
        M = squeeze(one_sample(:,:,5));
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
            [mask,M] = limo_clustering(M.^2,squeeze(one_sample(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('One Sample t-test \n correction by spatial-temporal cluster');
            else
                mytitle = sprintf('One Sample t-test \n correction by spatial-frequency cluster');
            end
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
            mytitle = sprintf('One Sample t-test \n correction by T max');
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
            mytitle = sprintf('One Sample t-test \n correction using TFCE');
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
        
        mask = two_samples(:,:,5) <= p;
        M = squeeze(two_samples(:,:,5));
        mytitle = sprintf('Two samples t-test: uncorrected threshold');
        
        
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
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Two Samples t-test \n correction by spatial-temporal cluster');
            else
                mytitle = sprintf('Two Samples t-test \n correction by spatial-frequency cluster');
            end
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
            mytitle = sprintf('Two Samples t-test \n correction by T max');
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
            mytitle = sprintf('Two Samples t-test \n correction using TFCE');
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
        mytitle = sprintf('Paired samples t-test: \n uncorrected threshold');
        
        
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
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('Paired Samples t-test \n correction by spatial-temporal cluster');
            else
                mytitle = sprintf('Paired Samples t-test \n correction by spatial-frequency cluster');
            end
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
            mytitle = sprintf('Paired Samples t-test \n correction by T max');
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
            mytitle = sprintf('Paired Samples t-test \n correction using TFCE');
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
        M = Rep_ANOVA_Interaction_with_gp(:,:,1); % get the F values
        PVAL = Rep_ANOVA_Interaction_with_gp(:,:,2);
    elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
        M = Rep_ANOVA_Gp_effect(:,:,1); % get the F values
        PVAL = Rep_ANOVA_Gp_effect(:,:,2);
    elseif strncmp(FileName,'Rep_ANOVA',9)
        M = Rep_ANOVA(:,:,1); % get the F values
        PVAL = Rep_ANOVA(:,:,2);
    end
    
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = PVAL <= p;
        M = PVAL;
        if strncmp(FileName,'Rep_ANOVA_Interaction',21)
            mytitle = sprintf('Rep ANOVA Interaction: \n uncorrected threshold');
        elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
            mytitle = sprintf('Rep ANOVA Gp effect: \n uncorrected threshold');
        elseif strncmp(FileName,'Rep_ANOVA',9)
            mytitle = sprintf('Rep ANOVA: \n uncorrected threshold');
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
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                bootT = squeeze(H0_Rep_ANOVA_Gp_effect(:,:,1,:));
                bootP = squeeze(H0_Rep_ANOVA_Gp_effect(:,:,2,:));
                if size(Rep_ANOVA_Gp_effect,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA_Gp_effect,2),size(H0_Rep_ANOVA_Gp_effect,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
            elseif strncmp(FileName,'Rep_ANOVA',9)
                bootT = squeeze(H0_Rep_ANOVA(:,:,1,:)); % get all F values under H0
                bootP = squeeze(H0_Rep_ANOVA(:,:,2,:)); % get all P values under H0
                if size(Rep_ANOVA,1) == 1
                    tmp = NaN(1,size(Rep_ANOVA,2),size(H0_Rep_ANOVA,4));
                    tmp(1,:,:) = bootT; bootT = tmp;
                    tmp(1,:,:) = bootP; bootP = tmp;
                    clear tmp
                end
            end
            
            [mask,M] = limo_clustering(M,PVAL,bootT,bootP,LIMO,MCC,p);
            if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                if strcmp(LIMO.Analysis,'Time')
                    mytitle = sprintf('Rep ANOVA Interaction: \n correction by spatial-temporal cluster');
                else
                    mytitle = sprintf('Rep ANOVA Interaction: \n correction by spatial-frequency cluster');
                end
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                if strcmp(LIMO.Analysis,'Time')
                    mytitle = sprintf('Rep ANOVA Gp effect: \n correction by spatial-temporal cluster');
                else
                    mytitle = sprintf('Rep ANOVA Gp effect: \n correction by spatial-frequency cluster');
                end
            elseif strncmp(FileName,'Rep_ANOVA',9)
                if strcmp(LIMO.Analysis,'Time')
                    mytitle = sprintf('Rep ANOVA: \n correction by spatial-temporal cluster');
                else
                    mytitle = sprintf('Rep ANOVA: \n correction by spatial-frequency cluster');
                end
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
            if strcmp(LIMO.Analysis,'Time')
                mytitle = sprintf('LI Map: one sample T values \n correction by spatial-temporal cluster');
            else
                mytitle = sprintf('LI Map: one sample T values \n correction by spatial-frequency cluster');
            end
            
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
end

%% subfunctions which do the actual thresholding
% function [mask,cluster_p] = limo_clustering(M,P,bootM,bootP,LIMO,MCC,p)
% % call field trip functions to do 2D or 1D clustering
% %
% % M = 2D matrix of observed F values (note for a single electrode the format is 1*time frames*trials)
% % P = 2D matrix of observed p values (note for a single electrode the format is 1*time frames*trials)
% % bootM = 3D matrix of F values for data bootstrapped under H0
% % bootP = 3D matrix of F values for data bootstrapped under H0
% % LIMO = LIMO structure - information requested is LIMO.data.chanlocs and LIMO.data.neighbouring_matrix
% % MCC = 2 (spatial-temporal clustering) or 3 (temporal clustering)
% % p = threshold to apply (note this applied to create clusters and to
% % threshold the cluster map)
% 
% if size(M,1) == 1
%     MCC = 21;
% end
% cluster_p = [];
% mask = [];
% 
% if MCC == 2
%     nboot = size(bootM,3);
%     U = round((1-p)*nboot); % bootstrap threshold
%     if size(bootM,1)>1 % many electrodes
%         minnbchan = 2;
%         expected_chanlocs = LIMO.data.chanlocs;
%         channeighbstructmat = LIMO.data.neighbouring_matrix;
%         boot_maxclustersum=zeros(nboot,1); % compute bootstrap clusters
%         disp('getting clusters under H0 boot ...');
%         parfor boot = 1:nboot
%             % boot_maxclustersum(boot) = limo_getclustersum(bootM(:,:,boot),bootP(:,:,boot),channeighbstructmat,minnbchan,p);
%             [posclusterslabelmat,nposclusters] = limo_findcluster((bootP(:,:,boot) <= p),channeighbstructmat,minnbchan);
%             
%             bootM_b = bootM(:,:,boot);
%             if nposclusters~=0
%                 tmp=zeros(1,nposclusters);
%                 for C = 1:nposclusters % compute sum for each cluster
%                     tmp(C) = sum( bootM_b(posclusterslabelmat==C) );
%                 end
%                 boot_maxclustersum(boot) = max(tmp(:)); % save max across clusters
%             else
%                 boot_maxclustersum(boot) = 0;
%             end
%         end
%         [mask, cluster_p] = limo_cluster_test(M,P,boot_maxclustersum,channeighbstructmat,minnbchan,p);
%     
%     elseif size(bootM,1)==1 % one electrode
%         th = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);
%         sigcluster = limo_ecluster_test(squeeze(M),squeeze(P),th,p);
%         mask = sigcluster.elec_mask; cluster_p = [];
%     end
%     
% elseif MCC == 21
%     nboot = size(bootM,3);
%     U = round((1-p)*nboot); % bootstrap threshold
%     th = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);
%     sigcluster = limo_ecluster_test(squeeze(M),squeeze(P),th,p);
%     mask = sigcluster.elec_mask; cluster_p = [];
%     
% end
% end

% function [mask,p_val] = limo_max_correction(M,bootM,p)
% correction for multiple testing using the max stat value
% note this works for bootstrapped data under H0 and for TFCE
%
% M = 2D matrix of observed values (note for a single electrode the format is 1*time frames*trials)
% bootM = 3D matrix of F values for data bootstrapped under H0
% p = threshold to apply

% nboot = size(bootM,3);
% for boot=1:nboot
%     data = squeeze(bootM(:,:,boot));
%     maxM(boot) = max(data(:)); % collect highest absolute value in space and time for each boot
% end
% 
% U=round((1-p).*nboot);
% sortmaxM = sort(maxM); % sort bootstraps
% maxF_th = sortmaxM(U); % get threshold for each parameter
% mask = squeeze(M) >= maxF_th;
% % figure; imagesc(mask)
% for row =1:size(M,1)
%     for column=1:size(M,2)
%         p_val(row,column) = 1-(sum(squeeze(M(row,column)) >=sortmaxM) / nboot);
%     end
% end
% 
% end



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
%               2 2D clustering
%               3 1D clusteriung
%               4 uses max values
%               5 uses TFCE
%         LIMO = LIMO.mat
%         choice = if MCC = 1 but bootstrap was used then we have the choice between empirical or theoretical p values
%         choice = 'use empirical p values' or 'use theoretical p values'
%         effect_nb = indicates wich effect or regressor to plot
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
% --------------------------------------------------
%  Copyright (C) LIMO Team 2010


Type      = varargin{1}; % type of plot
FileName  = varargin{2}; % Name of the file selected
p         = varargin{3}; % p value
MCC       = varargin{4}; % multiple comparison option
LIMO      = varargin{5}; % LIMO.mat
choice = [];
effect_nb = [];
if nargin >5
    choice = varargin{6}; % if MCC = 1 use empirical or theoretical p values
    if nargin == 7
        effect_nb = varargin{7}; % indicates wich effect nb to plot
    end
end

load (FileName);
M = []; mask =[]; mytitle=[];
c = clock; disp(' ');

if MCC ~= 1
    fprintf('limo_display_results %gh %gmin %gsec: computing statistical correction...\n',c(4),c(5),c(6));
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
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
        
        if strcmp(choice,'use empirical p values')
            try cd('H0');load(MCC_data); cd ..
                H0_F_values = squeeze(H0_R2(:,:,2,:)); clear H0_R2;
                sorted_values = sort(H0_F_values,3); clear H0_F_values
                U = round((1-p)*size(sorted_values,3));
                mask = (M >= sorted_values(:,:,U));
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
                    end
                end
                M = 1- (tmp ./ size(sorted_values,3)) ; % p values
                mytitle = sprintf('R^2 : uncorrected threshold \n using bootstraped F values');
            catch ME
                mask = squeeze(R2(:,:,3)) < p;
                M = squeeze(R2(:,:,3)); % p values
                mytitle = sprintf('R^2 : uncorrected threshold');
            end
        else
            mask = squeeze(R2(:,:,3)) < p;
            M = squeeze(R2(:,:,3)); % p values
            mytitle = sprintf('R^2 : uncorrected threshold');
        end
        
        
        % cluster correction for multiple testing
        % ---------------------------------------  
    elseif MCC == 2 || MCC == 3
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_R2(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_R2(:,:,2,:)); % get all P values under H0
            [mask,M] = local_clustering(M,squeeze(R2(:,:,3)),bootM,bootP,LIMO,MCC,p); % mask and cluster p values
            if MCC == 2
                mytitle = sprintf('R^2: correction by \spatial-temporal cluster');
            elseif MCC == 3
                mytitle = sprintf('R^2: correction by temporal cluster');
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
            [mask,M] = max_correction(M,bootM,p);
            mytitle = sprintf('R^2 : correction by F max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
         
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max
        try cd('TFCE'); load('tfce_R2.mat'); cd ..
            cd('H0');load('tfce_H0_R2.mat'); cd ..
            [mask,M] = max_correction(tfce_score,tfce_H0_score,p);
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
    
    effect_nb = eval( FileName(18:end-4));
    M = squeeze(Condition_effect(:,:,1)); % F values
    MCC_data = sprintf('H0_Condition_effect_%g',effect_nb);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
        
        if strcmp(choice,'use empirical p values')
            try cd('H0');load(MCC_data); cd ..
                H0_F_values = squeeze(H0_Condition_effect(:,:,1,:)); clear H0_Condition;
                sorted_values = sort(H0_F_values,3); clear H0_F_values
                U = round((1-p)*size(sorted_values,3));
                mask = M >= sorted_values(:,:,U);
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
                    end
                end
                M = 1- (tmp ./ size(sorted_values,3)) ; % p values
                mytitle = sprintf('Condition %g: uncorrected threshold \n using bootstraped F values',effect_nb);
            catch ME
                mask = squeeze(Condition_effect(:,:,2)) < p;
                M = squeeze(Condition_effect(:,:,2)); % p values
                mytitle = sprintf('Condition %g: uncorrected threshold',effect_nb);
            end
        else
            mask = squeeze(Condition_effect(:,:,2)) < p;
            M = squeeze(Condition_effect(:,:,2));
            mytitle = sprintf('Condition %g: uncorrected threshold',effect_nb);
        end
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2 || MCC == 3

        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Condition_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Condition_effect(:,:,2,:)); % get all P values under H0
            clear H0_Conditions;
            [mask,M] = local_clustering(M,squeeze(Condition_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
            if MCC == 2
                mytitle = sprintf('Condition %g: \n correction by spatial-temporal cluster',effect_nb);
            elseif MCC == 3
                mytitle = sprintf('Condition %g: \n correction by temporal cluster',effect_nb);
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
            clear H0_Condition_effect; [mask,M] = max_correction(M,bootM,p);
            mytitle = sprintf('Condition %g: \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            clear H0_Condition_effect; [mask,M] = max_correction(tfce_score,tfce_H0_score,p);
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
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
        
        if strcmp(choice,'use empirical p values')
            try cd('H0');load(MCC_data); cd ..
                H0_F_values = squeeze(H0_Covariate_effect(:,:,1,:)); clear H0_Covariates;
                sorted_values = sort(H0_F_values,3); clear H0_F_values
                U = round((1-p)*size(sorted_values,3));
                mask = M >= sorted_values(:,:,U);
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
                    end
                end
                M = 1- (tmp ./ size(sorted_values,3)) ; % p values
                mytitle = sprintf('Covariate %g: uncorrected threshold \n using bootstraped F values',effect_nb);
            catch ME
                mask = squeeze(Covariate_effect(:,:,2)) < p;
                M = squeeze(Covariate_effect(:,:,2)); % p values
                mytitle = sprintf('Covariate %g: uncorrected threshold ',effect_nb);
            end
        else
            mask = squeeze(Covariate_effect(:,:,2)) < p;
            M = squeeze(Covariate_effect(:,:,2)); % p values
            mytitle = sprintf('Covariate %g: uncorrected threshold ',effect_nb);
        end
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2 || MCC == 3
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Covariate_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Covariate_effect(:,:,2,:)); % get all P values under H0
            clear H0_Covariate_effect;
            [mask,M] = local_clustering(M,squeeze(Covariate_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
            if MCC == 2
                mytitle = sprintf('Covariate %g: \n correction by spatial-temporal cluster',effect_nb);
            elseif MCC == 3
                mytitle = sprintf('Covariate %g: \n correction by temporal cluster',effect_nb);
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
            [mask,M] = max_correction(M,bootM,p);
            mytitle = sprintf('Covariate %g: \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max
        
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            [mask,M] = max_correction(tfce_score,tfce_H0_score,p);
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
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
        
        if strcmp(choice,'use empirical p values')
            try cd('H0');load(MCC_data); cd ..
                H0_F_values = squeeze(H0_Interaction(:,:,1,:)); clear H0_Interaction;
                sorted_values = sort(H0_F_values,3); clear H0_F_values
                U = round((1-p)*size(sorted_values,3));
                mask = M >= sorted_values(:,:,U);
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
                    end
                end
                M = 1- (tmp ./ size(sorted_values,3)) ; % p values
                mytitle = sprintf('Interaction %g: uncorrected threshold \n using bootstraped F values',effect_nb);
            catch ME
                mask = squeeze(Interaction_effect(:,:,2)) < p;
                M = squeeze(Interaction_effect(:,:,2)); % p values
                mytitle = sprintf('Interaction %g: uncorrected threshold',effect_nb);
            end
        else
            mask = squeeze(Interaction_effect(:,:,2)) < p;
            M = squeeze(Interaction_effect(:,:,2));
            mytitle = sprintf('Interaction %g: uncorrected threshold',effect_nb);
        end
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2 || MCC == 3
        
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_Interaction_effect(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_Interaction_effect(:,:,2,:)); % get all P values under H0
            clear H0_Interaction;
            [mask,M] = local_clustering(M,squeeze(Interaction_effect(:,:,2)),bootM,bootP,LIMO,MCC,p);
            if MCC == 2
                mytitle = sprintf('Interaction %g: \n correction by spatial-temporal cluster',effect_nb);
            elseif MCC == 3
                mytitle = sprintf('Interaction %g: \n correction by temporal cluster',effect_nb);
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
            clear H0_Interaction_effect; [mask,M] = max_correction(M,bootM,p);
            mytitle = sprintf('Interaction %g: \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            clear H0_Interaction_effect; [mask,M] = max_correction(tfce_score,tfce_H0_score,p);
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
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
        
        if strcmp(choice,'use empirical p values')
            try cd('H0'); load(MCC_data); cd ..
                H0_F_values = squeeze(H0_semi_partial_coef(:,:,2,:)); % get F values
                clear H0_semi_partial_coef;
                sorted_values = sort(H0_F_values,3);
                clear H0_F_values
                U = round((1-p)*size(sorted_values,3));
                mask = M >= sorted_values(:,:,U);
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
                    end
                end
                M = 1- (tmp ./ size(sorted_values,3)) ; % p values
                mytitle = sprintf('Semi partial coef %g: uncorrected threshold \n using on bootstraped F values',effect_nb);
            catch ME
                mask = squeeze(semi_partial_coef(:,:,3)) < p; % simply threshold p values
                M = squeeze(semi_partial_coef(:,:,3));
                mytitle = sprintf('Semi partial coef %g: uncorrected threshold ',effect_nb);
            end
        else
            mask = squeeze(semi_partial_coef(:,:,3)) < p;
            M = squeeze(semi_partial_coef(:,:,3));
            mytitle = sprintf('Semi partial coef %g: uncorrected threshold ',effect_nb);
        end
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2 || MCC == 3
        try cd('H0');load(MCC_data); cd ..
            bootM = squeeze(H0_semi_partial_coef(:,:,1,:)); % get all F values under H0
            bootP = squeeze(H0_semi_partial_coef(:,:,2,:)); % get all P values under H0
            clear H0_semi_partial_coef;
            [mask,M] = local_clustering(M,squeeze(semi_partial_coef(:,:,3)),bootM,bootP,LIMO,MCC,p);
            if MCC == 2
                mytitle = sprintf('Semi partial coef %g: \n correction by spatial-temporal cluster',effect_nb);
            elseif MCC == 3
                mytitle = sprintf('Semi partial coef %g: \n correction by temporal cluster',effect_nb);
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
            [mask,M] = max_correction(M,bootM,p);
            mytitle = sprintf('Semi partial coef %g: \n correction by F max',effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            clear H0_semi partial coef; [mask,M] = max_correction(tfce_score,tfce_H0_score,p);
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
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
        
        if strcmp(choice,'use empirical p values')
            try cd H0; load(MCC_data);
                H0_T_values  = squeeze(boot_H0_con(:,:,2,:)); % T values under H0
                sorted_values = sort(H0_T_values,3);
                clear boot_H0_con H0_F_values
                U = round((1-p)*size(sorted_values,3));
                mask = M >= sorted_values(:,:,U);
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
                    end
                end
                M = min((tmp ./ size(sorted_values,3)), 1- (tmp ./ size(sorted_values,3))) ; % p values
            catch ME
                mask = con(:,:,5) <= p;
                M = con(:,:,5);
                mytitle = sprintf('Contrast T %g: uncorrected threshold',effect_nb);
            end
        else
            mask = con(:,:,5) <= p;
            M = con(:,:,5);
            mytitle = sprintf('Contrast T %g: uncorrected threshold',effect_nb);
        end

        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2 || MCC == 3
        try cd H0; load(MCC_data); % dim electrode, frames, param/t/p, nboot
            bootT = H0_con(:,:,2,:); % get all t values under H0
            bootP = H0_con(:,:,3,:); % get all p values under H0
            clear H0_con
            [mask,M] = local_clustering(M.^2,squeeze(con(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
            if MCC == 2
                mytitle = sprintf('Contrast T %g: correction by spatial-temporal cluster', effect_nb);
            elseif MCC == 3
                mytitle = sprintf('Contrast T %g: correction by temporal cluster', effect_nb);
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
            [mask,M] = max_correction(abs(M),abs(bootT),p); % absolute T values
            mytitle = sprintf('Contrast T %g: correction by T max', effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
        
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            [mask,M] = max_correction(tfce_score,tfce_H0_score,p);
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
    M = squeeze(ess(:,:,end-1)); 
    MCC_data = sprintf('H0_ess_%g',effect_nb);
        
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
        
        if strcmp(choice,'use empirical p values')
            try
                cd H0; load(MCC_data);
                % sort all F values
                sorted_values = sort(squeeze(boot_H0_ess(:,:,end-1,:)),3);
                U = round((1-p)*LIMO.design.nboot);
                mask = ess(:,:,end-1) >= sorted_values(:,:,U);
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(sorted_values(row,column,:)));
                    end
                end
                M = 1- (tmp ./ size(sorted_values,3)) ; % p values
                mytitle = sprintf('Contrast F %g: threshold using bootstrapped F values',effect_nb);
            catch ME
                mask = squeeze(ess(:,:,end)) < p;
                M = squeeze(ess(:,:,end));
                mytitle = sprintf('Contrast F %g: uncorrected threshold', effect_nb);
            end
        else
            mask = squeeze(ess(:,:,end)) < p;
            M = squeeze(ess(:,:,end));
            mytitle = sprintf('Contrast F %g: uncorrected threshold', effect_nb);
        end
        
       
        % 1D 2D cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2 || MCC == 3
        
        try cd H0; load(MCC_data);
                bootF = squeeze(H0_ess(:,:,end-1,:)); 
                bootP = squeeze(H0_ess(:,:,end,:)); 
            clear H0_ess
            [mask,M] = local_clustering(M,squeeze(ess(:,:,end)),bootF,bootP,LIMO,MCC,p); 
            if MCC == 2
                mytitle = sprintf('Contrast F %g: correction by spatial-temporal cluster', effect_nb);
            elseif MCC == 3
                mytitle = sprintf('Contrast F %g: correction by temporal cluster', effect_nb);
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
            [mask,M] = max_correction(M,bootF,p); % absolute T values
            mytitle = sprintf('Contrast F %g: correction by F max', effect_nb);
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end

        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max
        try cd('TFCE'); tfceName = ['tfce_' FileName]; load(tfceName); cd ..
            cd('H0'); tfceName = ['tfce_H0_' FileName]; load(tfceName); cd ..
            [mask,M] = max_correction(tfce_score,tfce_H0_score,p);
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
    
    effect_nb = eval(FileName(28:end-4));
    if size(one_sample,1)>1
        M = squeeze(one_sample(:,:,4)); % T values
    else
        M = one_sample(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName); 
     
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        try
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
            
            if strcmp(choice,'use empirical p values')
                lo = round((LIMO.design.bootstrap.*p)/2);
                hi = LIMO.design.bootstrap - lo;
                load(MCC_data)
                Tsorted  = sort(squeeze(H0_one_sample(:,:,1,:)),3); % sort T values under H0 along bootstrap dimension
                if size(one_sample,1) == 1
                    tmp = NaN(1,size(one_sample,2),size(H0_one_sample,4));
                    tmp(1,:,:) = Tsorted; Tsorted = tmp;
                end
                TCI(:,:,1) = Tsorted(:,:,lo+1);
                TCI(:,:,2) = Tsorted(:,:,hi);
                mask = M >= TCI(:,:,2) | M <= TCI(:,:,1);
                
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(Tsorted(row,column,:)));
                    end
                end
                M = min(tmp./size(Tsorted,3), 1-(tmp./size(Tsorted,3))) ; % p values
                mytitle = sprintf('One sample t-test \n thresholding using bootstrapped T values');
            else
                mask = one_sample(:,:,5) <= p;
                M = squeeze(one_sample(:,:,5));
                mytitle = sprintf('One sample t-test: uncorrected threshold');
            end
        catch ME
            mask = one_sample(:,:,5) <= p;
            M = squeeze(one_sample(:,:,5));
            mytitle = sprintf('One sample t-test: uncorrected threshold');
        end
        
        
   % 2D cluster and 1D correction for multiple testing
   % ------------------------------------------
   elseif MCC == 2 || MCC == 3
       if size(M,1) == 1; MCC =3; end
        try load(MCC_data);
            bootT = squeeze(H0_one_sample(:,:,1,:)); % get all T values under H0
            bootP = squeeze(H0_one_sample(:,:,2,:)); % get all P values under H0
            if size(one_sample,1) == 1
                tmp = NaN(1,size(one_sample,2),size(H0_one_sample,4));
                tmp(1,:,:) = bootT; bootT = tmp;
                tmp(1,:,:) = bootP; bootP = tmp; 
                clear tmp
            end
            [mask,M] = local_clustering(M.^2,squeeze(one_sample(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
            if MCC == 2
                mytitle = sprintf('One Sample t-test \n correction by spatial-temporal cluster');
            elseif MCC == 3
                mytitle = sprintf('One Sample t-test \n correction by temporal cluster');
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
            [mask,M] = max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('One Sample t-test \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
            
    % Correction using TFCE
    % -------------------------------------    
    elseif MCC == 5 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
        	load(MCC_tfce_data)
            [mask,M] = max_correction(tfce_one_sample, tfce_H0_one_sample,p);
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
    
    effect_nb = eval(FileName(29:end-4));
    if size(two_samples,1)>1
        M = squeeze(two_samples(:,:,4)); % T values
    else
        M = two_samples(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        try
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
            
            if strcmp(choice,'use empirical p values')
                load(MCC_data)
                lo = round((LIMO.design.bootstrap.*p)/2);
                hi = LIMO.design.bootstrap - lo;
                Tsorted  = sort(squeeze(H0_two_samples(:,:,1,:)),3); % sort T values under H0 along bootstrap dimension
                if size(two_samples,1) == 1
                    tmp = NaN(1,size(two_samples,2),size(H0_two_samples,4));
                    tmp(1,:,:) = Tsorted; Tsorted = tmp;
                end
                TCI(:,:,1) = Tsorted(:,:,lo+1);
                TCI(:,:,2) = Tsorted(:,:,hi);
                mask = M >= TCI(:,:,2) | M <= TCI(:,:,1);
                
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(Tsorted(row,column,:)));
                    end
                end
                M = min(tmp./size(Tsorted,3), 1-(tmp./size(Tsorted,3))) ; % p values
                mytitle = sprintf('Two samples t-test \n thresholding using bootstrapped T values');
            else
                mask = two_samples(:,:,5) <= p;
                M = squeeze(two_samples(:,:,5));
                mytitle = sprintf('Two samples t-test: uncorrected threshold');
            end
        catch ME
            mask = two_samples(:,:,5) <= p;
            M = squeeze(two_samples(:,:,5));
            mytitle = sprintf('Two samples t-test: uncorrected threshold');
        end
        
        
   % 2D cluster and 1D correction for multiple testing
   % ------------------------------------------
   elseif MCC == 2 || MCC == 3
       if size(M,1) == 1; MCC =3; end
        try load(MCC_data);
            bootT = squeeze(H0_two_samples(:,:,1,:)); % get all T values under H0
            bootP = squeeze(H0_two_samples(:,:,2,:)); % get all P values under H0
            if size(two_samples,1) == 1
                tmp = NaN(1,size(two_samples,2),size(H0_two_samples,4));
                tmp(1,:,:) = bootT; bootT = tmp;
                tmp(1,:,:) = bootP; bootP = tmp; 
                clear tmp
            end
            [mask,M] = local_clustering(M.^2,squeeze(two_samples(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
            if MCC == 2
                mytitle = sprintf('Two Samples t-test \n correction by spatial-temporal cluster');
            elseif MCC == 3
                mytitle = sprintf('Two Samples t-test \n correction by temporal cluster');
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
            [mask,M] = max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('Two Samples t-test \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
            
    % Correction using TFCE
    % -------------------------------------    
    elseif MCC == 5 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
        	load(MCC_tfce_data)
            [mask,M] = max_correction(tfce_two_samples, tfce_H0_two_samples,p);
            mytitle = sprintf('Two Samples t-test \n correction using TFCE');
        catch ME
            errordlg('no tfce bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end    
    
    % paired t-test
    % --------------------
    
elseif strncmp(FileName,'paired_samples',14)
    
    effect_nb = eval(FileName(32:end-4));
    if size(paired_samples,1)>1
        M = squeeze(paired_samples(:,:,4)); % T values
    else
        M = paired_samples(1,:,4);
    end
    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        try
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
            
            if strcmp(choice,'use empirical p values')
                load(MCC_data)
                lo = round((LIMO.design.bootstrap.*p)/2);
                hi = LIMO.design.bootstrap - lo;
                Tsorted  = sort(squeeze(H0_paired_samples(:,:,1,:)),3); % sort T values under H0 along bootstrap dimension
                if size(paired_samples,1) == 1
                    tmp = NaN(1,size(paired_samples,2),size(H0_paired_samples,4));
                    tmp(1,:,:) = Tsorted; Tsorted = tmp;
                end
                TCI(:,:,1) = Tsorted(:,:,lo+1);
                TCI(:,:,2) = Tsorted(:,:,hi);
                mask = M >= TCI(:,:,2) | M <= TCI(:,:,1);
                
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(Tsorted(row,column,:)));
                    end
                end
                M = min(tmp./size(Tsorted,3), 1-(tmp./size(Tsorted,3))) ; % p values
                mytitle = sprintf('Paired samples t-test \n thresholding using bootstrapped T values');
            else
                mask = paired_samples(:,:,5) <= p;
                M = squeeze(paired_samples(:,:,5));
                mytitle = sprintf('Paired samples t-test: \n uncorrected threshold');
            end
        catch ME
            mask = paired_samples(:,:,5) <= p;
            M = squeeze(paired_samples(:,:,5));
            mytitle = sprintf('Paired samples t-test: \n uncorrected threshold');
        end
        
        
   % 2D cluster and 1D correction for multiple testing
   % ------------------------------------------
   elseif MCC == 2 || MCC == 3
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
            [mask,M] = local_clustering(M.^2,squeeze(paired_samples(:,:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
            if MCC == 2
                mytitle = sprintf('Paired Samples t-test \n correction by spatial-temporal cluster');
            elseif MCC == 3
                mytitle = sprintf('Paired Samples t-test \n correction by temporal cluster');
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
            [mask,M] = max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('Paired Samples t-test \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
            
    % Correction using TFCE
    % -------------------------------------    
    elseif MCC == 5 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
        	load(MCC_tfce_data)
            [mask,M] = max_correction(tfce_paired_samples, tfce_H0_paired_samples,p);
            mytitle = sprintf('Paired Samples t-test \n correction using TFCE');
        catch ME
            errordlg('no tfce bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end    
        
    
    % Repeated measure ANOVA
    % --------------------
    
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
        try
        if LIMO.design.bootstrap == 1 && isempty(choice)
            choice = questdlg('an associated bootstrap should be there','p value choice','use theoretical p values','use empirical p values','use empirical p values');
        end
            
            if strcmp(choice,'use empirical p values')
                load(MCC_data)
                lo = round((LIMO.design.bootstrap.*p)/2);
                hi = LIMO.design.bootstrap - lo;
                
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    Tsorted  = sort(squeeze(Rep_ANOVA_Interaction_with_gp(:,:,1,:)),3); % sort T values under H0 along bootstrap dimension
                    if size(Rep_ANOVA_Interaction_with_gp,1) == 1
                        tmp = NaN(1,size(Rep_ANOVA_Interaction_with_gp,2),size(H0_Rep_ANOVA_Interaction_with_gp,4));
                        tmp(1,:,:) = Tsorted; Tsorted = tmp;
                    end
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    Tsorted  = sort(squeeze(Rep_ANOVA_Gp_effect(:,:,1,:)),3); % sort T values under H0 along bootstrap dimension
                    if size(Rep_ANOVA_Gp_effect,1) == 1
                        tmp = NaN(1,size(Rep_ANOVA_Gp_effect,2),size(H0_Rep_ANOVA_Gp_effect,4));
                        tmp(1,:,:) = Tsorted; Tsorted = tmp;
                    end
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    Tsorted  = sort(squeeze(Rep_ANOVA(:,:,1,:)),3); % sort T values under H0 along bootstrap dimension
                    if size(Rep_ANOVA,1) == 1
                        tmp = NaN(1,size(Rep_ANOVA,2),size(H0_Rep_ANOVA,4));
                        tmp(1,:,:) = Tsorted; Tsorted = tmp;
                    end
                end
                TCI(:,:,1) = Tsorted(:,:,lo+1);
                TCI(:,:,2) = Tsorted(:,:,hi);
                mask = M >= TCI(:,:,2) | M <= TCI(:,:,1);
                
                for row = 1:size(M,1)
                    for column = 1:size(M,2)
                        tmp(row,column) = sum(M(row,column)>squeeze(Tsorted(row,column,:)));
                    end
                end
                M = 2*min(tmp./size(Tsorted,3), 1-(tmp./size(Tsorted,3))) ; % p values
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    mytitle = sprintf('Rep ANOVA Interaction: \n thresholding using bootstrapped T values');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    mytitle = sprintf('Rep ANOVA Gp effect: \n thresholding using bootstrapped T values');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    mytitle = sprintf('Rep ANOVA: \n thresholding using bootstrapped T values');
                end
            else
                mask = PVAL <= p;
                M = PVAL;
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    mytitle = sprintf('Rep ANOVA Interaction: \n uncorrected threshold');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    mytitle = sprintf('Rep ANOVA Gp effect: \n uncorrected threshold');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    mytitle = sprintf('Rep ANOVA: \n uncorrected threshold');
                end
            end
        catch ME
            mask = PVAL <= p;
            M = PVAL;
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    mytitle = sprintf('Rep ANOVA Interaction: \n uncoorected threshold');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    mytitle = sprintf('Rep ANOVA Gp effect: \n uncoorected threshold');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    mytitle = sprintf('Rep ANOVA: \n uncoorected threshold');
                end
        end
        
        
        % cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2 || MCC == 3
       if size(M,1) == 1; MCC =3; end
       
       try load(MCC_data);
           if strncmp(FileName,'Rep_ANOVA_Interaction',21)
               bootT = H0_Rep_ANOVA_Interaction_with_gp(:,:,1,:);
               bootP = H0_Rep_ANOVA_Interaction_with_gp(:,:,2,:);
               if size(Rep_ANOVA_Interaction_with_gp,1) == 1
                   tmp = NaN(1,size(Rep_ANOVA_Interaction_with_gp,2),size(H0_Rep_ANOVA_Interaction_with_gp,4));
                   tmp(1,:,:) = bootT; bootT = tmp;
                   tmp(1,:,:) = bootP; bootP = tmp;
                   clear tmp
               end
           elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
               bootT = H0_Rep_ANOVA_Gp_effect(:,:,1,:);
               bootP = H0_Rep_ANOVA_Gp_effect(:,:,2,:);
               if size(Rep_ANOVA_Gp_effect,1) == 1
                   tmp = NaN(1,size(Rep_ANOVA_Gp_effect,2),size(H0_Rep_ANOVA_Gp_effect,4));
                   tmp(1,:,:) = bootT; bootT = tmp;
                   tmp(1,:,:) = bootP; bootP = tmp;
                   clear tmp
               end
           elseif strncmp(FileName,'Rep_ANOVA',9)
               bootT = H0_Rep_ANOVA(:,:,1,:); % get all F values under H0
               bootP = H0_Rep_ANOVA(:,:,2,:); % get all P values under H0
               if size(Rep_ANOVA,1) == 1
                   tmp = NaN(1,size(Rep_ANOVA,2),size(H0_Rep_ANOVA,4));
                   tmp(1,:,:) = bootT; bootT = tmp;
                   tmp(1,:,:) = bootP; bootP = tmp;
                   clear tmp
               end
           end
            
            [mask,M] = local_clustering(M,PVAL,bootT,bootP,LIMO,MCC,p); 
            if MCC == 2
                if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    mytitle = sprintf('Rep ANOVA Interaction: \n correction by spatial-temporal cluster');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    mytitle = sprintf('Rep ANOVA Gp effect: \n correction by spatial-temporal cluster');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    mytitle = sprintf('Rep ANOVA: \n correction by spatial-temporal cluster');
                end
            elseif MCC == 3
                 if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                    mytitle = sprintf('Rep ANOVA Interaction: \n correction by temporal cluster');
                elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                    mytitle = sprintf('Rep ANOVA Gp effect: \n correction by temporal cluster');
                elseif strncmp(FileName,'Rep_ANOVA',9)
                    mytitle = sprintf('Rep ANOVA: \n correction by temporal cluster');
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
           
           [mask,M] = max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
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
    elseif MCC == 5 % Stat tfce
        MCC_tfce_data = sprintf('H0%stfce_H0_%s', filesep, FileName);
        tfce_data = sprintf('tfce%stfce_%s',filesep, FileName);
        try load(tfce_data);
        	load(MCC_tfce_data)

            if strncmp(FileName,'Rep_ANOVA_Interaction',21)
                [mask,M] = max_correction(tfce_Rep_ANOVA_Interaction_with_gp, tfce_H0_Rep_ANOVA_Interaction_with_gp,p);
            elseif strncmp(FileName,'Rep_ANOVA_Gp_effect',19)
                [mask,M] = max_correction(tfce_Rep_ANOVA_Gp_effect, tfce_H0_Rep_ANOVA_Gp_effect,p);
            elseif strncmp(FileName,'Rep_ANOVA',9)
                [mask,M] = max_correction(tfce_Rep_ANOVA, tfce_H0_Rep_ANOVA,p);
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
        try
            MCC_data = sprintf('boot_%s',FileName); load(MCC_data);
            if isempty(choice)
                choice = questdlg('an associated bootstrap file has been found','p value choice','use theoretical p values','use empirical p values','use empirical p values');
            end
            
            if strcmp(choice,'use empirical p values')
                
                lo = round((LIMO.design.nboot.*p)/2);
                hi = LIMO.design.nboot - lo;
                Tsorted  = sort(squeeze(boot_LI_Map(:,:,1,:)),3); % sort T values under H1 along bootstrap dimension
                if size(M,1) > 1
                    TCI(:,:,1) = Tsorted(:,:,lo+1);
                    TCI(:,:,2) = Tsorted(:,:,hi);
                    mask = (TCI(:,:,1) > 0) + (TCI(:,:,2) < 0); % sig if does not include 0
                else
                    TCI(:,1) = Tsorted(:,lo+1);
                    TCI(:,2) = Tsorted(:,hi);
                    mask = (TCI(:,1)' >0) + (TCI(:,2)' < 0);
                end
                mytitle = sprintf('LI Map: one sample T values \n threshold based on bootstrap T values');
            else
                mask = LI_Map(:,:,5) <= p;
                mytitle = sprintf('LI Map: one sample T values \n threshold based on theoretical p values');
                
            end
        catch ME
            mask = LI_Map(:,:,5) <= p;
            mytitle = sprintf('LI Map: one sample T values \n threshold based on theoretical p values');
        end
        
        
        % 2D cluster correction for multiple testing
        % ---------------------------------------
    elseif MCC == 2  
        
        MCC_data = sprintf('boot_%s',FileName);
        try load(MCC_data);
            bootT = squeeze(boot_LI_Map(:,:,2,:)); % get all T values under H0
            bootP = squeeze(boot_LI_Map(:,:,3,:)); % get all P values under H0
            U = round((1-p)*LIMO.design.nboot); % bootstrap threshold
            
            if size(bootT,1)>1 % many electrodes
                
                minnbchan = 2; % we take at leat two channels to make a cluster
                expected_chanlocs = LIMO.data.chanlocs;
                channeighbstructmat = LIMO.data.neighbouring_matrix;
                boot_maxclustersum=zeros(LIMO.design.nboot,1); % compute bootstrap clusters
                for s=1:LIMO.design.nboot
                    boot_maxclustersum(s) = limo_getclustersum(bootT(:,:,s).^2,bootP(:,:,s),channeighbstructmat,minnbchan,p);
                end
                sort_boot_maxclustersum = sort(boot_maxclustersum,1);
                mask = limo_cluster_test(LI_Map(:,:,4).^2,LI_Map(:,:,5),sort_boot_maxclustersum(U),channeighbstructmat,minnbchan,p);
                mytitle = sprintf('LI Map: one sample T values \n correction by spatial-temporal cluster');
                
            elseif size(booT,1)==1 % one electrode
                th = limo_ecluster_make( squeeze(bootT).^2,squeeze(bootP),p );
                sigcluster = limo_ecluster_test( squeeze(one_sample(:,:,4)).^2,squeeze(one_sample(:,:,5)),th,p );
                mask = sigcluster.elec;
                mytitle = sprintf('LI Map: one sample T values \n correction by temporal cluster');
            end
        catch ME
            errordlg('no bootstrap file was found to compute the cluster correction','missing data')
            return
        end
        
        
    elseif MCC == 3
        
        MCC_data = sprintf('boot_%s',FileName);
        try load(MCC_data);
            bootT = squeeze(boot_LI_Map(:,:,2,:)); % get all T values under H0
            bootP = squeeze(boot_LI_Map(:,:,3,:)); % get all P values under H0
            U = round((1-p)*LIMO.design.nboot); % bootstrap threshold
            th = limo_ecluster_make( squeeze(bootT).^2,squeeze(bootP),p );
            sigcluster = limo_ecluster_test( squeeze(one_sample(:,:,4)).^2,squeeze(one_sample(:,:,5)),th,p );
            mask = sigcluster.elec;
            mytitle = sprintf('LI Map: one sample T values \n correction by temporal cluster');
        catch ME
            errordlg('no bootstrap file was found to compute the cluster correction','missing data')
            return
        end
        
        
        % T max correction for multiple testing
        % --------------------------------------
    elseif MCC == 4
        MCC_data = sprintf('boot_%s',FileName);
        try load(MCC_data);
            T  = squeeze(boot_LI_Map(:,:,2,:)); % take all T values under H0
            for s=1:LIMO.design.nboot
                if length(size(T)) == 2
                    maxT(s) = max(abs(squeeze(T(:,s)))); % take max across frames for each bootstrap
                else
                    maxT(s) = max(max(abs(squeeze(T(:,:,s))))); % take max across electrodes and frames for each bootstrap
                end
            end
            
            U=round((1-p).*LIMO.design.nboot);
            sortmaxT = sort(maxT);
            maxT_th = sortmaxT(U);
            mask = squeeze(abs(LI_Map(:,:,4))) >= maxT_th;
            mytitle = sprintf('One sample T values \n correction by T max');
        catch ME
            errordlg('no bootstrap file was found to compute the max distribution','missing data')
            return
        end
    end
    
else
    errordlg2('unidentified FileName - no thresholding done');
end
end

%% subfunctions which do the actual thresholding
function [mask,cluster_p] = local_clustering(M,P,bootM,bootP,LIMO,MCC,p)
% call field trip functions to do 2D or 1D clustering
%
% M = 3D matrix of observed F values (note for a single electrode the format is 1*time frames*trials)
% P = 3D matrix of observed p values (note for a single electrode the format is 1*time frames*trials)
% bootM = 3D matrix of F values for data bootstrapped under H0
% bootP = 3D matrix of F values for data bootstrapped under H0
% LIMO = LIMO structure - information requested is LIMO.data.chanlocs and LIMO.data.neighbouring_matrix
% MCC = 2 (spatial-temporal clustering) or 3 (temporal clustering)
% p = threshold to apply (note this applied to create clusters and to
% threshold the cluster map)

if size(M,1) == 1
    MCC = 3;
end
cluster_p = [];
mask = [];

if MCC == 2 
    nboot = size(bootM,3);
    U = round((1-p)*nboot); % bootstrap threshold
    if size(bootM,1)>1 % many electrodes
        minnbchan = 2;
        expected_chanlocs = LIMO.data.chanlocs;
        channeighbstructmat = LIMO.data.neighbouring_matrix;
        boot_maxclustersum=zeros(nboot,1); % compute bootstrap clusters
        for boot=1:nboot
            boot_maxclustersum(boot) = limo_getclustersum(bootM(:,:,boot),bootP(:,:,boot),channeighbstructmat,minnbchan,p);
        end
        [mask, cluster_p] = limo_cluster_test(M,P,boot_maxclustersum,channeighbstructmat,minnbchan,p);

    elseif size(bootM,1)==1 % one electrode
        th = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);
        sigcluster = limo_ecluster_test(squeeze(M),squeeze(P),th,p);
        mask = sigcluster.elec; cluster_p = [];
    end
    
elseif MCC == 3
    nboot = size(bootM,3);
    U = round((1-p)*nboot); % bootstrap threshold
    th = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);
    sigcluster = limo_ecluster_test(squeeze(M),squeeze(P),th,p);
    mask = sigcluster.elec;
    
end
end

function [mask,p_val] = max_correction(M,bootM,p)
% correction for multiple testing using the max stat value
% note this works for bootstrapped data under H0 and for TFCE
%
% M = 3D matrix of observed values (note for a single electrode the format is 1*time frames*trials)
% bootM = 3D matrix of F values for data bootstrapped under H0
% p = threshold to apply 

nboot = size(bootM,3);
for boot=1:nboot
    data = squeeze(bootM(:,:,boot));
    maxM(boot) = max(data(:)); % collect highest absolute value in space and time for each boot
end

U=round((1-p).*nboot);
sortmaxM = sort(maxM); % sort bootstraps 
maxF_th = sortmaxM(U); % get threshold for each parameter
mask = squeeze(M) >= maxF_th; 
% figure; imagesc(mask)
for row =1:size(M,1)
    for column=1:size(M,2)
        p_val(row,column) = 1-(sum(squeeze(M(row,column)) >=sortmaxM) / nboot);
    end
end 
 
end


    
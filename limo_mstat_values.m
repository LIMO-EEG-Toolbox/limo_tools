function [M, mask, mytitle] = limo_mstat_values(varargin)

% find corrected p values and mask from data under H0
% this file is more multivariate statistics only
%
% FORMAT [M, mask, mytitle] = limo_mstat_values(varargin)
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
%         M the stat values
%         mask the significant data
%         mytitle is the right title to use to make the figure for type 1 and 2
%
% see limo_display_results limo_stat_values
%
% Cyril Pernet v1 13-06-2013
% ------------------------------
%  Copyright (C) LIMO Team 2019


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
        effect_nb = varargin{7}; % for regressions only, indicates wich regressor to plot
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

    % ---------
    %% R2.mat 
    % ---------
    
if strcmp(FileName,'R2.mat')
    
    if strcmp(choice,'Roy')
        M = squeeze(R2(:,2)); % F values
    else
        M = squeeze(R2(:,4));
    end
    MCC_data = 'H0_R2.mat';
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if LIMO.design.bootstrap == 1
            try cd('H0');load(MCC_data); cd ..
                if strcmp(choice,'Roy')
                    H0_F_values = squeeze(H0_R2(:,2,:)); clear H0_R2;
                else
                    H0_F_values = squeeze(H0_R2(:,4,:)); clear H0_R2;
                end
                sorted_values = sort(H0_F_values,2); clear H0_F_values
                U = round((1-p)*size(sorted_values,2));
                mask = (M >= sorted_values(:,U));
                for column = 1:size(M,2)
                    tmp(column) = sum(M(column)>squeeze(sorted_values(column,:)));
                end
                M = 1- (tmp ./ size(sorted_values,2)) ; % p values
                mytitle = sprintf('R^2 : uncorrected threshold \n based on bootstrapped F values');
            catch ME
                if strcmp(choice,'Roy')
                    mask = R2(:,3) < p;
                    M = squeeze(R2(:,3)); % p values
                    mytitle = sprintf('R^2 : uncorrected threshold \n based on Roy test');
                else
                    mask = R2(:,5) < p;
                    M = squeeze(R2(:,5)); % p values
                    mytitle = sprintf('R^2 : uncorrected threshold \n based on Pillai test');
                end
            end
        else
            if strcmp(choice,'Roy')
                mask = R2(:,3) < p;
                M = squeeze(R2(:,3)); % p values
                mytitle = sprintf('R^2 : uncorrected threshold \n based on Roy test');
            else
                mask = R2(:,5) < p;
                M = squeeze(R2(:,5)); % p values
                mytitle = sprintf('R^2 : uncorrected threshold \n based on Pillai test');
            end
        end
        mask = mask' ; M = M'; 
        
        % cluster correction for multiple testing
        % ---------------------------------------  
    elseif MCC == 2 || MCC == 3

        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        
         
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max

    end 

    
    
    % ---------------------------------
    %% Condition_effect (from 1st level) 
    % ---------------------------------
    
elseif strncmp(FileName,'Condition_effect',16)
    
    effect_nb = eval(FileName(18:end-4));
    MCC_data = sprintf('H0_Condition_effect_%g', effect_nb); 
    
    if strcmp(choice,'Roy')
        M = squeeze(Condition_effect(:,1)); % F values Roy
        P = squeeze(Condition_effect(:,2)); % p values Roy
    else
        M = squeeze(Condition_effect(:,3)); % F values Pillai
        P = squeeze(Condition_effect(:,4)); % p values Pillai
    end    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1  
        if LIMO.design.bootstrap >= 1
            try cd('H0');load(MCC_data); cd ..
                if strcmp(choice,'Roy')
                    H0_F_values = squeeze(H0_Condition_effect(:,1,:)); clear H0_Condition_effect;
                else
                    H0_F_values = squeeze(H0_Condition_effect(:,3,:)); clear H0_Condition_effect;
                end
                sorted_values = sort(H0_F_values,2); clear H0_F_values
                U = round((1-p)*size(sorted_values,2));
                mask = (M >= sorted_values(:,U));
                for column = 1:size(M,2)
                    tmp(column) = sum(M(column)>squeeze(sorted_values(column,:)));
                end
                M = 1- (tmp ./ size(sorted_values,2)) ; % p values
                mytitle = sprintf('Condition effect : uncorrected threshold \n based on bootstrapped F values');
            catch ME
                if strcmp(choice,'Roy')
                    mask = Condition_effect(:,2) < p;
                    M = squeeze(Condition_effect(:,2)); % p values
                    mytitle = sprintf('Condition effect : uncorrected threshold \n based on Roy test');
                else
                    mask = Condition_effect(:,4) < p;
                    M = squeeze(Condition_effect(:,4)); % p values
                    mytitle = sprintf('Condition effect : uncorrected threshold \n based on Pillai test');
                end
            end
        else
            if strcmp(choice,'Roy')
                mask = Condition_effect(:,2) < p;
                M = squeeze(Condition_effect(:,2)); % p values
                mytitle = sprintf('Condition effect : uncorrected threshold \n based on Roy test');
            else
                mask = Condition_effect(:,4) < p;
                M = squeeze(Condition_effect(:,4)); % p values
                mytitle = sprintf('Condition effect : uncorrected threshold \n based on Pillai test');
            end
        end
        mask = mask' ; M = M'; 
        
        % cluster correction for multiple testing
        % ---------------------------------------  
    elseif MCC == 2 || MCC == 3
       try cd('H0'); load(MCC_data); cd ..
            if strcmp(choice,'Roy')
                bootM = squeeze(H0_Condition_effect(:,1,:)); % get all F values 
                bootP = squeeze(H0_Condition_effect(:,2,:)); % get all p-values     
            elseif strcmp(choice,'Pillai')
                bootM = squeeze(H0_Condition_effect(:,3,:)); % get all F values pillai
                bootP = squeeze(H0_Condition_effect(:,4,:)); % get all p-values pillai
            end
            % do the cluster correction:
            [mask, M] = limo_clustering(M', P', bootM, bootP, LIMO, MCC, p, 1); 
            mytitle = sprintf('Condition effect : correction by temporal clustering \n ');

        catch ME
            errordlg('no bootstrap file was found')
            return
        end
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        try cd('H0'); load(MCC_data); cd ..
            if strcmp(choice,'Roy')
                bootM = squeeze(H0_Condition_effect(:,1,:)); % get all F values 
                bootP = squeeze(H0_Condition_effect(:,2,:)); % get all p-values     
            elseif strcmp(choice,'Pillai')
                bootM = squeeze(H0_Condition_effect(:,3,:)); % get all F values pillai
                bootP = squeeze(H0_Condition_effect(:,4,:)); % get all p-values pillai
            end
            % correction:
            [mask, M] = limo_max_correction(M, bootM, p,1); 
            mask = mask';
            mytitle = sprintf('Condition effect : correction by max F \n ');

        catch ME
            errordlg('no bootstrap file was found')
            return
        end
         
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max

    end 
    % ---------------------------------
    %% Linear_Classification (from 1st level) 
    % ---------------------------------
    
elseif strncmp(FileName,'Linear_Classification',16)
    
    MCC_data = sprintf('H0_Linear_Classification'); 
    
    M = Linear_Classification(:,2); % observed CV classification accuracies
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1  
        if LIMO.design.bootstrap >= 1
            try cd('H0');load(MCC_data); cd ..
                H0_values = squeeze(H0_Linear_Classification);
                sorted_values = sort(H0_values,2); 
                U = round((1-p)*size(sorted_values,2));
                mask = (M >= sorted_values(:,U));
                for column = 1:size(M,1)
                    tmp(column) = sum(M(column)>squeeze(sorted_values(column,:))); % how many times H0 bootstraps higher than observed value
                end
                M = 1- (tmp ./ size(sorted_values,2)) ; % p values
                mytitle = sprintf('CV linear decoding accuracies ± 2SD: uncorrected threshold \n based on H0 bootstrapped classification');
            catch ME
                mask = Linear_Classification(:,2) >  1/LIMO.design.nb_conditions;
                M = NaN(size(Linear_Classification(:,2))); % p values
                mytitle = sprintf('CV linear decoding accuracies ± 2SD: no threshold, classification > chance');
            end
        else
                mask = Linear_Classification(:,2) >  1/LIMO.design.nb_conditions;
                M = NaN(size(Linear_Classification(:,2))); % p values
                mytitle = sprintf('CV linear decoding accuracies ± 2SD: no threshold, classification > chance');
        end
        mask = mask' ; M = M'; 
        
        % cluster correction for multiple testing
        % ---------------------------------------  
    elseif MCC == 2 || MCC == 3
       try cd('H0'); load(MCC_data); cd ..
            bootM = H0_Linear_Classification; 
            bootP = NaN(size(H0_Linear_Classification));            
            % do the cluster correction:
            % 1.make the clusters under H0
            b = size(bootM,2);
            U = round((1-p)*b);
            boot_values = zeros(b,1);
            for kk=1:b % bootstrap samples
                [L,NUM] = bwlabeln(squeeze(bootM(:,kk))>=1/LIMO.design.nb_conditions); % find clusters
                if NUM~=0
                    tmp=zeros(1,NUM);
                    for C = 1:NUM % compute sum for each cluster
                        tmp(C) = sum(squeeze(bootM(L==C,kk)) );
                    end
                    boot_values(kk) = max(tmp(:)); % save max across clusters
                else
                    boot_values(kk) = 0;
                end

            end % bootstrap loop
            sortSC = sort(boot_values);
            threshold = sortSC(U); % threshold 
            % 2.make clusters under H1 and threshold
            sigcluster = zeros(size(M,1),1);
            pval = NaN(size(M,1),1);
            [L,NUM] = bwlabeln(Linear_Classification(:,2) >= 1/LIMO.design.nb_conditions); % find clusters
            maxval = zeros(1,NUM);
            for C = 1:NUM % compute cluster sums & compare to bootstrap threshold
                maxval(C) = sum(M(L==C)); % sum of classification accuracies
                if maxval(C) >= threshold;
                    sigcluster(L==C)=1; % flag clusters above threshold
                    p = sum(boot_values >= sum(M(L==C))) / b; % corrrected p value
                    if p ==0
                        p = 1/b;
                    end
                    pval(L==C) = p;
                end
            end
            maxval = max(maxval);
            % return output
            mask = sigcluster';
            M = pval;
            %[mask, M] = limo_clustering(M', M', bootM, bootM, LIMO, MCC, 1/LIMO.design.nb_conditions, 1); 
            mytitle = sprintf('CV linear decoding accuracies ± 2SD: correction by temporal clustering \n ');

        catch ME
            errordlg('no bootstrap file was found')
            return
        end
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        try cd('H0'); load(MCC_data); cd ..
            bootM = H0_Linear_Classification; 
            % correction:
            [mask, M] = limo_max_correction(M, bootM, p,1); 
            mask = mask';
            mytitle = sprintf('CV linear decoding accuracies ± 2SD: correction by max statistic\n ');

        catch ME
            errordlg('no bootstrap file was found')
            return
        end
         
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max

    end    
    % ---------------------------------
    %% Covariate_effect (from 1st level) 
    % ---------------------------------
    
elseif strncmp(FileName,'Covariate_effect',16)
    
    if strcmp(choice,'Roy')
        M = squeeze(R2(:,2)); % F values
    else
        M = squeeze(R2(:,4));
    end
    MCC_data = 'H0_R2.mat';
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        
        if LIMO.design.bootstrap == 1
            try cd('H0');load(MCC_data); cd ..
                if strcmp(choice,'Roy')
                    H0_F_values = squeeze(H0_R2(:,2,:)); clear H0_R2;
                else
                    H0_F_values = squeeze(H0_R2(:,4,:)); clear H0_R2;
                end
                sorted_values = sort(H0_F_values,2); clear H0_F_values
                U = round((1-p)*size(sorted_values,2));
                mask = (M >= sorted_values(:,U));
                for column = 1:size(M,2)
                    tmp(column) = sum(M(column)>squeeze(sorted_values(column,:)));
                end
                M = 1- (tmp ./ size(sorted_values,2)) ; % p values
                mytitle = sprintf('R^2 : uncorrected threshold \n based on bootstrapped F values');
            catch ME
                if strcmp(choice,'Roy')
                    mask = R2(:,3) < p;
                    M = squeeze(R2(:,3)); % p values
                    mytitle = sprintf('R^2 : uncorrected threshold \n based on Roy test');
                else
                    mask = R2(:,5) < p;
                    M = squeeze(R2(:,5)); % p values
                    mytitle = sprintf('R^2 : uncorrected threshold \n based on Pillai test');
                end
            end
        else
            if strcmp(choice,'Roy')
                mask = R2(:,3) < p;
                M = squeeze(R2(:,3)); % p values
                mytitle = sprintf('R^2 : uncorrected threshold \n based on Roy test');
            else
                mask = R2(:,5) < p;
                M = squeeze(R2(:,5)); % p values
                mytitle = sprintf('R^2 : uncorrected threshold \n based on Pillai test');
            end
        end
        mask = mask' ; M = M'; 
        
        % cluster correction for multiple testing
        % ---------------------------------------  
    elseif MCC == 2 || MCC == 3
                
        % correction using the max
        % --------------------------
    elseif MCC == 4 % Stat max
        
         
        % correction using TFCE
        % --------------------------
    elseif MCC == 5 % Stat max

    end 
    % ------------------------------------------
    %% One sample t-test
    % ------------------------------------------
    
elseif strncmp(FileName,'one_sample',10)
    
    try effect_nb = eval(FileName(28:end-4)); end

    M = one_sample(:,4); % T values

    MCC_data = sprintf('H0%sH0_%s', filesep, FileName);
    
    % no correction for multiple testing
    % -----------------------------------
    if MCC == 1
        mask = one_sample(:,5) <= p;
        %M = squeeze(one_sample(:,5));
        mytitle = sprintf('One sample t-test on classification accuracies \n uncorrected threshold');
        
        
        % 2D cluster and 1D correction for multiple testing
        % ------------------------------------------
    elseif MCC == 2
        try load(MCC_data);
            bootT = squeeze(H0_one_sample(:,1,:)); % get all T values under H0
            bootP = squeeze(H0_one_sample(:,2,:)); % get all P values under H0
            [mask,M] = limo_clustering(M.^2,squeeze(one_sample(:,5)),bootT.^2,bootP,LIMO,MCC,p); % square T values
            mytitle = sprintf('One sample t-test results on classification accuracies \n correction by temporal cluster');

        catch ME
            errordlg('no bootstrap file was found to compute the cluster distribution','missing data')
            return
        end
                
        % T max correction for multiple testing
        % -------------------------------------
    elseif MCC == 4 % Stat max
        try load(MCC_data)
            bootT = squeeze(H0_one_sample(:,1,:)); % get all T values under H0
            [mask,M] = limo_max_correction(abs(M),abs(bootT),p); % threshold max absolute T values
            mytitle = sprintf('One sample t-test results on classification accuracies \n correction by T max');
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
        
else
    errordlg2('unidentified FileName - no thresholding done');
end
end





    
function limo_tfce_handling(varargin)

% routine for create tfce files commensurate to boostrapped files
%
% FORMAT  limo_tfce_handling(filename,'checkfile','yes')
%
% INPUTS filename is the stat file that need to be tfce (if H0 exist it is done too)
%        'checkfiles' is 'yes' by default - if 'no' it overwrites without asking
%
% OUTPUTS tfce_* files are created on the drive
% -----------------------------------------------
%  Copyright (C) LIMO Team 2020

% inputs
file = varargin{1};
[filepath,filename,ext] = fileparts(file);
if isempty(filepath); filepath = pwd; end
if isempty(ext)
    file = dir(fullfile(filepath,[filename '*']));
    filename = file.name;
end

if exist(fullfile(filepath,[filename ext]),'file')
    if ~exist(fullfile(filepath,'LIMO.mat'),'file')
        error('no LIMO.mat found next to %s',file)
    else
        filename = [filename ext];
        tmp  = load(fullfile(filepath,'LIMO.mat'));
        LIMO = tmp.LIMO;
        clear tmp
    end
else
    error('can''t find %s',varargin{1})
end

checkfile = 'yes';
for i=2:2:nargin
    if strcmpi(varargin{i},'checkfile')
        checkfile = varargin{i+1};
    end
end

tfce_file    = fullfile(LIMO.dir,['tfce' filesep 'tfce_' filename]);
H0_tfce_file = fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' filename]);

% quick user check
% ----------------
if strcmpi(checkfile,'yes')
    if exist(tfce_file,'file')
        answer = questdlg('tfce file already exist - overwrite?','data check','Yes','No','Yes');
        if strcmp(answer,'Yes')
            LIMO.design.tfce = 1;
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        else
            LIMO.design.tfce = 0;
            return
        end
    else
        LIMO.design.tfce = 1;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        if ~exist(fullfile(LIMO.dir,'tfce'),'dir')
            mkdir(fullfile(LIMO.dir,'tfce'));
        end
    end
end


% simply tfce the file
nboot = LIMO.design.bootstrap;
mkdir(fullfile(LIMO.dir,'tfce'))
fprintf('Thresholding %s using TFCE \n',filename);

% -------------------------------------------------------------------------
if contains(filename,'Covariate_effect_') % continuous regressors, any level
    if LIMO.Level ==1 && ~isfield(LIMO.data,'neighbouring_matrix')
        warning('no neighbouring matrix found, this is required for TFCE')
        [~, LIMO.data.neighbouring_matrix] = limo_expected_chanlocs;
        if isempty(LIMO.data.neighbouring_matrix)
            return
        else
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
        end
    end
    
    Covariate_effect = load(filename);
    Covariate_effect = Covariate_effect.(cell2mat(fieldnames(Covariate_effect)));
    if size(Covariate_effect,1) == 1
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score(1,:,:) = limo_tfce(2,squeeze(Covariate_effect(:,:,:,1)),[]);
        else
            tfce_score(1,:) = limo_tfce(1,squeeze(Covariate_effect(:,:,1)),LIMO.data.neighbouring_matrix);
        end
    else
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score = limo_tfce(3,squeeze(Covariate_effect(:,:,:,1)),LIMO.data.neighbouring_matrix);
        else
            tfce_score = limo_tfce(2,squeeze(Covariate_effect(:,:,1)),LIMO.data.neighbouring_matrix);
        end
    end
    save(tfce_file,'tfce_score','-v7.3');
    clear Covariate_effect tfce_score;
    
    H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
    if exist(H0filename,'file')
        fprintf('Applying TFCE to null data ... \n')
        H0_Covariate_effect = load(H0filename);
        H0_Covariate_effect = H0_Covariate_effect.(cell2mat(fieldnames(H0_Covariate_effect)));
        if size(H0_Covariate_effect,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),size(H0_Covariate_effect,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,:,b) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            else
                tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),size(H0_Covariate_effect,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            else
                tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            end
        end
    end
    save(H0_tfce_file,'tfce_H0_score','-v7.3');
    
    % -------------------------------------------------------------------------
elseif ~contains(filename,'Covariate_effect_') && LIMO.Level == 1
    % -------------------------------------------------------------------------

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
    
    % R2
    if contains(filename,'R2')
        R2 = load(fullfile(LIMO.dir,'R2.mat'));
        R2 = R2.(cell2mat(fieldnames(R2)));
        if size(R2,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_score(1,:,:) = limo_tfce(2, squeeze(R2(:,:,:,2)),[]); % no neighbouring, time-freq cluster
            else
                tfce_score(1,:)   = limo_tfce(1, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_score        = limo_tfce(3, squeeze(R2(:,:,:,2)),LIMO.data.neighbouring_matrix);
            else
                tfce_score        = limo_tfce(2, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
            end
        end
        save(tfce_file,'tfce_score','-v7.3');
        clear R2;
        
        H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
        if exist(H0filename,'file')
            fprintf('Applying TFCE to null data ... \n')
            H0_R2 = load(H0filename);
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
                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_H0_score = NaN(size(H0_R2,1),size(H0_R2,2),size(H0_R2,3),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_R2(:,:,:,2,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score = NaN(size(H0_R2,1),size(H0_R2,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            end
            save(H0_tfce_file,'tfce_H0_score','-v7.3');
            clear H0_R2 tfce_H0_score;
        end
    end
    
    % conditions
    if contains(filename,'Condition_effect_')
        Condition_effect = load(filename);
        Condition_effect = Condition_effect.(cell2mat(fieldnames(Condition_effect)));
        if size(Condition_effect,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_score(1,:,:) = limo_tfce(2, squeeze(Condition_effect(:,:,:,1)),[]);
            else
                tfce_score(1,:) = limo_tfce(1, squeeze(Condition_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_score      = limo_tfce(3, squeeze(Condition_effect(:,:,:,1)),LIMO.data.neighbouring_matrix);
            else
                tfce_score      = limo_tfce(2, squeeze(Condition_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            end
        end
        save(tfce_file,'tfce_score','-v7.3');
        clear Condition_effect tfce_score;
        
        H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
        if exist(H0filename,'file')
            fprintf('Applying TFCE to null data ... \n')
            H0_Condition_effect = load(H0filename);
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
                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Condition_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_Condition_effect,2),size(H0_Condition_effect,3),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_Condition_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_Condition_effect,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Condition_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            end
            save(H0_tfce_file,'tfce_H0_score','-v7.3');
            clear H0_Condition_effect tfce_H0_score;
        end
    end
    
    % interactions
    if contains(filename,'Interaction_effect_')
        Interaction_effect = load(filename);
        Interaction_effect = Interaction_effect.(cell2mat(fieldnames(Interaction_effect)));
        if size(Interaction_effect,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_score(1,:,:) = limo_tfce(2,squeeze(Interaction_effect(:,:,:,1)),[]);
            else
                tfce_score(1,:) = limo_tfce(1,squeeze(Interaction_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_score = limo_tfce(3,squeeze(Interaction_effect(:,:,:,1)),LIMO.data.neighbouring_matrix);
            else
                tfce_score = limo_tfce(2,squeeze(Interaction_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            end
        end
        save(tfce_file,'tfce_score','-v7.3');
        clear Interaction_effect tfce_score;
        
        H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
        if exist(H0filename,'file')
            fprintf('Applying TFCE to null data ... \n')
            H0_Interaction_effect = load(H0filename);
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
                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Interaction_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            else
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tfce_H0_score = NaN(size(H0_Interaction_effect,1),size(H0_Interaction_effect,2),size(H0_Interaction_effect,3),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_Interaction_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score = NaN(size(H0_Interaction_effect,1),size(H0_Interaction_effect,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            end
            save(H0_tfce_file,'tfce_H0_score','-v7.3');
            clear H0_Interaction_effect tfce_H0_score;
        end
    end
    
    if contains(filename,'con')
        con = load(filename);
        con = con.(cell2mat(fieldnames(con)));
        if strcmp(LIMO.Analysis ,'Time-Frequency')
            tfce_score = limo_tfce(3,squeeze(con(:,:,4)),LIMO.data.neighbouring_matrix);
        else
            tfce_score = limo_tfce(2,squeeze(con(:,:,4)),LIMO.data.neighbouring_matrix);
        end
        save(tfce_file,'tfce_score','-v7.3');
        
        H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
        if exist(H0filename,'file')
            fprintf('Applying TFCE to null data ... \n')
            H0_con = load(H0filename);
            H0_con = H0_con.(cell2mat(fieldnames(H0_con)));
            if strcmp(LIMO.Analysis ,'Time-Frequency')
                tfce_H0_score = limo_tfce(3,squeeze(H0_con(:,:,2,:)),LIMO.data.neighbouring_matrix);
            else
                tfce_H0_score = limo_tfce(2,squeeze(H0_con(:,:,2,:)),LIMO.data.neighbouring_matrix);
            end
            save(H0_tfce_file,'tfce_H0_score','-v7.3');
        end
    end
    
    if contains(filename,'ess')
        ess = load(filename);
        ess = ess.(cell2mat(fieldnames(ess)));
        if strcmp(LIMO.Analysis ,'Time-Frequency')
            tfce_score = limo_tfce(3,squeeze(ess(:,:,end-1)),LIMO.data.neighbouring_matrix);
        else
            tfce_score = limo_tfce(2,squeeze(ess(:,:,end-1)),LIMO.data.neighbouring_matrix);
        end
        save(tfce_file,'tfce_score','-v7.3');
        
        H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
        if exist(H0filename,'file')
            fprintf('Applying TFCE to null data ... \n')
            H0_ess = load(H0filename);
            H0_ess = H0_ess.(cell2mat(fieldnames(H0_ess)));
            if strcmp(LIMO.Analysis ,'Time-Frequency')
                tfce_H0_score = limo_tfce(3,squeeze(H0_ess(:,:,end-1,:)),LIMO.data.neighbouring_matrix);
            else
                tfce_H0_score = limo_tfce(2,squeeze(H0_ess(:,:,end-1,:)),LIMO.data.neighbouring_matrix);
            end
            save(H0_tfce_file,'tfce_H0_score','-v7.3');
        end
    end
    
    % -------------------------------------------------------------------------
else % 2nd level specific tests 
    % -------------------------------------------------------------------------
    
    % t-tests
    % ---------
    if contains(LIMO.design.name,'One sample','IgnoreCase',true) || ...
            contains(LIMO.design.name,'Two samples','IgnoreCase',true) || ...
            contains(LIMO.design.name,'Paired','IgnoreCase',true) 

        data = load(fullfile(filepath,filename));
        data = data.(cell2mat(fieldnames(data)));
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            if size(data,1) == 1
                tfce_data = limo_tfce(2,squeeze(data(:,:,:,4)),[]); % cluster in freq-time
            else
                tfce_data = limo_tfce(3,squeeze(data(:,:,:,4)),LIMO.data.neighbouring_matrix);
            end
        else
            if size(data,1) == 1
                tfce_data = limo_tfce(1,squeeze(data(:,:,4)),[]); % cluster in time or freq
            else
                tfce_data = limo_tfce(2,squeeze(data(:,:,4)),LIMO.data.neighbouring_matrix);
            end
        end
        
        if contains(LIMO.design.name,'One sample','IgnoreCase',true)
            tfce_one_sample = tfce_data;
            save(tfce_file, 'tfce_one_sample', '-v7.3');
            clear tfce_one_sample
        elseif contains(LIMO.design.name,'Two samples','IgnoreCase',true)
            tfce_two_samples = tfce_data;
            save(tfce_file, 'tfce_two_samples', '-v7.3');
            clear tfce_two_samples
        elseif contains(LIMO.design.name,'Paired','IgnoreCase',true)
            tfce_paired_samples = tfce_data;
            save(tfce_file, 'tfce_paired_samples', '-v7.3');
            clear tfce_paired_samples
        end
        clear data tfce_data ; disp('.. done');
        
        % do tfce for the data under H0
        H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
        if exist(H0filename,'file')
            fprintf('Applying TFCE to null data ... \n')
            H0_data = load(H0filename);
            H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
            if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                if size(H0_data,1) == 1
                    parfor b=1:LIMO.design.bootstrap
                        tfce_H0_data(:,:,b) = limo_tfce(2,squeeze(H0_data(1,:,:,1,b)),[],0);
                    end
                else
                    parfor b=1:LIMO.design.bootstrap
                        tfce_H0_data(:,:,:,b) = limo_tfce(3,squeeze(H0_data(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            else
                if size(H0_data,1) == 1
                    parfor b=1:LIMO.design.bootstrap
                        tfce_H0_data(:,:,b) = limo_tfce(1,squeeze(H0_data(1,:,1,b)),[],0);
                    end
                else
                    parfor b=1:LIMO.design.bootstrap
                        tfce_H0_data(:,:,b) = limo_tfce(2,squeeze(H0_data(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            end
            
            if contains(LIMO.design.name,'One sample','IgnoreCase',true)
                tfce_H0_one_sample = tfce_H0_data;
                save(H0_tfce_file, 'tfce_H0_one_sample', '-v7.3');
                clear tfce_H0_one_sample
            elseif contains(LIMO.design.name,'Two samples','IgnoreCase',true)
                tfce_H0_two_samples = tfce_H0_data;
                save(H0_tfce_file, 'tfce_H0_two_samples', '-v7.3');
                clear tfce_H0_two_samples
            elseif contains(LIMO.design.name,'Paired','IgnoreCase',true)
                tfce_H0_paired_samples = tfce_H0_data;
                save(H0_tfce_file, 'tfce_H0_paired_samples', '-v7.3');
                clear tfce_H0_paired_samples
            end
            clear tfce_H0_data; disp(' .. done')
        end
        
    % AN(C)OVA
    % ---------
    elseif contains(LIMO.design.name,'ANOVA','IgnoreCase',true) && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
            contains(LIMO.design.name,'ANCOVA','IgnoreCase',true)
   
        Condition = load(filename);
        Condition = Condition.(cell2mat(fieldnames(Condition)));
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            if size(Condition,1) == 1
                tfce_Condition = limo_tfce(2,squeeze(Condition(1,:,:,1)),[]); % use bwlabel to cluster freq*time map
            else
                tfce_Condition = limo_tfce(3,squeeze(Condition(:,:,:,1)),LIMO.data.neighbouring_matrix);
            end
        else
            if size(Condition,1) == 1
                tfce_Condition = limo_tfce(1,squeeze(Condition(:,:,1)),[]);
            else
                tfce_Condition = limo_tfce(2,squeeze(Condition(:,:,1)),LIMO.data.neighbouring_matrix);
            end
        end
        save(tfce_file, 'tfce_Condition', '-v7.3');
        clear tfce_Condition
        
        H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
        if exist(H0filename,'file')
            fprintf('Applying TFCE to null data ... \n')
            boot_H0_Condition = load(H0filename);
            boot_H0_Condition = boot_H0_Condition.(cell2mat(fieldnames(boot_H0_Condition)));
            if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                if size(boot_H0_Condition,1) == 1
                    tfce_H0_Condition = limo_tfce(2,squeeze(boot_H0_Condition(1,:,:,1,:)),[]);
                else
                    tfce_H0_Condition = NaN(size(boot_H0_Condition,1),size(boot_H0_Condition,2),size(boot_H0_Condition,3),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_Condition(:,:,:,b) = limo_tfce(3,squeeze(boot_H0_Condition(:,:,:,1,b)),LIMO.data.neighbouring_matrix);
                    end
                end
            else
                if size(boot_H0_Condition,1) == 1
                    tfce_H0_Condition = limo_tfce(1,squeeze(boot_H0_Condition(1,:,1,:)),[]);
                else
                    tfce_H0_Condition = NaN(size(boot_H0_Condition,1),size(boot_H0_Condition,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_Condition(:,:,b) = limo_tfce(2,squeeze(boot_H0_Condition(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                end
            end
            save(H0_tfce_file, 'tfce_H0_Condition', '-v7.3');
            clear tfce_H0_Condition
        end
    
    % Repeated Measure ANOVA
    % -----------------------
    elseif contains(LIMO.design.name,'Repeated measures ANOVA','IgnoreCase',true)
        
        if contains(filename,'Main_effect')
            % -----------------------------
            Rep_ANOVA = load(filename);
            Rep_ANOVA = Rep_ANOVA.(cell2mat(fieldnames(Rep_ANOVA)));
            if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                if size(Rep_ANOVA,1) == 1
                    tfce_Rep_ANOVA = limo_tfce(2,squeeze(Rep_ANOVA(1,:,:,1)),[]); % use bwlabel to cluster freq*time map
                else
                    tfce_Rep_ANOVA = limo_tfce(3,squeeze(Rep_ANOVA(:,:,:,1)),LIMO.data.neighbouring_matrix);
                end
            else
                if size(Rep_ANOVA,1) == 1
                    tfce_Rep_ANOVA = limo_tfce(1,squeeze(Rep_ANOVA(:,:,1)),[]);
                else
                    tfce_Rep_ANOVA = limo_tfce(2,squeeze(Rep_ANOVA(:,:,1)),LIMO.data.neighbouring_matrix);
                end
            end
            save(tfce_file, 'tfce_Rep_ANOVA', '-v7.3');
            clear tfce_Rep_ANOVA
            
            H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
            if exist(H0filename,'file')
                fprintf('Applying TFCE to null data ... \n')
                boot_H0_Rep_ANOVA = load(H0filename);
                boot_H0_Rep_ANOVA = boot_H0_Rep_ANOVA.(cell2mat(fieldnames(boot_H0_Rep_ANOVA)));
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    if size(boot_H0_Rep_ANOVA,1) == 1
                        tfce_H0_Rep_ANOVA = limo_tfce(2,squeeze(boot_H0_Rep_ANOVA(1,:,:,1,:)),[]);
                    else
                       tfce_H0_Rep_ANOVA = NaN(size(boot_H0_Rep_ANOVA,1),size(boot_H0_Rep_ANOVA,2),size(boot_H0_Rep_ANOVA,3),LIMO.design.bootstrap);
                       parfor b=1:nboot
                            tfce_H0_Rep_ANOVA(:,:,:,b) = limo_tfce(3,squeeze(boot_H0_Rep_ANOVA(:,:,:,1,b)),LIMO.data.neighbouring_matrix);
                        end
                    end
                else
                    if size(boot_H0_Rep_ANOVA,1) == 1
                        tfce_H0_Rep_ANOVA = limo_tfce(1,squeeze(boot_H0_Rep_ANOVA(1,:,1,:)),[]);
                    else
                        tfce_H0_Rep_ANOVA = NaN(size(boot_H0_Rep_ANOVA,1),size(boot_H0_Rep_ANOVA,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_Rep_ANOVA(:,:,b) = limo_tfce(2,squeeze(boot_H0_Rep_ANOVA(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                        end
                    end
                end
                save(H0_tfce_file, 'tfce_H0_Rep_ANOVA', '-v7.3');
                clear tfce_H0_Rep_ANOVA
            end
        
               
        elseif contains(filename,'ess')
            % -----------------------------
            ess = load(filename);
            ess = ess.(cell2mat(fieldnames(ess)));
            if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                if size(ess,1) == 1
                    tfce_ess = limo_tfce(2,squeeze(ess(1,:,:,4)),[]); % use bwlabel to cluster freq*time map
                else
                    tfce_ess = limo_tfce(3,squeeze(ess(:,:,:,4)),LIMO.data.neighbouring_matrix);
                end
            else
                if size(ess,1) == 1
                    tfce_ess = limo_tfce(1,squeeze(ess(:,:,4)),[]);
                else
                    tfce_ess = limo_tfce(2,squeeze(ess(:,:,4)),LIMO.data.neighbouring_matrix);
                end
            end
            save(tfce_file, 'tfce_ess', '-v7.3');
            clear tfce_ess
            
            H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
            if exist(H0filename,'file')
                fprintf('Applying TFCE to null data ... \n')
                H0_ess = load(H0filename);
                H0_ess = H0_ess.(cell2mat(fieldnames(H0_ess)));
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    if size(H0_ess,1) == 1
                        tfce_H0_ess = limo_tfce(2,squeeze(H0_ess(1,:,:,1,:)),[]);
                    else
                       tfce_H0_ess = NaN(size(H0_ess,1),size(H0_ess,2),size(H0_ess,3),LIMO.design.bootstrap);
                       parfor b=1:nboot
                            tfce_H0_ess(:,:,:,b) = limo_tfce(3,squeeze(H0_ess(:,:,:,1,b)),LIMO.data.neighbouring_matrix);
                        end
                    end
                else
                    if size(H0_ess,1) == 1
                        tfce_H0_ess = limo_tfce(1,squeeze(H0_ess(1,:,1,:)),[]);
                    else
                        tfce_H0_ess = NaN(size(H0_ess,1),size(H0_ess,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_ess(:,:,b) = limo_tfce(2,squeeze(H0_ess(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                        end
                    end
                end
                save(H0filename, 'tfce_H0_ess', '-v7.3');
                clear tfce_H0_ess
            end
        
        
        elseif contains(filename,'Gp') && ~contains(filename,'Interaction','IgnoreCase',true)
            Rep_ANOVA_Gp_effect = load(filename);
            Rep_ANOVA_Gp_effect = Rep_ANOVA_Gp_effect.(cell2mat(fieldnames(Rep_ANOVA_Gp_effect)));
            if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                if size(Rep_ANOVA_Gp_effect,1) == 1
                    tfce_Rep_ANOVA_Gp_effect = limo_tfce(2,limo_tf_4d_reshape(squeeze(1,Rep_ANOVA_Gp_effect(:,:,1))),[]);
                else
                    tfce_Rep_ANOVA_Gp_effect = limo_tfce(3,limo_tf_4d_reshape(squeeze(2,Rep_ANOVA_Gp_effect(:,:,1))),LIMO.data.neighbouring_matrix);
                end
            else
                if size(Rep_ANOVA_Gp_effect,1) == 1
                    tfce_Rep_ANOVA_Gp_effect = limo_tfce(1,squeeze(Rep_ANOVA_Gp_effect(:,:,1)),[]);
                else
                    tfce_Rep_ANOVA_Gp_effect = limo_tfce(2,squeeze(Rep_ANOVA_Gp_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                end
            end
            save(tfce_file, 'tfce_Rep_ANOVA_Gp_effect', '-v7.3');
            clear tfce_Rep_ANOVA_Gp_effect
            
            H0filename = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);
            if exist(H0filename,'file')
                fprintf('Applying TFCE to null data ... \n')
                H0_tfce_Rep_ANOVA_Gp_effect = load(H0filename);
                H0_tfce_Rep_ANOVA_Gp_effect = H0_tfce_Rep_ANOVA_Gp_effect.(cell2mat(fieldnames(H0_tfce_Rep_ANOVA_Gp_effect)));
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    if size(H0_tfce_Rep_ANOVA_Gp_effect,1) == 1
                        tfce_H0_Gp_effect = limo_tfce(2,squeeze(H0_tfce_Rep_ANOVA_Gp_effect(1,:,:,1,:)),[]);
                    else
                       tfce_H0_Gp_effect = NaN(size(H0_tfce_Rep_ANOVA_Gp_effect,1),size(H0_tfce_Rep_ANOVA_Gp_effect,2),size(H0_tfce_Rep_ANOVA_Gp_effect,3),LIMO.design.bootstrap);
                       parfor b=1:nboot
                            tfce_H0_Gp_effect(:,:,:,b) = limo_tfce(3,squeeze(H0_tfce_Rep_ANOVA_Gp_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix);
                        end
                    end
                else
                    if size(H0_tfce_Rep_ANOVA_Gp_effect,1) == 1
                        tfce_H0_Gp_effect = limo_tfce(1,squeeze(H0_tfce_Rep_ANOVA_Gp_effect(1,:,1,:)),[]);
                    else
                        tfce_H0_Gp_effect = NaN(size(H0_tfce_Rep_ANOVA_Gp_effect,1),size(H0_tfce_Rep_ANOVA_Gp_effect,2),LIMO.design.bootstrap);
                        parfor b=1:nboot
                            tfce_H0_Gp_effect(:,:,b) = limo_tfce(2,squeeze(H0_tfce_Rep_ANOVA_Gp_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                        end
                    end
                end
                save(H0filename, 'tfce_H0_Gp_effect', '-v7.3');
                clear tfce_H0_Gp_effect
            end
            
        elseif contains(filename,'Interaction','IgnoreCase',true) 
            for i=1:size(tmp_Rep_ANOVA_Interaction_with_gp,3)
                fprintf('analyzing interaction effect %g \n',i)
                tfce_name = sprintf('tfce%s%s',filesep,IRep_filenames{i});
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    if size(tmp_Rep_ANOVA_Interaction_with_gp,1) == 1
                        tfce_Rep_ANOVA_Interaction_with_gp = limo_tfce(2,limo_tf_4d_reshape(squeeze(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,1))),[]);
                    else
                        tfce_Rep_ANOVA_Interaction_with_gp = limo_tfce(3,limo_tf_4d_reshape(squeeze(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,1))),LIMO.data.neighbouring_matrix);
                    end
                else
                    if size(tmp_Rep_ANOVA_Interaction_with_gp,1) == 1
                        tfce_Rep_ANOVA_Interaction_with_gp = limo_tfce(1,squeeze(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,1)),[]);
                    else
                        tfce_Rep_ANOVA_Interaction_with_gp = limo_tfce(2,squeeze(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,1)),LIMO.data.neighbouring_matrix);
                    end
                end
                save(tfce_name, 'tfce_Rep_ANOVA_Interaction_with_gp');
                clear tfce_Rep_ANOVA_Interaction_with_gp
            end
        end
    end
end
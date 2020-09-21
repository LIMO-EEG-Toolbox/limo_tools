function limo_tfce_handling(varargin)

% routine for create tfce files commensurate to boostrapped files
%
% FORMAT  limo_tfce_handling(filename,'checkfile','yes')
%
% INPUTS filename is the stat file that need to be tfced (if H0 exist it is done too)
%        'checkfiles' is 'yes' by default - if 'no' anf tfce files already exist,
%                     it overwrites without asking overwise user is prompted
%
% OUTPUTS tfce_* files are saved on the drive in a tfce folder 
%         H0_tfce_* files are saved on the drive in the H0 folder
% ------------------------------------------------------------------
% Cyril R. Pernet - Copyright (C) LIMO Team 2020

%% Check inputs

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


%% quick user check
% ----------------

% files to create
tfce_file    = fullfile(LIMO.dir,['tfce' filesep 'tfce_' filename]);
H0_tfce_file = fullfile(LIMO.dir,['H0' filesep 'tfce_H0_' filename]);
% given filename input, we expect H0 to be
H0filename   = fullfile(LIMO.dir,['H0' filesep 'H0_' filename]);

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

if isfield(LIMO.design,'bootstrap')
    nboot = LIMO.design.bootstrap;
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


% create tfce folder
if ~exist(fullfile(LIMO.dir,'tfce'),'dir')
    mkdir(fullfile(LIMO.dir,'tfce'))
end
fprintf('Thresholding %s using TFCE \n',filename);

% -------------------------------------------------------------------------------
if contains(filename,'R2') || ...
        contains(filename,'semi_partial') % these files last dimension is R2, F, p
    % ----------------------------------------------------------------------------
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
    clear R2 tfce_score ;
    
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

    % -------------------------------------------------------------------------------------------------------------
elseif contains(filename,'con') || contains(filename,'ess') || ...
        contains(LIMO.design.name,'One sample','IgnoreCase',true) || ...
        contains(LIMO.design.name,'Two samples','IgnoreCase',true) || ...
        contains(LIMO.design.name,'Paired','IgnoreCase',true)  % these file last dimension is mean, se, df, t and p
    % --------------------------------------------------------------------------------------------------------------
    tval = load(filename);
    tval = tval.(cell2mat(fieldnames(tval)));
    if size(tval,1) == 1
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score(1,:,:) = limo_tfce(2, squeeze(tval(:,:,:,end-1)),[]); % no neighbouring, time-freq cluster
        else
            tfce_score(1,:)   = limo_tfce(1, squeeze(tval(:,:,end-1)),LIMO.data.neighbouring_matrix);
        end
    else
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score        = limo_tfce(3, squeeze(tval(:,:,:,end-1)),LIMO.data.neighbouring_matrix);
        else
            tfce_score        = limo_tfce(2, squeeze(tval(:,:,4)),LIMO.data.neighbouring_matrix);
        end
    end
    save(tfce_file,'tfce_score','-v7.3');
    clear tval tfce_score ;
    
    if exist(H0filename,'file')
        fprintf('Applying TFCE to null data ... \n')
        H0_tval = load(H0filename);
        H0_tval = H0_tval.(cell2mat(fieldnames(H0_tval)));
        if size(H0_tval,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(1,size(H0_tval,2),size(H0_tval,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,:,b) = limo_tfce(2,squeeze(H0_tval(:,:,:,end-1,b)),[],0);
                end
            else
                tfce_H0_score = NaN(1,size(H0_tval,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_tval(:,:,end-1,b)),LIMO.data.neighbouring_matrix,0);
                end
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(size(H0_tval,1),size(H0_tval,2),size(H0_tval,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,:,b) = limo_tfce(3,squeeze(H0_tval(:,:,:,end-1,b)),LIMO.data.neighbouring_matrix,0);
                end
            else
                tfce_H0_score = NaN(size(H0_tval,1),size(H0_tval,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_tval(:,:,end-1,b)),LIMO.data.neighbouring_matrix,0);
                end
            end
        end
        save(H0_tfce_file,'tfce_H0_score','-v7.3');
        clear H0_tval tfce_H0_score;
    end
    
    % ------------------------------------------
else % anything else last dimension is F and p
    % ------------------------------------------    
    Fval = load(filename);
    Fval = Fval.(cell2mat(fieldnames(Fval)));
    if size(Fval,1) == 1
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score(1,:,:) = limo_tfce(2,squeeze(Fval(:,:,:,1)),[]);
        else
            tfce_score(1,:) = limo_tfce(1,squeeze(Fval(:,:,1)),LIMO.data.neighbouring_matrix);
        end
    else
        if strcmpi(LIMO.Analysis,'Time-Frequency')
            tfce_score = limo_tfce(3,squeeze(Fval(:,:,:,1)),LIMO.data.neighbouring_matrix);
        else
            tfce_score = limo_tfce(2,squeeze(Fval(:,:,1)),LIMO.data.neighbouring_matrix);
        end
    end
    save(tfce_file,'tfce_score','-v7.3');
    clear Fval tfce_score;
    
    if exist(H0filename,'file')
        fprintf('Applying TFCE to null data ... \n')
        H0_Fval = load(H0filename);
        H0_Fval = H0_Fval.(cell2mat(fieldnames(H0_Fval)));
        if size(H0_Fval,1) == 1
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(1,size(H0_Fval,2),size(H0_Fval,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,:,b) = limo_tfce(1,squeeze(H0_Fval(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            else
                tfce_H0_score = NaN(1,size(H0_Fval,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Fval(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            end
        else
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                tfce_H0_score = NaN(size(H0_Fval,1),size(H0_Fval,2),size(H0_Fval,3),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,:,b) = limo_tfce(2,squeeze(H0_Fval(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            else
                tfce_H0_score = NaN(size(H0_Fval,1),size(H0_Fval,2),LIMO.design.bootstrap);
                parfor b=1:nboot
                    tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Fval(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                end
            end
        end
    end
    save(H0_tfce_file,'tfce_H0_score','-v7.3');
    clear H0_Fval tfce_H0_score   
end

function limo_batch_design_matrix(LIMOfile)

% this function wrap aound the LIMOfile structure from the batch and
% calls limo_design_matrix - of special interest is the ability to
% rearrange data per components based on the study structure clustering
% output - this allowing group statistics at the 2nd level. 
%
% Cyril Pernet and Ramon Martinez-Cancino, October 2014 updates for EEGLIMOLAB STUDY
% see also limo_batch 
%% -----------------------------
% Copyright (C) LIMO Team 2015

global EEGLIMO
load(LIMOfile);
if exist('EEGLIMO','var')
    if ~strcmp([LIMO.data.data_dir filesep LIMO.data.data],[EEGLIMO.filepath filesep EEGLIMO.filename])
        cd (LIMO.data.data_dir);
        disp('reloading data ..');
        EEGLIMO=pop_loadset([LIMO.data.data_dir filesep LIMO.data.data]);
    end
else
    disp('reloading data ..');
    EEGLIMO=pop_loadset([LIMO.data.data_dir filesep LIMO.data.data]);
end

if strcmp(LIMO.Analysis,'Time')
    if strcmp(LIMO.Type,'Components')
        if isfield(EEGLIMO.etc.datafiles,'icaerp')
            if ~iscell(EEGLIMO.etc.datafiles.icaerp) && strcmp(EEGLIMO.etc.datafiles.icaerp(end-3:end),'.mat')
                Y = load(EEGLIMO.etc.datafiles.icaerp);
                if isstruct(Y)
                    Y = getfield(Y,cell2mat(fieldnames(Y)));
                end
            else
                try
                    signal = load('-mat',EEGLIMO.etc.datafiles.icaerp);
                    if isstruct(signal); signal  = limo_struct2mat(signal); end
                catch
                    for d=1:length(EEGLIMO.etc.datafiles.icaerp)
                        signal{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.icaerp(d)));
                        if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
                    end
                end
                signal = limo_concatcells(signal);
            end
        else
            signal = eeg_getdatact(EEGLIMO,'component',[1:size(EEGLIMO.icawinv,2)]);
        end
        Y = signal(:,LIMO.data.trim1:LIMO.data.trim2,:); clear signal
        if isfield(LIMO.data,'cluster')
            try
                STUDY = evalin('base','STUDY');
            catch
                error('to run component clustering, you need the EEGLAB study loaded in the workspace with the clustering computed and saved')
            end
            nb_clusters = size(STUDY.cluster(1).child,2);
            nb_subjects = length({STUDY.datasetinfo.subject}); % length(unique({STUDY.datasetinfo.subject}));
            Cluster_matrix = parse_clustinfo(STUDY,STUDY.cluster(1).name);
            dsetinfo = rel2fullpath(STUDY.filepath,{STUDY.datasetinfo.filepath}');
            data_dir = rel2fullpath(STUDY.filepath,LIMO.data.data_dir(1:end));
            current_subject = find(cellfun(@strcmp, dsetinfo',repmat({data_dir},nb_subjects,1)));
            subject_name = {STUDY.datasetinfo(current_subject(1)).subject};
            newY = NaN(nb_clusters,size(Y,2),size(Y,3));
            for c=1:nb_clusters
                n = length(Cluster_matrix.clust(c).subj);
                tmp = find(cellfun(@strcmp,Cluster_matrix.clust(c).subj',repmat(subject_name,n,1)));
                if ~isempty(tmp)
                    which_ics = unique(Cluster_matrix.clust(c).ics(tmp));
                    if length(which_ics)==1
                        newY(c,:,:) =  Y(which_ics,:,:);
                    else
                        newY(c,:,:) = limo_combine_components(Y,EEGLIMO.icaweights,EEGLIMO.icawinv,which_ics);
                    end
                end
            end
            Y = newY; clear newY;
        end
    else % channels
        try isfield(EEGLIMO.etc, 'datafiles') && isfield(EEGLIMO.etc.datafiles,'daterp')
            if ~iscell(EEGLIMO.etc.datafiles.daterp) && strcmp(EEGLIMO.etc.datafiles.daterp(end-3:end),'.mat')
                Y = load(EEGLIMO.etc.datafiles.daterp);
                if isstruct(Y)
                    Y = getfield(Y,cell2mat(fieldnames(Y)));
                end
            else
                for d=1:length(EEGLIMO.etc.datafiles.daterp)
                    Y{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.daterp(d)));
                    if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                end
                Y = limo_concatcells(Y);
            end
        catch 
            disp('the field EEG.etc.datafiles.daterp pointing to the data is missing - using EEGLIMO.data')
            Y = EEGLIMO.data; %Y = EEGLIMO.data(:,LIMO.data.trim1:LIMO.data.trim2,:); 
        end
        Y = Y(:,LIMO.data.trim1:LIMO.data.trim2,:); 
        clear EEGLIMO
    end
    
elseif strcmp(LIMO.Analysis,'Frequency')
    
    if strcmp(LIMO.Type,'Components')
        if isfield(EEGLIMO.etc.datafiles,'icaspec')
            if ~iscell(EEGLIMO.etc.datafiles.icaspec) && strcmp(EEGLIMO.etc.datafiles.icaspec(end-3:end),'.mat')
                Y = load(EEGLIMO.etc.datafiles.icaspec);
                if isstruct(Y)
                    Y = getfield(Y,cell2mat(fieldnames(Y)));
                end
            else
                try
                    signal = load('-mat',EEGLIMO.etc.datafiles.icaspec);
                    if isstruct(signal); signal  = limo_struct2mat(signal); end
                catch
                    for d=1:length(EEGLIMO.etc.datafiles.icaspec)
                        signal{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.icaspec(d)));
                        if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
                    end
                end
                signal = limo_concatcells(signal);
            end
        else
            signal = eeg_getdatact(EEGLIMO,'component',[1:length(EEGLIMO.icawinv)]);
        end
        Y = signal(:,LIMO.data.trim1:LIMO.data.trim2,:); clear signal
        if isfield(LIMO.data,'cluster')
            try
                STUDY = evalin('base','STUDY');
            catch
                error('to run component clustering, you need the EEGLAB study loaded in the workspace with the clustering computed and saved')
            end
            nb_clusters = size(STUDY.cluster(1).child,2);
            nb_subjects = length({STUDY.datasetinfo.subject}); % length(unique({STUDY.datasetinfo.subject}));
            Cluster_matrix = parse_clustinfo(STUDY,STUDY.cluster(1).name);
            dsetinfo = rel2fullpath(STUDY.filepath,{STUDY.datasetinfo.filepath}');
            data_dir = rel2fullpath(STUDY.filepath,LIMO.data.data_dir(1:end));
            current_subject = find(cellfun(@strcmp, dsetinfo',repmat({data_dir},nb_subjects,1)));
            subject_name = {STUDY.datasetinfo(current_subject(1)).subject};
            newY = NaN(nb_clusters,size(Y,2),size(Y,3));
            for c=1:nb_clusters
                n = length(Cluster_matrix.clust(c).subj);
                tmp = find(cellfun(@strcmp,Cluster_matrix.clust(c).subj',repmat(subject_name,n,1)));
                if ~isempty(tmp)
                    which_ics = unique(Cluster_matrix.clust(c).ics(tmp));
                    if length(which_ics)==1
                        newY(c,:,:) =  Y(which_ics,:,:);
                    else
                        newY(c,:,:) = limo_combine_components(Y,EEGLIMO.icaweights,EEGLIMO.icawinv,which_ics);
                    end
                end
            end
            Y = newY; clear newY;
        end
    else % channels
        if isfield(EEGLIMO.etc, 'datafiles') && isfield(EEGLIMO.etc.datafiles,'datspec')
            if ~iscell(EEGLIMO.etc.datafiles.datspec) && strcmp(EEGLIMO.etc.datafiles.datspec(end-3:end),'.mat')
                Y = load(EEGLIMO.etc.datafiles.datspec);
                if isstruct(Y)
                    Y = getfield(Y,cell2mat(fieldnames(Y)));
                end
            else
                for d=1:length(EEGLIMO.etc.datafiles.datspec)
                    Y{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.datspec(d)));
                    if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
                end
                Y = limo_concatcells(Y); clear EEGLIMO
            end
        else
            error('the field EEG.etc.datspec pointing to the data is missing')
        end
        Y = Y(:,LIMO.data.trim1:LIMO.data.trim2,:);
    end
    
elseif strcmp(LIMO.Analysis,'Time-Frequency')
    disp('Time-Frequency implementation - loading tf data, be patient ...');
    
    if strcmp(LIMO.Type,'Components')
        if ~iscell(EEGLIMO.etc.datafiles.datspec) && isfield(EEGLIMO.etc.datafiles,'icatimef')
            try
                signal = load('-mat',EEGLIMO.etc.datafiles.icatimef);
                if isstruct(signal); signal  = limo_struct2mat(signal); end
            catch
                for d=1:length(EEGLIMO.etc.datafiles.icaspec)
                    signal{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.icatimef(d)));
                    if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
                end
            end
            signal = limo_concatcells(signal);
        else
            signal = eeg_getdatact(EEGLIMO,'component',[1:length(EEGLIMO.icawinv)]);
        end
        Y = signal(:,LIMO.data.trim_low_f:LIMO.data.trim_high_f,LIMO.data.trim1:LIMO.data.trim2,:); clear signal
        if isfield(LIMO.data,'cluster')
            try
                STUDY = evalin('base','STUDY');
            catch
                error('to run component clustering, you need the EEGLAB study loaded in the workspace with the clustering computed and saved')
            end
            nb_clusters = size(STUDY.cluster(1).child,2);
            nb_subjects = length(unique({STUDY.datasetinfo.subject}));
            Cluster_matrix = parse_clustinfo(STUDY,STUDY.cluster(1).name);
            dsetinfo = rel2fullpath(STUDY.filepath,{STUDY.datasetinfo.filepath}');
            data_dir = rel2fullpath(STUDY.filepath,LIMO.data.data_dir(1:end));
            current_subject = find(cellfun(@strcmp, dsetinfo',repmat({data_dir},nb_subjects,1)));
            subject_name = {STUDY.datasetinfo(current_subject).subject};
            newY = NaN(nb_clusters,size(Y,2),size(Y,3));
            for c=1:nb_clusters
                n = length(Cluster_matrix.clust(c).subj);
                tmp = find(cellfun(@strcmp,Cluster_matrix.clust(c).subj',repmat(subject_name,n,1)));
                if ~isempty(tmp)
                    which_ics = unique(Cluster_matrix.clust(c).ics(tmp));
                    if length(which_ics)==1
                        newY(c,:,:,:) =  Y(which_ics,:,:);
                    else
                        newY(c,:,:,:) = limo_combine_components(Y,EEGLIMO.icaweights,EEGLIMO.icawinv,which_ics);
                    end
                end
            end
            Y = newY; clear newY;
        end
    else % channels
        if isfield(EEGLIMO.etc, 'datafiles') && isfield(EEGLIMO.etc.datafiles,'dattimef')
            for d=1:length(EEGLIMO.etc.datafiles.dattimef)
                Y{d} = load('-mat',cell2mat(EEGLIMO.etc.datafiles.dattimef(d)));
                if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
            end
            Y = limo_concatcells(Y);
        elseif isfield(EEGLIMO.etc.datafiles,'datersp')
            [~,~,ext]=fileparts(EEGLIMO.etc.datafiles.datersp);
            if strcmp(ext,'.dattimef')
                [Y,TimeVect,FreqVect] = limo_struct2mat(EEGLIMO.etc.datafiles.datersp);
            else
                Y = load(EEGLIMO.etc.datafiles.datersp);
                if isstruct(Y)
                    Y = getfield(Y,cell2mat(fieldnames(Y)));
                end
            end
        else
            error('no data found, the field EEG.etc.dattimef or EEGLIMO.etc.datersp pointing to the data is missing')
        end
        Y = Y(:,LIMO.data.trim_low_f:LIMO.data.trim_high_f,LIMO.data.trim1:LIMO.data.trim2,:);
    end
    clear EEGLIMO
    LIMO.data.size4D= size(Y);
    LIMO.data.size3D= [LIMO.data.size4D(1) LIMO.data.size4D(2)*LIMO.data.size4D(3) LIMO.data.size4D(4)];
end

clear ALLEEGLIMO
cd (LIMO.dir) ; save LIMO LIMO

% make the design matrix
disp('computing design matrix');
if strcmp(LIMO.Analysis,'Time-Frequency') % use limo_design_matrix_tf
    [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
        LIMO.design.nb_continuous] = limo_design_matrix_tf(Y, LIMO,0);
else  % for time or power use limo_design_matrix
    [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
        LIMO.design.nb_continuous] = limo_design_matrix(Y, LIMO,0);
end

% update LIMO.mat
LIMO.design.name  = 'batch processing';
LIMO.design.status = 'to do';
save LIMO LIMO; clear Y

end
% -------------------------------------------------------------------------
function file_fullpath = rel2fullpath(studypath,filepath)
% Return full path if 'filepath' is a relative path. The output format will
% fit the one of 'filepath'. That means that if 'filepath' is a cell array,
% then the output will a cell array too, and the same if is a string.

nit = 1; if iscell(filepath), nit = length(filepath);end

for i = 1:nit
    if iscell(filepath),pathtmp = filepath{i}; else pathtmp = filepath; end
    if strfind(pathtmp(end),filesep), pathtmp = pathtmp(1:end-1); end % Getting rid of filesep at the end
    if strfind(pathtmp(1:2),['.' filesep])
        if iscell(filepath),
            file_fullpath{i} = fullfile(studypath,pathtmp(3:end));
        else
            file_fullpath = fullfile(studypath,pathtmp(3:end));
        end
    else
        if iscell(filepath),
            file_fullpath{i} = pathtmp;
        else
            file_fullpath = pathtmp;
        end
    end
end
end
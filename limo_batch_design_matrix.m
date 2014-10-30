function limo_batch_design_matrix(LIMOfile)
global EEG

load(LIMOfile);
if exist('EEG','var')
    if ~strcmp([LIMO.data.data_dir filesep LIMO.data.data],[EEG.filepath filesep EEG.filename])
        cd (LIMO.data.data_dir);
        disp('reloading data ..');
        EEG=pop_loadset(LIMO.data.data);
    end
else
    disp('reloading data ..');
    EEG=pop_loadset(LIMO.data.data);
end

if strcmp(LIMO.Analysis,'Time')
    if strcmp(LIMO.Type,'Components')
        if isfield(EEG.etc.datafiles,'icaerp')
            for d=1:length(EEG.etc.datafiles.icaerp)
                signal{d} = load('-mat',cell2mat(EEG.etc.datafiles.icaerp(d)));
                if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
            end
            signal = limo_concatcells(signal);
        else
            signal = eeg_getdatact(EEG,'component',[1:length(EEG.icawinv)]);
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
            current_subject = find(cellfun(@strcmp, {STUDY.datasetinfo.filepath}',repmat({LIMO.data.data_dir(1:end)},nb_subjects,1)));
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
                        newY(c,:,:) = limo_combine_components(Y,EEG.icaweights,EEG.icawinv,which_ics);
                    end
                end
            end
            Y = newY; clear newY;
        end
    else % channels
        if isfield(EEG.etc.datafiles,'daterp')
            for d=1:length(EEG.etc.datafiles.daterp)
                Y{d} = load('-mat',cell2mat(EEG.etc.datafiles.daterp(d)));
                if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
            end
            Y = limo_concatcells(Y);
        else
            disp('the field EEG.etc.datafiles.daterp pointing to the data is missing - using EEG.data')
            Y = EEG.data(:,LIMO.data.trim1:LIMO.data.trim2,:); 
        end
        clear EEG
    end
    
elseif strcmp(LIMO.Analysis,'Frequency')
    
    if strcmp(LIMO.Type,'Components')
        if isfield(EEG.etc.datafiles,'icaspec')
            for d=1:length(EEG.etc.datafiles.icaspec)
                signal{d} = load('-mat',cell2mat(EEG.etc.datafiles.icaspec(d)));
                if isstruct(signal{d}); signal{d}  = limo_struct2mat(signal{d}); end
            end
            signal = limo_concatcells(signal);
        else
            signal = eeg_getdatact(EEG,'component',[1:length(EEG.icawinv)]);
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
            current_subject = find(cellfun(@strcmp, {STUDY.datasetinfo.filepath}',repmat({LIMO.data.data_dir(1:end)},nb_subjects,1)));
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
                        newY(c,:,:) = limo_combine_components(Y,EEG.icaweights,EEG.icawinv,which_ics);
                    end
                end
            end
            Y = newY; clear newY;
        end
    else % channels
        if isfield(EEG.etc.datafiles,'datspec')
            for d=1:length(EEG.etc.datafiles.datspec)
                Y{d} = load('-mat',cell2mat(EEG.etc.datafiles.datspec(d)));
                if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
            end
            Y = limo_concatcells(Y); clear EEG
        else
            error('the field EEG.etc.datspec pointing to the data is missing')
        end
    end
    
elseif strcmp(LIMO.Analysis,'Time-Frequency')
    disp('Time-Frequency implementation - loading tf data...');
    
    if strcmp(LIMO.Type,'Components')
        if isfield(EEG.etc.datafiles,'icatimef')
            for d=1:length(EEG.etc.datafiles.icatimef)
                signal{d} = load('-mat',cell2mat(EEG.etc.datafiles.icatimef(d)));
                if isstruct(signal{d}); signal{d} = limo_struct2mat(signal{d}); end
            end
            signal = limo_concatcells(signal);
        else
            signal = eeg_getdatact(EEG,'component',[1:length(EEG.icawinv)]);
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
            current_subject = find(cellfun(@strcmp, {STUDY.datasetinfo.filepath}',repmat({LIMO.data.data_dir(1:end)},nb_subjects,1)));
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
                        newY(c,:,:,:) = limo_combine_components(Y,EEG.icaweights,EEG.icawinv,which_ics);
                    end
                end
            end
            Y = newY; clear newY;
        end
    else % channels
        if isfield(EEG.etc.datafiles,'dattimef')
            for d=1:length(EEG.etc.datafiles.dattimef)
                Y{d} = load('-mat',cell2mat(EEG.etc.datafiles.dattimef(d)));
                if isstruct(Y{d}); Y{d}  = limo_struct2mat(Y{d}); end
            end
            Y = limo_concatcells(Y);
            clear EEG
        else
            error('the field EEG.etc.dattimef pointing to the data is missing')
        end
    end
    
    LIMO.data.size4D= size(Y);
    LIMO.data.size3D= [LIMO.data.size4D(1) LIMO.data.size4D(2)*LIMO.data.size4D(3) LIMO.data.size4D(4)];
end

clear ALLCOM ALLEEG CURRENTSET CURRENTSTUDY LASTCOM STUDY
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



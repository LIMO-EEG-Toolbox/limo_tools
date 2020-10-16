%% Call LIMO functions from FieldTrip files

%% (DATA SPECIFIC) Set the paths and names
PATH_TO_RAW_EEG = 'D:\__EEG-data';
PATH_TO_PROCESSED_EEG = 'D:\__EEG-data\EEG_Erika_format\EEG';
PATH_TO_TEMPLATE_ELEC = fullfile('D:\FieldTrip\template\electrode','GSN-HydroCel-64_1.0.sfp');
PATH_TO_TEMPLATE_NEIGHBOURS = 'D:\FieldTrip\template\neighbours\biosemi64_neighb.mat';
PATH_TO_LIMO = 'C:\Users\luca-\OneDrive - UMONS\_PhD\_Matlab\3-Dynamic-Analysis\limo_tools';
PATH_TO_ROOT = 'D:\__EEG-data\BIDS_files';
PATH_TO_DERIV = 'D:\__EEG-data\BIDS_files\derivatives';
raw_eeg_common_name = 'ARC_J*';
processed_eeg_common_name = 'J_*';
processed_eeg_mat = 'clean_eeg.mat';
task_name = 'semantic-priming';

%% (DATA SPECIFIC) Create derivatives files and set model.cat_files
template_elec  = ft_read_sens(PATH_TO_TEMPLATE_ELEC);
e1 = find(strcmp(template_elec.label,'E1'));

dinfo = dir(fullfile(PATH_TO_PROCESSED_EEG,processed_eeg_common_name));
subj = {dinfo.name};

model.cat_files = {};
subj_ID = 1;
for subj_name = drange(subj)
    processed_eeg = load(fullfile(PATH_TO_PROCESSED_EEG,subj_name{1},processed_eeg_mat));
    processed_eeg = struct2cell(processed_eeg);
    processed_eeg = processed_eeg{1};
    deriv_mat  = processed_eeg{1};
    
    model.cat_files{subj_ID,1} = ones(length(deriv_mat.trial),1);
    for i = 2:length(processed_eeg)
        tmp = processed_eeg{i};
        for j = 1:length(tmp.trial)
            deriv_mat.trial{end+1} = tmp.trial{j};
            deriv_mat.time{end+1} = tmp.time{j};
            deriv_mat.trialinfo(end+1) = tmp.trialinfo(j);
            deriv_mat.sampleinfo(end+1,:) = tmp.sampleinfo(j,:);
        end
        model.cat_files{subj_ID,1} = [model.cat_files{subj_ID,1}; i*ones(length(tmp.trial),1)];
    end
    
    for i = e1:length(template_elec.chanpos)
        deriv_mat.chanlocs(i-e1+1).labels = template_elec.label(i);
        deriv_mat.chanlocs(i-e1+1).X = template_elec.chanpos(i,1);
        deriv_mat.chanlocs(i-e1+1).Y = template_elec.chanpos(i,2);
        deriv_mat.chanlocs(i-e1+1).Z = template_elec.chanpos(i,3);
    end
    
    if subj_ID >= 10
        subfolder = 'sub-0%d';
    else
        subfolder = 'sub-00%d';
    end
    derivatives_path = fullfile(PATH_TO_DERIV,sprintf(subfolder,subj_ID),'eeg', [sprintf(subfolder,subj_ID) '_task-' task_name '_eeg.mat']);
    [root,name,ext] = fileparts(derivatives_path);
    if ~exist(root,'dir')
        mkdir(root)
    end
    save(derivatives_path,'deriv_mat');
    
    subj_ID = subj_ID + 1;
end

%% (NEED EEG_JSON AND CHANNELS.TSV FILES FOR SUB-001 TO BE MANUALLY CREATED) Create BIDS files
json_path = fullfile(PATH_TO_ROOT,sprintf('sub-00%d',1),'eeg',[sprintf('sub-00%d',1) '_task-' task_name '_eeg.json']);
channel_path = fullfile(PATH_TO_ROOT,sprintf('sub-00%d',1),'eeg', [sprintf('sub-00%d',1) '_task-' task_name '_channels.tsv']);
if ~exist(json_path,'file') || ~exist(channel_path,'file')
    disp('ERROR ! NEED EEG_JSON AND CHANNELS.TSV FILES FOR SUB-001 TO BE MANUALLY CREATED')
    return;
end

dinfo = dir(fullfile(PATH_TO_RAW_EEG,raw_eeg_common_name));
fileNames = { dinfo.name };
for iFile = 1 : numel( dinfo )
    if iFile >= 10
        subfolder = 'sub-0%d';
    else
        subfolder = 'sub-00%d';
    end
    output_path = fullfile(PATH_TO_ROOT,sprintf(subfolder,iFile),'eeg', [sprintf(subfolder,iFile) '_task-' task_name '_eeg.bdf']);
    [root,name,ext] = fileparts(output_path);
    if ~exist(root,'dir')
        mkdir(root)
    end
    copyfile(fullfile(PATH_TO_RAW_EEG,dinfo(iFile).name),output_path);

    % create json & channels file
    if iFile ~= 1
        copyfile(json_path,[output_path(1:end-3) 'json']);
        copyfile(channel_path,[output_path(1:end-7) 'channels.tsv']);
    end
    
    % create events file
    event = ft_read_event(fullfile(PATH_TO_RAW_EEG,dinfo(iFile).name));
    type = {event.type}';
    sample = {event.sample}';
    value = {event.value}';
    offset = {event.offset}';
    duration = {event.duration}';
    T = table(type,sample,value,offset,duration);
    event_path = [output_path(1:end-7) 'events.txt'];
    writetable(T,event_path,'Delimiter','\t');
    movefile(event_path, [event_path(1:end-3) 'tsv'])
end

%% (DATA SPECIFIC) Call limo_batch function
option = 'both';

% model is a structure that specifiy information to build a model
dinfo = dir(fullfile(PATH_TO_DERIV,'sub-*'));
model.set_files = {};
for i = 1 : numel( dinfo )
    if i >= 10
        subfolder = 'sub-0%d';
    else
        subfolder = 'sub-00%d';
    end
    model.set_files{i,1} = fullfile(PATH_TO_DERIV,sprintf(subfolder,i),'eeg',[sprintf(subfolder,i) '_task-' task_name '_eeg.mat']);
end
model.cont_files = {};
% model.defaults: specifiy the parameters to use for each subject
model.defaults.type = 'Channels'; %or 'Components'
model.defaults.analysis = 'Time'; %'Frequency' or 'Time-Frequency'
model.defaults.method = 'OLS'; %'IRLS' 'WLS'
model.defaults.type_of_analysis = 'Mass-univariate'; %or 'Multivariate'
model.defaults.fullfactorial = 0; %or 1
model.defaults.zscore = 0; %or 1
model.defaults.start = -0.2; %starting time in ms
model.defaults.end = 0.5; %ending time in ms
model.defaults.bootstrap = 0; %or 1
model.defaults.tfce = 0; %or 1
neighbours = load(PATH_TO_TEMPLATE_NEIGHBOURS);
neighbours = neighbours.neighbours;
neighb_mat = zeros(64,64);
for i = 1:64
    idx = [];
    for neighlab = convertCharsToStrings((neighbours(i).neighblabel)')
        tmp = find(strcmp(convertCharsToStrings({neighbours.label}),neighlab{1}));
        idx = [idx tmp];
    end
    neighb_mat(i,idx) = 1;
end
model.defaults.neighbouring_matrix = neighb_mat; %neighbouring matrix use for clustering (necessary if bootstrap = 1)
%neighbouring matrix format: [n_chan x n_chan] of 1/0 (neighbour or not)

contrast.mat = [1 0 0 0 1 0 0 0 0;
                0 1 0 0 0 1 0 0 0;
                0 0 1 0 1 0 0 0 0;
                0 0 0 1 0 1 0 0 0];

save(fullfile(PATH_TO_DERIV,'model.mat'),'model')
save(fullfile(PATH_TO_DERIV,'contrast.mat'),'contrast')
% 
% % %uncomment if you want to load an existing model/contrast
% % model = load(fullfile(PATH_TO_DERIV,'model.mat'));
% % model = model.model;
% % contrast = load(fullfile(PATH_TO_DERIV,'contrast.mat'));
% % contrast = contrast.contrast;
% 
[LIMO_files, procstatus] = limo_batch(option, model,contrast);

%% call limo_random_select
LIMOfiles = fullfile(LIMO_files.LIMO,'Beta_files_GLM_OLS_Time_Channels.txt');
for i = e1:length(template_elec.chanpos)
    expected_chanlocs.expected_chanlocs(i-e1+1).labels = template_elec.label(i);
    expected_chanlocs.expected_chanlocs(i-e1+1).X = template_elec.chanpos(i,1);
    expected_chanlocs.expected_chanlocs(i-e1+1).Y = template_elec.chanpos(i,2);
    expected_chanlocs.expected_chanlocs(i-e1+1).Z = template_elec.chanpos(i,3);
end
expected_chanlocs.channeighbstructmat = model.defaults.neighbouring_matrix;

LIMOPath = limo_random_select('Repeated Measures ANOVA',expected_chanlocs,'LIMOfiles',...
    LIMOfiles,'analysis_type','Full scalp analysis','parameters',{[1 2];[3 4]},...
    'factor names',{'relation'},'type','Channels','nboot',10,'tfce',1,'skip design check','yes');

function limo_batch_import_data(EEG_DATA,cat,cont,defaults)

% routine to import 
%
% FORMAT limo_batch_import_data(setfile,cat,cont,defaults)
% 
% INPUT setfile is a an EEG files .set to be loaded
%       cat and cont are either numeric or txt or mat files
%               corresponding to the regressors in the model
%       defaults is a structure specifying all the parameters
%                to use in the GLM (ie set in LIMO.mat)
% 
% OUTPUT create a LIMO.mat with the relevant info in the subject
%        directory specified in the defaults -- importantly if 
%        some info are not in the default, it tries to read it from
%        the EEG.set file, from EEG.etc
%
% see also limo_batch 
% ------------------------------
%  Copyright (C) LIMO Team 2021

global EEGLIMO

EEGLIMO                      = load('-mat',EEG_DATA);
EEGLIMO                      = EEGLIMO.(cell2mat(fieldnames(EEGLIMO)));
if ~isfield(EEGLIMO,'filepath')
    [root,name,ext]              = fileparts(EEG_DATA);
    EEGLIMO.filepath             = root;
    EEGLIMO.filename             = [name ext];
end
LIMO.dir                     = defaults.name;
LIMO.data.data               = [name ext];
LIMO.data.data_dir           = root;
if strcmp(ext,'.set') %EEGLAB
    LIMO.data.sampling_rate      = EEGLIMO.srate;
elseif strcmp(ft_datatype(EEGLIMO),'raw') %FieldTrip
    LIMO.data.sampling_rate      = EEGLIMO.fsample;
    if ~isfield(EEGLIMO,'elec') || (isfield(EEGLIMO,'elec') && ~isfield(EEGLIMO,'chanlocs'))
        EEGLIMO = limo_get_ft_chanlocs(EEGLIMO, defaults);
    end
elseif strcmp(ft_datatype(EEGLIMO),'source') %FieldTrip source
    if ~isfield(EEGLIMO,'chanlocs')
        EEGLIMO = limo_get_ft_chanlocs(EEGLIMO, defaults);
    end
    LIMO.data.sampling_rate = length(EEGLIMO.time)/(EEGLIMO.time(end)-EEGLIMO.time(1));
    % adapt the time field
    tmp = EEGLIMO.time;
    clear EEGLIMO.time
    EEGLIMO.time = {};
    EEGLIMO.time{1} = tmp;
else
    error('ERROR in limo_batch_import_data: neither EEGLAB nor FieldTrip data')
end

LIMO.Analysis                = defaults.analysis;
LIMO.Type                    = defaults.type;
LIMO.design.zscore           = defaults.zscore;
LIMO.design.method           = defaults.method;
LIMO.design.type_of_analysis = defaults.type_of_analysis;
LIMO.design.fullfactorial    = defaults.fullfactorial;
LIMO.design.bootstrap        = defaults.bootstrap;
LIMO.design.tfce             = defaults.tfce;
LIMO.design.status           = 'to do';
LIMO.Level                   = 1;

% optional fields for EEGLAB study
if isfield(defaults,'icaclustering')
    LIMO.data.cluster = defaults.icaclustering;
end

if isfield(defaults,'chanlocs')
    LIMO.data.chanlocs = defaults.chanlocs;
else
    LIMO.data.chanlocs = EEGLIMO.chanlocs;
end

if isfield(defaults,'neighbouring_matrix')
    LIMO.data.neighbouring_matrix = defaults.neighbouring_matrix;
end

if isfield(defaults,'studyinfo')
    LIMO.data.studyinfo = defaults.studyinfo; % same as STUDY.design(design_index).variable;
end

% update according to the type of data
if strcmp(defaults.analysis,'Time') 
  
    if isfield(EEGLIMO,'etc') && isfield(EEGLIMO.etc,'timeerp') %EEGLAB
        timevect = EEGLIMO.etc.timeerp;
    elseif isfield(EEGLIMO,'time') % FieldTrip
        timevect = EEGLIMO.time{1}*1000; %convert in ms
    else
        warning('the field EEG.etc.timeerp is missing');
        if isfield(EEGLIMO,'times')
            timevect = EEGLIMO.times;
        end
    end
    
    % start
    if isempty(defaults.start) || defaults.start < min(timevect)
        LIMO.data.start = timevect(1);
        LIMO.data.trim1 = 1;
    else
        [~,position]    = min(abs(timevect - defaults.start));
        LIMO.data.start = timevect(position);
        LIMO.data.trim1 = position;
    end
    
    % end
    if isempty(defaults.end) || defaults.end > max(timevect)
        LIMO.data.end   = timevect(end);
        LIMO.data.trim2 = length(timevect);
    else
        [~,position]    = min(abs(timevect - defaults.end));
        LIMO.data.end   = timevect(position);
        LIMO.data.trim2 = position;
    end
    
    LIMO.data.timevect  = timevect(LIMO.data.trim1:LIMO.data.trim2);
    
% elseif strcmp(defaults.analysis,'Frequency') 
%     
%     if ~isfield(EEGLIMO.etc,'freqspec')
%         disp('the fied EEG.etc.freqspec is missing - reloading single trials');
%         data     = load('-mat',EEGLIMO.etc.freqspec);
%         freqvect = data.freqs; clear data;
%     else
%         freqvect = EEGLIMO.etc.freqspec;
%     end
% 
%     % start
%     if isempty(defaults.lowf) || defaults.lowf < freqvect(1)
%         LIMO.data.start = freqvect(1);
%         LIMO.data.trim1 = 1;
%     else
%         [~,position]    = min(abs(freqvect-defaults.lowf));
%         LIMO.data.start = freqvect(position);
%         LIMO.data.trim1 = position; 
%     end
%     
%     % end
%     if isempty(defaults.highf) || defaults.highf > freqvect(end)
%         LIMO.data.end   = freqvect(end);
%         LIMO.data.trim2 = numel(freqvect);
%     else
%         [~,position]    = min(abs(freqvect-defaults.highf));
%         LIMO.data.end   = freqvect(position);
%         LIMO.data.trim2 = position; 
%     end
%     
%     LIMO.data.freqlist  = freqvect(LIMO.data.trim1:LIMO.data.trim2);
% 
% elseif strcmp(defaults.analysis,'Time-Frequency')
%     
%     if ~isfield(EEGLIMO.etc,'timeersp') || ~isfield(EEGLIMO.etc,'freqersp')
%         disp('ersp fied in EEG.etc absent or impcomplete, reloading the single trials')
%         data = load('-mat',EEGLIMO.etc.timef,'times','freqs');
%         timevect = data.times;
%         freqvect = data.freqs;
%     else
%         timevect = EEGLIMO.etc.timeersp;
%         freqvect = EEGLIMO.etc.freqersp;
%     end
%        
%     % start
%     if isempty(defaults.start) || defaults.start < min(timevect)
%         LIMO.data.start = timevect(1);
%         LIMO.data.trim1 = 1;    
%     else
%         [~,position]    = min(abs(timevect - defaults.start));
%         LIMO.data.start = timevect(position);
%         LIMO.data.trim1 =  find(timevect == LIMO.data.start);
%     end
%     
%     % end
%     if isempty(defaults.end) || defaults.end > max(timevect)
%         LIMO.data.end   = timevect(end);
%         LIMO.data.trim2 = length(timevect);    
%     else
%         [~,position]    = min(abs(timevect - defaults.end));
%         LIMO.data.end   = timevect(position);
%         LIMO.data.trim2 =  position;
%     end
% 
%     LIMO.data.tf_times  = timevect(LIMO.data.trim1:LIMO.data.trim2);
% 
%     % start
%     if isempty(defaults.lowf) || defaults.lowf < freqvect(1)
%         LIMO.data.lowf = freqvect(1);
%         LIMO.data.trim_lowf = 1;
%     else
%         [~,position] = min(abs(freqvect-defaults.lowf));
%         LIMO.data.lowf = freqvect(position);
%         LIMO.data.trim_lowf = position; 
%     end
%     
%     % end
%     if isempty(defaults.highf) || defaults.highf > freqvect(end)
%         LIMO.data.highf = freqvect(end);
%         LIMO.data.trim_highf = length(freqvect);
%     else
%         [~,position] = min(abs(freqvect-defaults.highf));
%         LIMO.data.highf = freqvect(position);
%         LIMO.data.trim_highf = position; 
%     end
%     
%     LIMO.data.tf_freqs = freqvect(LIMO.data.trim_lowf:LIMO.data.trim_highf);
end

% deal with categorical and continuous regressors
if isnumeric(cat)
    LIMO.data.Cat = cat;
else
    if strcmp(cat(end-3:end),'.txt')
        LIMO.data.Cat = load(cat);
    elseif strcmp(cat(end-3:end),'.mat')
        name = load(cat); f = fieldnames(name);
        LIMO.data.Cat = getfield(name,f{1});
    else
        error('ERROR cat')
    end
end

if isnumeric(cont)
    LIMO.data.Cont = cont;
else
    if strcmp(cont(end-3:end),'.txt')
        LIMO.data.Cont = load(cont);
    elseif strcmp(cont(end-3:end),'.mat')
        [~,name,~] = fileparts(cont);
        load(cont); LIMO.data.Cont = eval(name);
    else
        error('ERROR cont')
    end
end

if ~exist('LIMO.dir','dir')
    mkdir(LIMO.dir)
end
cd(LIMO.dir); 
save LIMO LIMO; 
cd ..



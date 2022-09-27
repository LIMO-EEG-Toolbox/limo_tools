function limo_batch_import_data(EEG_DATA,cat,cont,defaults)

% routine to import 
%
% FORMAT limo_batch_import_data(EEG_DATA,cat,cont,defaults)
% 
% INPUT EEG_DATA is a an EEG file (temporary .set or .mat) to be loaded
%                from EEGLAB or FieldTrip
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
    [root,name,ext]          = fileparts(EEG_DATA);
    EEGLIMO.filepath         = root;
    EEGLIMO.filename         = [name ext];
end
LIMO.dir                     = defaults.name;
LIMO.data.data               = EEGLIMO.filename;
LIMO.data.data_dir           = EEGLIMO.filepath;
LIMO.Type                    = defaults.type;
LIMO.design.zscore           = defaults.zscore;
LIMO.design.method           = defaults.method;
LIMO.design.type_of_analysis = defaults.type_of_analysis;
LIMO.design.fullfactorial    = defaults.fullfactorial;
LIMO.design.bootstrap        = defaults.bootstrap;
LIMO.design.tfce             = defaults.tfce;
LIMO.design.status           = 'to do';
LIMO.Level                   = 1;

if strcmp(LIMO.data.data(end-3:end),'.set') % EEGLAB
    LIMO.Analysis                = defaults.analysis;
    LIMO.data.sampling_rate      = EEGLIMO.srate;
elseif ~strcmp(ft_datatype(EEGLIMO),'unknown')
    % this seems a data structure according to FieldTrip specs
    dtype = ft_datatype(EEGLIMO);
    switch dtype
        case 'raw'
            error('ERROR in limo_batch_import_data: FieldTrip ''raw'' data structures are not supported, convert to a ''timelock'' representation first');
        case 'timelock'
            % check whether the data has a 'trial' field
            if ~isfield(EEGLIMO, 'trial')
                error('ERROR in limo_batch_import_data: FieldTrip ''timelock'' data structures need a ''trial'' field');
            end
            if ~isfield(EEGLIMO, 'chanlocs')
                EEGLIMO = limo_get_ft_chanlocs(EEGLIMO, defaults);
            end
            LIMO.Analysis = 'Time';
            EEGLIMO.timevect = EEGLIMO.time*1000;
            LIMO.data.sampling_rate = mean(diff(EEGLIMO.time));
        case 'freq'
          % check whether the data has a 'powspctrm' field
          if ~isfield(EEGLIMO, 'powspctrm')
              error('ERROR in limo_batch_import_data: FieldTrip ''freq'' data structures require a ''powspctrm'' field');
          end
          if ~isfield(EEGLIMO, 'chanlocs')
              EEGLIMO = limo_get_ft_chanlocs(EEGLIMO, defaults);
          end
          if isfield(EEGLIMO, 'time')
            % time-freq
            LIMO.Analysis = 'Time-Frequency';
            EEGLIMO.tf_times = EEGLIMO.time*1000;
            EEGLIMO.tf_freqs = EEGLIMO.freq;
          else
            LIMO.Analysis = 'Frequency';
            EEGLIMO.freqvect = EEGLIMO.freq;
          end
          
        case 'source'
            error('ERROR in limo_batch_import_data: FieldTrip ''source'' data structures are not (yet) supported');
      otherwise
          error('ERROR in limo_batch_import_data: FieldTrip ''%s'' data structures are not supported', dtype); 
    end
else
    error('ERROR in limo_batch_import_data: neither EEGLAB nor FieldTrip data')
end


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
    
    if ~isfield(EEGLIMO.etc,'timeerp')
        try
            data     = load('-mat',EEGLIMO.etc.timeerp);
        catch erperr
            warning(erperr.identifier,'no timerp %s, trying something else',erperr.message)
            fprintf('the fied EEG.etc.timeerp is missing - looking for %s\n',EEGLIMO.data);
            if exist(fullfile(EEGLIMO.filepath,EEGLIMO.data),'file')
                data = load('-mat',fullfile(EEGLIMO.filepath,EEGLIMO.data));
            else
                data = load('-mat',fullfile(pwd,EEGLIMO.data));
            end
        end
        timevect = data.times; clear data;
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
    
elseif strcmp(LIMO.Analysis,'Frequency') 
    
    if ~isfield(EEGLIMO.etc,'freqspec')
        try
            data     = load('-mat',EEGLIMO.etc.freqspec);
        catch freqerr
            warning(freqerr.identifier,'no freqspec %s, trying something else',freqerr.message)
            fprintf('the fied EEG.etc.freqspec is missing - looking for %s\n',EEGLIMO.data);
            if exist(fullfile(EEGLIMO.filepath,EEGLIMO.data),'file')
                data = load('-mat',fullfile(EEGLIMO.filepath,EEGLIMO.data));
            else
                data = load('-mat',fullfile(pwd,EEGLIMO.data));
            end
        end
    end
    
    if isfield(EEGLIMO, 'freqvect')
        freqvect = EEGLIMO.freqvect;
    elseif ~isfield(EEGLIMO.etc,'freqspec')
        disp('the fied EEG.etc.freqspec is missing - reloading single trials');
        data     = load('-mat',EEGLIMO.etc.freqspec);
    else
        freqvect = EEGLIMO.etc.freqspec;
    end

    % start
    if isempty(defaults.lowf) || defaults.lowf < freqvect(1)
        LIMO.data.start = freqvect(1);
        LIMO.data.trim1 = 1;
    else
        [~,position]    = min(abs(freqvect-defaults.lowf));
        LIMO.data.start = freqvect(position);
        LIMO.data.trim1 = position; 
    end
    
    % end
    if isempty(defaults.highf) || defaults.highf > freqvect(end)
        LIMO.data.end   = freqvect(end);
        LIMO.data.trim2 = numel(freqvect);
    else
        [~,position]    = min(abs(freqvect-defaults.highf));
        LIMO.data.end   = freqvect(position);
        LIMO.data.trim2 = position; 
    end
    
    LIMO.data.freqlist  = freqvect(LIMO.data.trim1:LIMO.data.trim2);

elseif strcmp(LIMO.Analysis,'Time-Frequency')
    
    if ~isfield(EEGLIMO.etc,'timeersp') || ~isfield(EEGLIMO.etc,'freqersp')
        try
            data = load('-mat',EEGLIMO.etc.timef,'times','freqs');
        catch timeferr
            warning(timeferr.identifier,'error loading data %s\n, trying simeting else',timeferr.message)
            fprintf('ersp fied in EEG.etc absent or impcomplete, - looking for %s\n',EEGLIMO.data);
            if exist(fullfile(EEGLIMO.filepath,EEGLIMO.data),'file')
                data = load('-mat',fullfile(EEGLIMO.filepath,EEGLIMO.data));
            else
                data = load('-mat',fullfile(pwd,EEGLIMO.data));
            end
        end

    if isfield(EEGLIMO, 'tf_times') && isfield(EEGLIMO, 'tf_freqs')
        timevect = EEGLIMO.tf_times;
        freqvect = EEGLIMO.tf_freqs;
    elseif ~isfield(EEGLIMO.etc,'timeersp') || ~isfield(EEGLIMO.etc,'freqersp')
        disp('ersp fied in EEG.etc absent or impcomplete, reloading the single trials')
        data = load('-mat',EEGLIMO.etc.timef,'times','freqs');
    else
        timevect = EEGLIMO.etc.timeersp;
        freqvect = EEGLIMO.etc.freqersp;
    end
       
    % start
    if isempty(defaults.start) || defaults.start < min(timevect)
        LIMO.data.start = timevect(1);
        LIMO.data.trim1 = 1;    
    else
        [~,position]    = min(abs(timevect - defaults.start));
        LIMO.data.start = timevect(position);
        LIMO.data.trim1 =  find(timevect == LIMO.data.start);
    end
    
    % end
    if isempty(defaults.end) || defaults.end > max(timevect)
        LIMO.data.end   = timevect(end);
        LIMO.data.trim2 = length(timevect);    
    else
        [~,position]    = min(abs(timevect - defaults.end));
        LIMO.data.end   = timevect(position);
        LIMO.data.trim2 =  position;
    end

    LIMO.data.tf_times  = timevect(LIMO.data.trim1:LIMO.data.trim2);

    % start
    if isempty(defaults.lowf) || defaults.lowf < freqvect(1)
        LIMO.data.lowf = freqvect(1);
        LIMO.data.trim_lowf = 1;
    else
        [~,position] = min(abs(freqvect-defaults.lowf));
        LIMO.data.lowf = freqvect(position);
        LIMO.data.trim_lowf = position; 
    end
    
    % end
    if isempty(defaults.highf) || defaults.highf > freqvect(end)
        LIMO.data.highf = freqvect(end);
        LIMO.data.trim_highf = length(freqvect);
    else
        [~,position] = min(abs(freqvect-defaults.highf));
        LIMO.data.highf = freqvect(position);
        LIMO.data.trim_highf = position; 
    end
    
    LIMO.data.tf_freqs = freqvect(LIMO.data.trim_lowf:LIMO.data.trim_highf);
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
save(fullfile(LIMO.dir, 'LIMO.mat'), 'LIMO'); 

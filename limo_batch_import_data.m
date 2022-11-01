function limo_batch_import_data(setfile,cat,cont,defaults)

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
%  Copyright (C) LIMO Team 2019

global EEGLIMO

EEGLIMO                      = load('-mat',setfile);
EEGLIMO                      = EEGLIMO.EEG;
[root,name,ext]              = fileparts(setfile); 
LIMO.dir                     = defaults.name;
LIMO.data.data               = [name ext];
LIMO.data.data_dir           = root;
LIMO.data.sampling_rate      = EEGLIMO.srate;
LIMO.Analysis                = defaults.analysis;
LIMO.Type                    = defaults.type;
LIMO.design.zscore           = defaults.zscore;
LIMO.design.method           = defaults.method;
LIMO.design.type_of_analysis = defaults.type_of_analysis;
LIMO.design.fullfactorial    = defaults.fullfactorial;
LIMO.design.bootstrap        = defaults.bootstrap;
LIMO.design.tfce             = defaults.tfce;
LIMO.Level                   = 1;
LIMO.design.status           = 'to do';
if isfield(defaults, 'labels')
    LIMO.design.labels          = defaults.labels;
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
        data     = alternative_load('EEG.etc.timeerp',setfile,LIMO);
        timevect = data.times; clear data;
    else
        timevect = EEGLIMO.etc.timeerp;
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
    if isempty(defaults.end) || defaults.end > max(EEGLIMO.times)
        LIMO.data.end   = timevect(end);
        LIMO.data.trim2 = length(timevect);
    else
        [~,position]    = min(abs(EEGLIMO.times - defaults.end));
        LIMO.data.end   = timevect(position);
        LIMO.data.trim2 = position;
    end
    
    LIMO.data.timevect  = timevect(LIMO.data.trim1:LIMO.data.trim2);
    
elseif strcmp(defaults.analysis,'Frequency') 
    
    if ~isfield(EEGLIMO.etc,'freqspec')
        data     = alternative_load('EEG.etc.freqspec',setfile,LIMO);
        freqvect = data.freqs; clear data;
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

elseif strcmp(defaults.analysis,'Time-Frequency')
    
    if ~isfield(EEGLIMO.etc,'timeersp') || ~isfield(EEGLIMO.etc,'freqersp')
        data     = alternative_load('EEG.etc.timeersp .freqersp',setfile,LIMO);
        timevect = data.times;
        freqvect = data.freqs;
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
    else
        strcmp(cat(end-3:end),'.mat')
        name = load(cat); f = fieldnames(name);
        LIMO.data.Cat = getfield(name,f{1});
    end
end

if isnumeric(cont)
    LIMO.data.Cont = cont;
else
    if strcmp(cont(end-3:end),'.txt')
        LIMO.data.Cont = load(cont);
    else
        strcmp(cont(end-3:end),'.mat')
        [~,name,~] = fileparts(cont);
        load(cont); LIMO.data.Cont = eval(name);
    end
end

if ~exist('LIMO.dir','dir')
    mkdir(LIMO.dir)
end
cd(LIMO.dir); 
save LIMO LIMO; 
cd ..

% if using Import rather than STUDY we need to find the data another way
function  data = alternative_load(EEGfieldname,setfile,LIMO)

try
    data = pop_loadset(setfile);
catch fielderr
    warning(fielderr.identifier,'no data found %s, trying something else',fielderr.message)
    fprintf('the field %s is missing - looking for %s\n',EEGfieldname,LIMO.data.data);
    if exist(fullfile(LIMO.data.data_dir,LIMO.data.data),'file')
        try
            data = pop_loadset(fullfile(LIMO.data.data_dir,LIMO.data.data));
        catch
            data = load('-mat',fullfile(LIMO.data.data_dir,LIMO.data.data));
            data = data.EEG;
        end
    else
        try
            data = pop_loadset(fullfile(pwd,LIMO.data.data));
        catch
            data = load('-mat',fullfile(pwd,LIMO.data.data));
            data = data.EEG;
        end
    end
end


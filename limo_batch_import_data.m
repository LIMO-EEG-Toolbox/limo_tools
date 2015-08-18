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
%% -----------------------------
% Copyright (C) LIMO Team 2015

global EEGLIMO

EEGLIMO                      = pop_loadset(setfile);
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
    
    if ~isfield(EEGLIMO.etc,'timeerp')
        disp('the fied EEG.etc.timeerp is missing using the default EEG.times');
        timevect = EEGLIMO.times;
    else
        timevect = EEGLIMO.etc.timeerp;
    end
    
    % start
    if isempty(defaults.start) || defaults.start < min(timevect)
        LIMO.data.start = timevect(1);
        LIMO.data.trim1 = 1;    
    else
        [~,position]=min(abs(timevect - defaults.start));
        LIMO.data.start = timevect(position);
        LIMO.data.trim1 = position;
    end
    
    % end
    if isempty(defaults.end) || defaults.end > max(EEGLIMO.times)
        LIMO.data.end   = timevect(end);
        LIMO.data.trim2 = length(timevect);    
    else
        [~,position]=min(abs(EEGLIMO.times - defaults.end));
        LIMO.data.end = timevect(position);
        LIMO.data.trim2 = position;
    end
    
elseif strcmp(defaults.analysis,'Frequency') 
    
    if ~isfield(EEGLIMO.etc,'freqspec')
        error('can''t import spectrum the fied EEG.etc.freqspec (all freq computed) is missing')
    else
        freqvect = EEGLIMO.etc.freqspec;
    end

    % start
    if isempty(defaults.lowf) || defaults.lowf < freqvect(1)
        LIMO.data.start = freqvect(1);
        LIMO.data.trim1 = 1;
    else
        [~,position] = min(abs(freqvect-defaults.lowf));
        LIMO.data.start = freqvect(position);
        LIMO.data.trim1 = position; 
    end
    
    % end
    if isempty(defaults.highf) || defaults.highf > freqvect(end)
        LIMO.data.end = freqvect(end);
        LIMO.data.trim2 = numel(freqvect);
    else
        [~,position] = min(abs(freqvect-defaults.highf));
        LIMO.data.end = freqvect(position);
        LIMO.data.trim2 = position; 
    end
    
    LIMO.data.freqlist = freqvect(LIMO.data.trim1:LIMO.data.trim2);

elseif strcmp(defaults.analysis,'Time-Frequency')
    
    if ~isfield(EEGLIMO.etc,'freqersp')
        error('can''t import ersp the fied EEG.etc.freqsersp (all freq computed) is missing')
    elseif ~isfield(EEGLIMO.etc,'timeersp')
        error('can''t import ersp the fied EEG.etc.timeersp (all times computed) is missing')
    else
        timevect = EEGLIMO.etc.timeersp;
        freqvect = EEGLIMO.etc.freqersp;
    end
       
    % start
    if isempty(defaults.start) || defaults.start < min(timevect)
        LIMO.data.start = timevect(1);
        LIMO.data.trim1 = 1;    
    else
        [~,position]=min(abs(timevect - defaults.start));
        LIMO.data.start = timevect(position);
        LIMO.data.trim1 =  find(timevect == LIMO.data.start);
    end
    
    % end
    if isempty(defaults.end) || defaults.end > max(timevect)
        LIMO.data.end = timevect(end);
        LIMO.data.trim2 = length(timevect);    
    else
        [~,position]=min(abs(timevect - defaults.end));
        LIMO.data.end = timevect(position);
        LIMO.data.trim2 =  position;
    end

    LIMO.data.tf_times = timevect(LIMO.data.trim1:LIMO.data.trim2);

    % start
    if isempty(defaults.lowf) || defaults.lowf < freqvect(1)
        LIMO.data.lowf = freqvect(1);
        LIMO.data.trim_low_f = 1;
    else
        [~,position] = min(abs(freqvect-defaults.lowf));
        LIMO.data.lowf = freqvect(position);
        LIMO.data.trim_low_f = position; 
    end
    
    % end
    if isempty(defaults.highf) || defaults.highf > freqvect(end)
        LIMO.data.highf = freqvect(end);
        LIMO.data.trim_high_f = length(freqvect);
    else
        [~,position] = min(abs(freqvect-defaults.highf));
        LIMO.data.highf = freqvect(position);
        LIMO.data.trim_high_f = position; 
    end
    
    LIMO.data.tf_freqs = freqvect(LIMO.data.trim_low_f:LIMO.data.trim_high_f);
end

% deal with categorical and continuous regressors
if isnumeric(cat)
    LIMO.data.Cat = cat;
else
    if strcmp(cat(end-3:end),'.txt')
        LIMO.data.Cat = load(cat);
    else strcmp(cat(end-3:end),'.mat')
        load(cat); LIMO.data.Cat = eval(cat(1:end-4));
    end
end

if isnumeric(cont)
    LIMO.data.Cont = cont;
else
    if strcmp(cont(end-3:end),'.txt')
        LIMO.data.Cont = load(cont);
    else strcmp(cont(end-3:end),'.mat')
        load(cont); LIMO.data.Cont = eval(cont(1:end-4));
    end
end


cd(LIMO.dir); save LIMO LIMO; cd ..

end

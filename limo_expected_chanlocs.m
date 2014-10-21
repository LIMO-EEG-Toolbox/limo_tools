function channeighbstructmat = limo_expected_chanlocs(varargin)

% This function loads an EEG dataset to create a file with the
% location of all expected channels and to create a neighbourhood
% distance matrix used to control for multiple comparisons.
%
% FORMAT
% limo_expected_chanlocs
% limo_expected_chanlocs(full path data set name)
% limo_expected_chanlocs(data set name, path)
% limo_expected_chanlocs(data set name, path,neighbour distance)
%
% INPUT
% data set name is the name of a eeglab.set
% path is the location of that file
% neighbour distance is the distance between channels to buid the
% neighbourhood matrix
%
% channeighbstructmat is the neighbourhood matrix
%
% See also LIMO_NEIGHBOURDIST LIMO_GET_CHANNEIGHBSTRUCMAT
%
% Guillaume Rousselet v1 11 June 2010
% Rewritten by Cyril Pernet so we don't have to know which subject has the
% largest channel description v2 16 July 2010 - update 18 July 2012 to get
% output channeighbstructmat so we can update subjects for tfce
% Marianne Latinus May 2014 - update to create a cap with a minimum number
% of subjects per electrodes ; loop through all subjects
% ----------------------------------
%  Copyright (C) LIMO Team 2014

%% variables set as defaults
min_subjects = 3; % we want at least 3 subjects per electrode

global EEG
current_dir = pwd;

%% ask if data are from one subject or a set then get data
% ---------------------------------------------------------
if nargin == 0
    quest = questdlg('Make the Expected Chanlocs file from 1 subject or search throughout a set of subjects?','Selection','Set','One','Cancel','Set');
    if strcmp(quest,'Cancel')
        return
    else
        FileName = []; PathName = []; FilterIndex = [];
    end
elseif nargin == 1
    quest = 'One';  FilterIndex = 1;
    [PathName,f,e] = fileparts(varargin{1}); 
    FileName = [f e]; 
elseif nargin >= 2
    FileName = varargin{1}; 
    PathName = varargin{2}; 
    if size(FileName,1) == 1
        quest = 'One';
    else
        quest = 'Skip';
        for n=1:size(FileName,1)
            [Paths{n},name,ext] = fileparts(FileName{n});
            Names{n} = [name ext];
            Files{n} = [Path{n} fielsep Names{n}];
        end
    end
    FilterIndex = 1;
else
    error('wrong number of arguments')
end

neighbourdist = [];
if nargin == 3
    neighbourdist = varargin{3};
end

if isempty(neighbourdist)
    neighbourdist = eval(cell2mat(inputdlg('enter neighbourhood distance','neighbourhood distance'))); % 0.37 for biosemi 128;
end

%% from 1 subject
% -----------------------
if strcmp(quest,'One')
    
    if isempty(FileName)
        [FileName,PathName,FilterIndex]=uigetfile('*.set','EEGLAB EEG dataset before electrode removal');
        if FilterIndex == 0
            return
        end
    end
    
    try
        if strcmp([PathName FileName],[EEG.filepath filesep EEG.filename])
            disp('Using Global variable EEG')
        else
            try
                EEG=pop_loadset([PathName filesep FileName]);
            catch
                EEG=pop_loadset([PathName FileName]);
            end
        end
        expected_chanlocs = EEG.chanlocs;
        fprintf('Data set %s loaded \n',FileName);
        [neighbours,channeighbstructmat] = limo_get_channeighbstructmat(EEG,neighbourdist);
        if sum(channeighbstructmat(:)) == 0
            msg = sprintf('the neighbouring matrix is empty, it''s likely a distance issue \n see imo_ft_neighbourselection.m');
            error(msg)
        end
        cd (current_dir);
        if nargout == 0
            save expected_chanlocs expected_chanlocs channeighbstructmat % save all in one file
            fprintf('expected_chanlocs & channeighbstructmatfile saved\n');
        end
    catch ME
        errordlg('pop_loadset eeglab issue','error');
    end
    
    
elseif strcmp(quest,'Set')   % from a set of subjects
    % -------------------------------------------
    
    %% get data
    [name,path,filt]=uigetfile({'LIMO.mat';'*.txt'; '*.mat'; '*.*'}, 'Pick a LIMO.mat (subject 1) or list', 'MultiSelect', 'on');
    if filt == 0
        return
    else
        [~,file,~]=fileparts(name);
        if strcmp(file,'LIMO.mat') % we go for multiple LIMO.mat by hand
            Names{index} = name;
            Paths{index} = path;
            Files{index} = sprintf('%s\%s',path,name);
            go = 1; 
            cd(current_dir); % go back to pwd...
        else
            go = 0;
            if strcmp(name(end-3:end),'.txt')
                name = importdata(name);
            elseif strcmp(name(end-3:end),'.mat')
                name = load([path name]);
                name = getfield(name,cell2mat(fieldnames(name)));
            end
            
            for f=1:size(name,1)
                if ~exist(name{f},'file')
                    errordlg(sprintf('%s \n file not found',FileName{f}));
                    return
                else
                    Files{f} = name{f};
                    [Paths{f}, n,e] = fileparts(name{f});
                    Names{f} = [n e];
                end
            end
        end
    end
    
    index = 2;
    while go == 1
    [name,path] = uigetfile('LIMO.mat',['select LIMO file subject ',num2str(index)]);
        if name == 0
            go = 0;
        else
            if ~strcmp(name,'LIMO.mat')
                error(['you selected the file ' name ' but a LIMO.mat file is expected']);
            else
                Names{index} = name;
                Paths{index} = path;
                Files{index} = sprintf('%s\%s',path,name);
                cd(current_dir)
                index = index + 1;
            end
        end
    end
    
    %% retreive all chanlocs and make up a cap where we have a least 3 subjects
    chanlocs = cell(length(Paths),1);
    size_chanlocs = zeros(length(Paths),1);
    
    % retreive all chanlocs
    for i=1:length(Paths)
        load(Files{i})
        chanlocs{i} = LIMO.data.chanlocs;
        size_chanlocs(i) = size(LIMO.data.chanlocs,2);
        clear LIMO
        for c = 1:size_chanlocs(i)
            chan_labs{i,c} = chanlocs{i}(c).labels;
        end
    end
    
    % take the largest set as reference
    [nm,ref] = max(size_chanlocs);
    load(Files{ref})
    EEG.xmin = LIMO.data.start;
    EEG.xmax = LIMO.data.end;
    EEG.pnts = length(LIMO.data.start:1000/LIMO.data.sampling_rate:LIMO.data.end); % note only for LIMO v2 in msec
    EEG.chanlocs = LIMO.data.chanlocs;
    EEG.srate = LIMO.data.sampling_rate;
    EEG.trials = size(LIMO.design.X,1);
    clear LIMO
    for c = 1:nm
        ref_chan_labs{c,1} = chan_labs{ref,c};
        counter(c) = 1;
    end
    
    % loop on subjects
    for i = 1:size(chan_labs,1)
        if i ~= ref % skip reference subject
            n = size_chanlocs(i);
            for c = 1:n
                tmp{c} = chan_labs{i,c};
            end
            
            new_chans = setdiff(tmp, ref_chan_labs);
            if isempty(new_chans)
                counter  = counter + ismember(ref_chan_labs, tmp);
            else
                ref_chan_labs = [ref_chan_labs;new_chans']; % add channel
                try
                    counter = [counter;zeros(length(new_chans),1)] + ismember(ref_chan_labs, tmp);
                catch dim_issue
                    counter = [counter';zeros(length(new_chans),1)] + ismember(ref_chan_labs, tmp);
                end
                load(Files{i}) % load LIMO to get chanlocs of chans to add
                for j = 1:length(LIMO.data.chanlocs)
                    if ismember(LIMO.data.chanlocs(j).labels, new_chans)
                        EEG.chanlocs = [EEG.chanlocs LIMO.data.chanlocs(j)];
                    end
                end
                
            end
        end
    end
    
    % extra-check to remove external channel
    index = 1; remove = 0;
    for i=1:size(EEG.chanlocs,2)
        if strncmp(EEG.chanlocs(i).labels,'EX',2) || strncmp(EEG.chanlocs(i).labels,'ex',2)
            fprintf('likely external channel detected %s\n',EEG.chanlocs(i).labels)
            answer = input('Do you want to remove it [Y/N]: ','s');
            if strncmp(answer,'Y',1) || strncmp(answer,'y',1)
                remove(index) = i;
                index = index +1;
            end
        end
    end
    if remove ~=0
        EEG.chanlocs(remove) = [];
    end
    
    % remove low count
    EEG.chanlocs(find(counter < min_subjects)) = [];
    expected_chanlocs = EEG.chanlocs;
    
    % make up fake data
    EEG.nbchan = length(EEG.chanlocs);
    EEG.data = zeros(EEG.nbchan, EEG.pnts, EEG.trials);
    cd (current_dir);

    % now we have 1 cap we can do as if we had a single subject to process
    [neighbours,channeighbstructmat] = limo_get_channeighbstructmat(EEG, neighbourdist);
    if sum(channeighbstructmat(:)) == 0
        msg = sprintf('the neighbouring matrix is empty, it''s likely a distance issue \n see imo_ft_neighbourselection.m');
        error(msg)
    end
    
    if nargout == 0
        save expected_chanlocs expected_chanlocs channeighbstructmat % save all in one file
        fprintf('expected_chanlocs & channeighbstructmatfile saved\n');
    end
end


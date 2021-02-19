function filepath = limo_random_select_multivariate(type,varargin)

% This function is used to combine timevectors (classification accuracies/Fvalues/..) 
% computed at the 1st level using limo_mglm. 
% Whereas in limo_mglm observations are assumed independent
% (i.e. N-way ANOVA/ANCOVA or Regression), limo_random_effect_multivariate distinguishes
% independents and non-independent (repeated) measures.
%
% Note that no statistical test is done in limo_random_select_multivariate, only the grouping
% and organization of the data (hence the name select) - once data are
% selected and re-organized they are send to limo_random_robust_multivariate which deals
% with the data structures and call the stat functions
%
% FORMAT
% limo_random_select_multivariate(type)
% limo_random_select_multivariate(type,'nboot',1000,'tfce',1)
% limo_random_select_multivariate(type,'nboot',nbootval,'tfce',tfceval);
%
% INPUT
% type = 1 for a one sample t-test
% type = 2 for a two-samples t-test
% type = 3 for a paired t-test
% type = 4 for a regression
% type = 5 for an ANOVA
%
% Optional inputs:
%  'type'           - 'ica' or 'chan'
%  'nboot'          - the number of bootstrap to do (default = 0)
%  'tfce'           - 0/1 indicates to computes tfce or not (default = 0)
%  'limofiles'      - Cell array with the full paths to the subject file or file
%                     who contains the list of path to a group of sets. The
%                     dimensions of the cell correspond with group, factor and
%                     level respectively. (no default)
%  'regfile'        - Full path to Regressor file in case type = 4 (no default)
%  'folderprefix'   - Prefix for the results folder.
%  'folderpath'     - Path to save the results. Default is current
%                     directory
% OUTPUT
% filepath - Path to the contrast result file. Mainly for EEGALB functionality to
%            allow loading test directly.
%
% Cyril Pernet - The University of Edinburgh
% Ramon Martinez-Cancino - UC San Diego
% Iege Bassez 
% ---------------------------------------------------------
% Copyright (C) LIMO Team 2015

%% take the inputs and load some parameters

g = finputcheck(varargin, { 'nboot'          'integer'  0                             0          ;     % Bootstrap
    'tfce'           'integer'  [0 1]                          0          ;     % tfce
    'limofiles'      'cell'     {}                             {}         ;     % Path to subject file or group file Cell array with dimensions {group,,level}
    'regfile'        'string'   ''                             ''         ;     % Path to regressor files
    'folderprefix'   'string'   ''                             ''         ;     % Prefix for folder to save
    'folderpath'     'string'   ''                             ''         ;     % Path to folder to save
    'type'           'string'   {'Channels','Components'}      'Channels' ;});     % Type of measure ['ica', 'chan']
if isstr(g), error(g); end;


% check chanlocs and g.nboot
global limo
if ~isempty(g.folderpath)
    cd(g.folderpath);
end
limo.dir = pwd;

limo.design.bootstrap = g.nboot;
limo.design.tfce = g.tfce;
limo.Level = 2;
limo.Type = g.type;
limo.Analysis = 'time';

% ----------------------------------
%%  One sample t-test and regression
% ----------------------------------
if type == 1 || type == 4
    
    % get files
    % ---------
    if isempty(g.limofiles)
        [Names,Paths,limo.data.data] = limo_get_files;
        % Case for path to the files
    elseif size(g.limofiles{1},1) == 1
        [Names,Paths,limo.data.data] = limo_get_files([],[],[],g.limofiles{1});
        % Case when all paths are provided
    elseif size(g.limofiles{1},1) > 1
        [Names,Paths,limo.data.data] = breaklimofiles(g.limofiles{1});
    end
    
    if isempty(Names)
        disp('no files selected')
        return
    end
    limo.data.data_dir = Paths;
    
    % get info
    % ------------
    [subj_chanlocs, timevec, k] = get_info(Paths);
    limo.timevector = timevec(1,:);
    limo.nb_conditions_fl = k(1);
    
    % get data 
    % -----------------------------
    disp('gathering data ...'); 
    for i=1:size(Paths,2) % for each subject
        lc = load(limo.data.data{i});
        %tmp = eval(str2mat(Names{1}(1:end-4)));
        tmp = lc.Linear_Classification;
        tmp = squeeze(tmp(:,2,:)); % pick cv accuracies
        data(:,i) = tmp;
    end
        clear tmp
end
        
% one-sample t-test
% -----------------
if type == 1
    limo.design.X = ones(size(data,2),1);
    LIMO = limo;

    % clear some memory
    %clear Names Paths  expected_chanlocs limo subj_chanlocs
    cd(limo.dir);

    LIMO.design.method = 'Trimmed means'; save LIMO LIMO
    Yr = data; save Yr Yr, clear Yr % just to be consistent with name
    filepath = limo_random_robust_multivariate(type,data,1,g.nboot,g.tfce);

    % regression
    % -------------
elseif type == 4

    % ------------------------------
    %%  Two samples t-test
    % -----------------------------
elseif type == 2
    
    % ------------------------------
    %%  Paired t-test
    % -----------------------------
elseif type == 3
    
    % -----------------------------------
    %%  Various sorts of ANOVAs/ANCOVAs
    % -----------------------------------
elseif type == 5

end % closes the function


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% getting info from subjects subfunction
function [subj_chanlocs, subj_timevec, subj_k] = get_info(Paths)

disp('match frames between subjects ...')
% check Paths format
if iscell(Paths{1})
    tmp = Paths; clear Paths
    index = 1;
    for gp=1:size(tmp,2)
        for s=1:size(tmp{gp},2)
            Paths{index} = tmp{gp}(s);
            index = index + 1;
        end
    end
end

% now loop loading the LIMO.mat for each subject to collect information
for i=1:size(Paths,2)
    try
        cd (Paths{i});
    catch
        cd (cell2mat(Paths{i}))
    end
    load LIMO;
    
    if i==1
        Analysis = LIMO.Analysis;
    else
        if ~strcmp(LIMO.Analysis,Analysis)
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
        end
    end
    
    sampling_rate(i)          = LIMO.data.sampling_rate;
    if strcmpi(LIMO.Type,'Channels')
        subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
        subj_timevec(i,:) = linspace(LIMO.data.start, LIMO.data.end, LIMO.data.trim2 - LIMO.data.trim1 + 1);
        subj_k(i) = LIMO.design.nb_conditions;
        try
            channeighbstructmat = LIMO.data.channeighbstructmat;
        catch ME
        end
        
    else
        subj_chanlocs(i).chanlocs = [];
    end
end

    % quick check things are ok
    if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
        error('data have different sampling rates')
    end
end

function [Names,Paths,Files] = breaklimofiles(cellfiles)
    for ifiles = 1:size(cellfiles,1)
        [Paths{ifiles} filename ext] = fileparts(cellfiles{ifiles});
        Names{ifiles} = [filename ext];
        Files{ifiles} = fullfile(Paths{ifiles},[filename ext]);
    end
end
end
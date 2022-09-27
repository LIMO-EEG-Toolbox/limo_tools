function limo_check_neighbourghs(LIMO)

% routine function to check is a chanloc and neighbouring matrix are
% present - also proposes cleanup
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019

% is the neighbouring matrix present
% -----------------------------------
if ~isfield(LIMO.data,'neighbouring_matrix')
    answer = questdlg('load or compute neighbouring matrix?','channel neighbouring definition','Load','Compute','Compute');
    if strcmp(answer,'Load')
        [file,newpath,whatsup] = uigetfile('*.mat','select neighbourghing matrix (or expected chanloc file)');
        if whatsup == 0
            disp('selection aborded');
            return
        else
            neighbouringmatrix  = load([newpath filesep file]);
            fields              = fieldnames(neighbouringmatrix);
            for f=1:length(fields)
                if isstruct(neighbouringmatrix.(fields{f}))
                    expected_chanlocs = neighbouringmatrix.(fields{f});
                elseif ismatrix(neighbouringmatrix.(fields{f}))
                    channeighbstructmat = neighbouringmatrix.(fields{f});
                end
            end
        end
    else
        if ~exist(fullfile(LIMO.data.data_dir,LIMO.data.data),'file')
            % typically tmp file from STUDY
            tmp = dir([LIMO.data.data_dir filesep '*.set']);
            channeighbstructmat = limo_expected_chanlocs(tmp(1).name, LIMO.data.data_dir);
        else
            channeighbstructmat = limo_expected_chanlocs(LIMO.data.data, LIMO.data.data_dir);
        end
    end
    % update and save
    LIMO.data.expected_chanlocs   = expected_chanlocs;
    if ~all(unique(channeighbstructmat) == [0 1]')
        error('the neighbouring matrix is not binary')
    else
        LIMO.data.neighbouring_matrix = channeighbstructmat;
    end
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
end

% do we need cleanup of ex channels
% ----------------------------------
index = 1; remove = 0;
for i=1:size(LIMO.data.chanlocs,2)
    if strncmp(LIMO.data.chanlocs(i).labels,'EX',2) || strncmp(LIMO.data.chanlocs(i).labels,'ex',2)
        fprintf('likely external channel detected %s\n',LIMO.data.chanlocs(i).labels)
        answer = input('Do you want to remove it [Y/N]: ','s');
        if strncmp(answer,'Y',1) || strncmp(answer,'y',1)
            remove(index) = i;
            index = index +1;
        end
    end
end

if remove ~=0
    LIMO.data.chanlocs(remove) = [];
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    
    % does neighbourgh and chan_loc match
    % -----------------------------------
    if ~any(length(LIMO.data.neighbouring_matrix) == size(LIMO.data.chanlocs))
        warndlg2('there is a mismatch between the number of channels and the neighbouring matrix - edit now')
    end
end

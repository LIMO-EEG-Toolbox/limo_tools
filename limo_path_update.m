function limo_path_update(filesin)

% Routine to update all LIMO path in LIMO.mat files
%
% FORMAT: limo_path_update
%         limo_path_update(filesin)
%
% INPUT: filesin can be empty in which case user is prompted
%        filesin can a .txt file list of all the LIMO.mat files to update
%        filesin can be a cell array of all the LIMO.mat files to update
%
% By default it updates the LIMO.mat to it's new path and look for the  associated data,
% either in the same directory or above for 1st level LIMO files. In cases where data 
% are not found, a question dialogue pops out asking the user if he/she want to select
% the relevant data set
%
% Cyril Pernet v3
% ------------------------------
%  Copyright (C) LIMO Team 2020

source_dir = pwd;
if nargin == 0
    filesin = [];
end

% 1st get the LIMO.mat files
% --------------------------
go    = 1;
index = 1;
if isempty(filesin)
    while go == 1
        [name,path] = uigetfile('LIMO.mat',['select LIMO file subject ',num2str(index),go]);
        if name == 0
            go = 0;
        else
            Names{index} = name;
            Paths{index} = path;
            Files{index} = sprintf('%s\%s',path,name);
        end
        cd(fileparts(Paths{index}));
        index = index + 1;
    end
    cd(source_dir)
else
    if ischar(filesin) 
        if strcmp(filesin(end-3:end),'.mat') % single file update
            [Paths{1},Names{1},ext] = fileparts(filesin);
            Files{1}                = fullfile(Paths{1},[Names{1} ext]);
        elseif strcmp(filesin(end-3:end),'.txt') % Case for path to the files .txt file
            [Names,Paths,Files]     = limo_get_files([],[],[],filesin);
        end
    else % cell-array
        N     = size(filesin,2);
        Paths = cell(1,N);
        Names = cell(1,N);
        Files = cell(1,N);
        for ifiles = 1:N
            [Paths{ifiles}, filename, ext] = fileparts(filesin{ifiles});
            Names{ifiles}                  = [filename ext];
            Files{ifiles}                  = fullfile(Paths{ifiles},[filename ext]);
        end
    end
end

% 2nd update the path and look for .set
% -------------------------------------
for i=size(Paths,2):-1:1
    LIMO = load(fullfile(Paths{i},'LIMO.mat'));
    LIMO = LIMO.LIMO;
    
    % same directory for LIMO.dir and data
    if strcmp(LIMO.dir,LIMO.data.data_dir) || strcmp(LIMO.dir,LIMO.data.data_dir(1:end-1)) % happens because from EEGLAB there is an extra /
        LIMO.dir = pwd;
        LIMO.data.data_dir = pwd;
        if LIMO.Level == 1
            name = dir(LIMO.data.data);
            if ~isempty(name)
                file_found(i) = 1;
                fprintf('LIMO file subject %g successfully updated\n',i);
            else
                file_found(i) = 0;
            end
            save LIMO LIMO
        end
        
    % different directory
    elseif LIMO.Level == 1 && length(LIMO.dir) > length(LIMO.data.data_dir)
        test = LIMO.dir(1:length(LIMO.data.data_dir)) == LIMO.data.data_dir;
        if sum(test) == length(LIMO.data.data_dir)  % data are located in a directory above LIMO
            current_dir = pwd; LIMO.dir = current_dir;
            [LIMO.data.data_dir,name,ext]=fileparts(current_dir); 
            cd(LIMO.data.data_dir)
            name = dir(LIMO.data.data);
            if ~isempty(name)
                file_found(i) = 1;
                fprintf('LIMO file subject %g successfully updated\n',i);
            else
                file_found(i) = 0;
            end
        end
        cd(LIMO.dir); save LIMO LIMO
                
    else
        file_found(i) = 0;
    end
end


% 3rd update missing .set
% ------------------------
if LIMO.Level == 1 && sum(file_found) ~= length(file_found)
    q = questdlg('.set missing - do you want to select data or simply update the LIMO file','association between file missing','Only update LIMO','Let me select .set','Only update LIMO');
    if strcmp(q,'Let me select .set')
        for i=1:length(file_found)
            if file_found(i) == 0
                cd(Paths{i}); load LIMO
                [name,path] = uigetfile('*.set','data set missing - please select data' );
                % check
                if ~strcmp(name,LIMO.data.data)
                    errordlg('a different set was previously assicoated to this LIMO.mat file, please check')
                else
                    cd(path); LIMO.data.data_dir = pwd;
                    cd(Paths{i}); save LIMO LIMO
                end
            end
        end
    end
end



function limo_path_update(filesin)

% Routine to update all LIMO path in LIMO.mat files
%
% FORMAT: limo_path_update
%         limo_path_update(filesin,newpath)
%
% INPUT: filesin can be empty in which case user is prompted
%        filesin can a .txt file list of all the LIMO.mat files to update
%        filesin can be a cell array of all the LIMO.mat files to update
%        newpath is the new path to use
%
% By default it updates the LIMO.mat to it's new path and look for the associated data,
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
        [filesin,fpath] = uigetfile({'*.mat;*.txt'},'select LIMO file or list.txt of LIMO files to update');
        if isempty(filesin) 
            warning('selection aborded'); return
        else
            if strcmp(filesin,'LIMO.mat')
                Paths{1} = fpath;
                newpath  = fpath;
            else
                [~,Paths] = limo_get_files([],[],[],fullfile(fpath,filesin));
                newpath = uigetdir('select the new directory where all data are located');
                if isempty(newpath)
                    warning('selection aborded'); return
                end
            end
        end
else
    if ischar(filesin) 
        if strcmp(filesin(end-3:end),'.mat') % single file update
            Paths{1}  = fileparts(filesin);
        elseif strcmp(filesin(end-3:end),'.txt') % Case for path to the files .txt file
            [~,Paths] = limo_get_files([],[],[],filesin); fpath = Paths;
        end
    else % cell-array
        N     = size(filesin,2);
        Paths = cell(1,N);
        for ifiles = 1:N
            [Paths{ifiles}, filename, ext] = fileparts(filesin{ifiles});
        end
    end
    
    if ~exist(newpath,'dir')
        error('the newpath doesn''t exist?')
    end
end

% 2nd update the path and look for .set
% -------------------------------------
for i=size(Paths,2):-1:1
    % let's see if we can find common ground
    common_str = newpath(max(strfind(newpath,filesep))+1:end);
    if ~contains(Paths{i},common_str)
        warning('subject %g not updated, could not find a string in common (here:''%s'') between the new and old path',i,common_str)
    else
        tmp = fullfile(newpath,Paths{i}(strfind(Paths{i},common_str)+length(common_str):end));
        if ~strcmpi(Paths{i},newpath)  % if LIMO.mat picked-up of a single subject in the new path
            if exist(fullfile(tmp,'LIMO.mat'),'file')
                LIMO     = load(fullfile(tmp,'LIMO.mat'));
                LIMO     = LIMO.LIMO;
                LIMO.dir = tmp;
            else
                warning('cannot locate %s',fullfile(tmp,'LIMO.mat'))
            end
        else
            LIMO     = load(fullfile(Paths{i},'LIMO.mat'));
            LIMO     = LIMO.LIMO;
            LIMO.dir = Paths{i};
        end
        
        % same directory for LIMO.dir and data
        if strcmp(Paths{i},LIMO.data.data_dir) || strcmp(Paths{i},LIMO.data.data_dir(1:end-1)) % happens because from EEGLAB there is an extra /
            LIMO.data.data_dir = Paths{i};
            if LIMO.Level == 1
                name = dir(LIMO.data.data);
                if ~isempty(name)
                    file_found(i) = 1;
                    fprintf('LIMO file subject %g successfully updated\n',i);
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                else
                    file_found(i) = 0;
                    pathtested{i} = LIMO.dir;
                end
            elseif LIMO.Level == 2
                if ~isempty(dir(fullfile(LIMO.dir,'Yr*.mat')))
                    file_found(i) = 1;
                    fprintf('2nd level LIMO file %g successfully updated\n',i);
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end
            end
            
        else % different directory
            tmp = fullfile(newpath,LIMO.data.data_dir(strfind(LIMO.data.data_dir,common_str)+length(common_str):end));
            if ~isempty(dir([tmp filesep '*.set']))  % data are located in the same directory as LIMO
                LIMO.data.data_dir = tmp; file_found(i) = 1;
            elseif ~isempty(dir([fileparts(tmp) filesep '*.set']))  % could be on directory up
                LIMO.data.data_dir = fileparts(tmp); file_found(i) = 1;
            elseif ~isempty(dir(fullfile(LIMO.dir,'Yr*.mat')))
                LIMO.data.data_dir = fileparts(tmp); file_found(i) = 1;
            else
                file_found(i) = 0;
                pathtested{i} = LIMO.dir;
            end
            
            if file_found(i) == 1
                fprintf('LIMO file subject %g successfully updated\n',i);
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
            end
        end
    end
end


% 3rd update missing .set
% ------------------------
if sum(file_found) ~= length(file_found)
    q = questdlg('.set missing - do you want to select data or simply update the LIMO file','association between file missing','Only update LIMO','Let me select .set','Only update LIMO');
    if strcmp(q,'Let me select .set')
        for i=1:length(file_found)
            if file_found(i) == 0
                cd(pathtested{i}); LIMO = load(fullfile(pathtested{i},'LIMO.mat')); LIMO = LIMO.LIMO;
                [name,path] = uigetfile('*.set','data set missing - please select data' );
                % check
                if ~strcmp(name,LIMO.data.data)
                    errordlg('a different set was previously associated to this LIMO.mat file, please check')
                else
                    LIMO.data.data_dir = path;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end
            end
        end
    end
end

% 4th update txt files
if strcmp(filesin(end-3:end),'.txt')
    alltxt = dir([fpath '*.txt']);
    for f=1:size(alltxt,1)
        [Names,Paths,Files] = limo_get_files([],[],[],fullfile(alltxt(f).folder,alltxt(f).name));
        for p=1:length(Paths)
            Files{p} = fullfile(newpath,[Paths{i}(strfind(Paths{i},common_str)+length(common_str):end) filesep Names{p}]);
        end
        cell2csv(fullfile(alltxt(f).folder,alltxt(f).name),Files')
    end
end


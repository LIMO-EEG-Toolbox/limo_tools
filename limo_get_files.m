function [Names,Paths,Files] = limo_get_files(gp,filter)

% routine to get multifiles from different directories
%
% FORNMAT [Names,Paths,Files] = limo_get_files(gp,filter)
%
% INPUT can be left empty, in the case ask for .mat or .txt
%       gp  is a simple numerical value, so the selection question reminds
%           the user which group is getting selected (useful for eg anova
%           analyses)
%       filter allows to specify the type of files to load,
%              for instance {'*.txt;*.set'} means load mat or txt
%              supported formats are .mat .txt .set .study
%              default is {'*.mat;*.txt'}
%              e.g. [Names,Paths,Files] = limo_get_files([],{'*.set'})
%
% OUTPUT Names , Paths, Full File names are returned as cells
%
% Cyril Pernet 25-08-09
% Ncolas Chauveau 07-04-11 - allows txt files
% CP update for filter and additional format 01-21-2015
% -------------------------------------------------
%  Copyright (C) LIMO Team 2015

if nargin == 0
    gp = [];
    filter = {'*.mat;*.txt'};
end

go = 1; index = 1;
while go == 1
    if ~isempty(gp)
        [name,path] = uigetfile(filter,['select a subject file',num2str(index),' ',gp,' or list file']);
    else
        [name,path] = uigetfile(filter,['select a subject file or list file']);
    end
    
    if name == 0
        go = 0; % exit
    elseif strcmp(name(end-2:end),'mat') || strcmp(name(end-2:end),'set') % select mat files
        Names{index} = name;
        Paths{index} = path;
        Files{index} = [path name];
        cd(path); cd ..
        index = index + 1;
    elseif strcmp(name(end-4:end),'study')  % select study file
        load('-mat', name);       
        for f=1:size(STUDY.datasetinfo,2)
            Names{f} = STUDY.datasetinfo(f).filename;
            Paths{f} = STUDY.datasetinfo(f).filepath;
            Files{f}=[Paths{f} filesep Names{f}];
        end
        index = f; go = 0;
    elseif strcmp(name(end-2:end),'txt')
        group_files=textread([path name],'%s'); % ,'whitespace','');  % select a txt file listing all files
        for f=1:size(group_files,1)
            Files{f} = group_files{f};
            [Paths{f},NAME,EXT] = fileparts(group_files{f});
            Names{f}=[NAME EXT];
        end
        index = f; go = 0;
    else
        errordlg('format not supported'); go = 0;
    end
end

if index == 1
    Names = []; Paths = []; Files = [];
end
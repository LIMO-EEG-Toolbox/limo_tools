function [Names,Paths,Files] = limo_get_files(varargin)

% routine to get multifiles from different directories
%
% FORNMAT [Names,Paths,Files] = limo_get_files(gp,filter,title)
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
%      title a default question dialogue is 'select a subject file or list
%            file' with possible the gp number inserted, but this can be
%            customized here
%
% OUTPUT Names , Paths, Full File names are returned as cells
%
% Cyril Pernet 25-08-09
% Nicolas Chauveau 07-04-11 - allows txt files
% CP update for filter and additional format 01-21-2015
% Arnaud Delorme fixed delimiters and added study option June 2015
% --------------------------------------------------------------
%  Copyright (C) LIMO Team 2015

%% defaults and inputs
gp        = [];
title     = ['select a subject file or list file'];
filter    = {'*.mat;*.txt'};
path2file = [];

if nargin >= 1; gp        = varargin{1}; end
if nargin >= 2; filter    = varargin{2}; end
if nargin >= 3; title     = varargin{3}; end
if nargin == 4; path2file = varargin{4}; end

go = 1; index = 1;
while go == 1
    if isempty(path2file)
        if ~isempty(gp)
            [name,path] = uigetfile(filter,['select a subject file',num2str(index),' ',gp,' or list file']);
        else
            [name,path] = uigetfile(filter,title);
        end
    else
        if exist(path2file,'file') == 2
            [path,filename,filext] = fileparts(path2file);
            name = [filename filext]; clear filename filext;
        else
            error('A valid path to the file must be provided ');
        end
    end
    
    if name == 0
        go = 0; % exit
    elseif strcmp(name(end-2:end),'mat') || strcmp(name(end-2:end),'set') % select mat files
        Names{index} = name;
        Paths{index} = path;
        Files{index} = fullfile(path,name);
        cd(path); cd ..
        index = index + 1;
    elseif strcmp(name(end-4:end),'study')  % select study file
        load('-mat', name);       
        for f=1:size(STUDY.datasetinfo,2)
            Names{f} = STUDY.datasetinfo(f).filename;
            Paths{f} = STUDY.datasetinfo(f).filepath;
            Files{f} = fullfile(Paths{f},Names{f});
        end
        index = f; go = 0;
    elseif strcmp(name(end-2:end),'txt')
        group_files = textread(fullfile(path,name),'%s','delimiter','');  % select a txt file listing all files
        for f=1:size(group_files,1)
            Files{f}            = group_files{f};
            [Paths{f},NAME,EXT] = fileparts(group_files{f});
            Names{f}            = [NAME EXT];
        end
        index = f; go = 0;
    else
        errordlg('format not supported'); go = 0;
    end
end

if index == 1
    Names = []; Paths = []; Files = [];
end
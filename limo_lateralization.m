function [LI_Map,channels] = limo_lateralization(varargin)

% computes a lateralization index and then it run a one_sample ttest on the
% LI_map to know if it is significanlty different from 0
%
% FORMAT [LI_Map, electrodepairs] = limo_lateralization(data,LIMO)
%        LI_Map = limo_lateralization(data,chanlocs,start,end,sampling_rate,channels)
%
% INPUT
% data is a 3D matrix [channels, time frames, subjects]
% - LIMO is a LIMO.mat which include the position of channels, times, etc
% - chanlocs is the expected chanel localtions
% - start and end are the time information
% - sampling_rate is the sampling rate in Hz
% - channels is a channel*2 matrix of channels to contrast (computed
% automaticaaly if a LIMO.mat is given instead)
%
% example for a 30 channels cap (veog eog removed)
% left =  [1 3 4 8  9  12 13 16 17 20 21 28 26]';
% right = [2 7 6 11 10 15 14 19 18 24 23 29 27]';
% channels = [left right];
%
% OUTPUT
% LI_map is a 3D matrix [channels, time frames, 5]
% the last dimension reports [mean value, std, nb_subjects, T values, p values])
% The lateralization index is computed as (left-right/left+right) * 100
% (or the opposite on the right sided channels to make a symetric map)
%
% boot_LI_map is a 4D matrix [channels, time frames, [T values under H1, T values under H0, p values under H0], nboot]
% this is used when looking for significant lateralization effects
%
% see also limo_pair_channels.m limo_LI.m
% 
% Cyril Pernet v1 25-05-2012
% ------------------------------
%  Copyright (C) LIMO Team 2021

%% check inputs
data = varargin{1};
if ischar(data)
    data = load(data);
    data = data.(cell2mat(fieldnames(data)));
end
[e,t,s]= size(data);

if nargin == 2
    tmp_LIMO = varargin{2};
    mkdir('lateralization'); 
    cd('lateralization');
    tmp_LIMO.dir = pwd;
    channels = limo_pair_channels(varargin{2}.data.chanlocs);
elseif nargin == 6
    tmp_LIMO.data.chanlocs            = varargin{2}.chanlocs;
    tmp_LIMO.data.neighbouring_matrix = varargin{2}.neighbouring_matrix;
    tmp_LIMO.data.start               = varargin{3};
    tmp_LIMO.data.end                 = varargin{4};
    tmp_LIMO.data.sampling_rate       = varargin{5};
    tmp_LIMO.design.electrode         = [];
    channels                          = varargin{6};
else
    error('wrong number of arguments in')
end
  

%% compute LI maps per subject

Maps = zeros(e,t,s);
if size(data,1) == 1
    if exist('errordlg2','file')
        errordlg2('data seem to have only one channel, no lateralization can be conputed','dimension issue')
    else
        errordlg('data seem to have only one channel, no lateralization can be conputed','dimension issue')
    end
end

for subject = 1:s
    tmp   = data(:,:,subject);
    left  = tmp(channels(:,1),:);
    right = tmp(channels(:,2),:);
    Maps(channels(:,1),:,subject) = ((left-right)./(left+right)).*100;
    Maps(channels(:,2),:,subject) = ((right-left)./(left+right)).*100;
end

%% compute the statistic for paired channels

% rename and copy the t-test file
save(fullfile(tmp_LIMO.dir,'LIMO.mat'),'tmp_LIMO')
limo_random_robust(1,Maps,1,tmp_LIMO);
LI_Map = load('one_sample_ttest_parameter_1.mat');
LI_Map = LI_Map.(cell2mat(fieldnames(LI_Map)));
save(fullfile(fileparts(tmp_LIMO.dir),'LI_Map.mat'),'LI_Map','-v7.3')
LI_Map = fullfile(fileparts(tmp_LIMO.dir),'LI_Map.mat');

% rename and copy the bootstrap file
if tmp_LIMO.design.bootstrap ~= 0
    H0_LI_Map = load(['H0' filesep 'H0_one_sample_ttest_parameter_1.mat']);
    H0_LI_Map = H0_LI_Map.(cell2mat(fieldnames(H0_LI_Map)));
    save(fullfile(fileparts(tmp_LIMO.dir),'H0','H0_LI_Map'),'H0_LI_Map','-v7.3')
end

% rename and copy the tfce files
if tmp_LIMO.design.tfce ~= 0
    tfce_LI_Map = load(['tfce' filesep 'tfce_one_sample_ttest_parameter_1.mat']);
    tfce_LI_Map = tfce_LI_Map.(cell2mat(fieldnames(tfce_LI_Map)));
    save(fullfile(fileparts(tmp_LIMO.dir),'tfce','tfce_LI_Map.mat'),'tfce_LI_Map','-v7.3')
    tfce_H0_LI_Map = load(['H0' filesep 'H0_one_sample_ttest_parameter_1.mat']);
    tfce_H0_LI_Map = tfce_H0_LI_Map.(cell2mat(fieldnames(tfce_H0_LI_Map)));
    save(fullfile(fileparts(tmp_LIMO.dir),'H0','tfce_H0_LI_Map.mat'),'tfce_H0_LI_Map','-v7.3')
end
    
cd(fileparts(tmp_LIMO.dir));
try
    rmdir(tmp_LIMO.dir,'s')
catch Er
    warning('temporary lateralization folder not removed - matlab error :%s\n',Er.message) %#ok<MEXCEP>
end

%% make a LIMO structure to be used in limo_display_results
if nargin > 2
    LIMO       = tmp_LIMO;
    LIMO.Level = 'LI';
    LIMO.dir   = pwd;
    LIMO.data.LIelectrodes = channels;
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
end
clear tmp_LIMO

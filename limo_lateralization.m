function [LI_Map, electrodes] = limo_lateralization(varargin)

% computes a lateralization index and then it run a one_sample ttest on the
% LI_map to know if it is significanlty different from 0
%
% FORMAT [LI_Map, electrodepairs] = limo_lateralization(data,LIMO)
%        LI_Map = limo_lateralization(data,chanlocs,start,end,sampling_rate,electrodes)
%
% INPUT
% data is a 3D matrix [electrodes, time frames, subjects]
% - LIMO is a LIMO.mat which include the position of electrodes, times, etc
% - chanlocs is the expected chanel localtions
% - start and end are the time information
% - sampling_rate is the sampling rate in Hz
% - electrodes is a electrode*2 matrix of electrodes to contrast (computed
% automaticaaly if a LIMO.mat is given instead)
%
% example for a 30 channels cap (veog eog removed)
% left =  [1 3 4 8  9  12 13 16 17 20 21 28 26]';
% right = [2 7 6 11 10 15 14 19 18 24 23 29 27]';
% electrodes = [left right];
%
% OUTPUT
% LI_map is a 3D matrix [electrodes, time frames, 5]
% the last dimension reports [mean value, std, nb_subjects, T values, p values])
% The lateralization index is computed as (left-right/left+right) * 100
% (or the opposite on the right sided electrodes to make a symetric map)
%
% boot_LI_map is a 4D matrix [electrodes, time frames, [T values under H1, T values under H0, p values under H0], nboot]
% this is used when looking for significant lateralization effects
%
% Cyril Pernet v1 25-05-2012
% ------------------------------
%  Copyright (C) LIMO Team 2019

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
    % check LIMO.data.chanlocs ; get theta and radius to be certain of position 
    radius = arrayfun(@(x) x.radius, varargin{2}.data.chanlocs);
    theta  = arrayfun(@(x) x.theta, varargin{2}.data.chanlocs);
    % remove midline electrodes
    remove = [find(theta==0) find(theta ==180)]; 
    radius(remove) = NaN;
    theta(remove)  = NaN;
    % create the matrix electrode
    Rtheta     = round(theta(theta>0));
    Rradius    = radius(find(theta>0));
    Lradius    = radius(find(theta<0));
    N          = size(Rtheta,2);
    electrodes = NaN(N,2);
    for n=1:N
        v = Rtheta(n);
        index = find(Rtheta == v); % always n and some other value(s)
        if length(index) > 1 % few electrodes at same angle
            position             = find(index == n);
            index2               = index(position);
            R_value              = find(round(theta) == v);
            electrodes(index2,2) = R_value(position);
            R_value              = Rradius(index);
            R_value              = R_value(position);
            [~,position]         = min(rem(Lradius(index),R_value)); % get electrode at 'same' radius in left side
            L_value              = find(round(theta) == -v);
            electrodes(index2,1) = L_value(position);
        else
            electrodes(index,2)  = find(round(theta) == v);
            electrodes(index,1)  = find(round(theta) == -v);
        end
    end
elseif nargin == 6
    tmp_LIMO.data.chanlocs = varargin{2}.chanlocs;
    tmp_LIMO.data.neighbouring_matrix = varargin{2}.neighbouring_matrix;
    tmp_LIMO.data.start = varargin{3};
    tmp_LIMO.data.end  = varargin{4};
    tmp_LIMO.data.sampling_rate = varargin{5};
    tmp_LIMO.design.electrode = [];
    electrodes = varargin{6};
else
    error('wrong number of arguments in')
end
  

%% compute LI maps per subject

Maps = zeros(e,t,s);
for subject = 1:s
    tmp = data(:,:,subject);
    left = tmp(electrodes(:,1),:);
    right = tmp(electrodes(:,2),:);
    Maps(electrodes(:,1),:,subject) = ((left-right)./(left+right)).*100;
    Maps(electrodes(:,2),:,subject) = ((right-left)./(left+right)).*100;
end

%% compute the statistic

% rename and copy the t-test file
if tmp_LIMO.design.tfce ~= 0
    LIMO = tmp_LIMO;
    LIMO.dir = pwd;
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
end
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
    tfce_LI_Map    = load(['tfce' filesep 'tfce_one_sample_ttest_parameter_1.mat']);
    tfce_LI_Map    = tfce_LI_Map.(cell2mat(fieldnames(tfce_LI_Map)));
    save(fullfile(fileparts(tmp_LIMO.dir),'tfce','tfce_LI_Map.mat'),'tfce_LI_Map','-v7.3')
    tfce_H0_LI_Map = load(['H0' filesep 'tfce_H0_one_sample_ttest_parameter_1.mat']);
    tfce_H0_LI_Map = tfce_H0_LI_Map.(cell2mat(fieldnames(tfce_H0_LI_Map)));
    save(fullfile(fileparts(tmp_LIMO.dir),'H0','tfce_H0_LI_Map.mat'),'tfce_H0_LI_Map','-v7.3')
end

cd(fileparts(tmp_LIMO.dir));
rmdir(tmp_LIMO.dir,'s')

%% make a LIMO structure to be used in limo_display_results
if nargin > 2
    LIMO       = tmp_LIMO;
    LIMO.Level = 'LI';
    LIMO.dir   = pwd;
    LIMO.data.LIelectrodes = electrodes;
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
end
clear tmp_LIMO

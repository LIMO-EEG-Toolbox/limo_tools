function clusters = limo_get_summary(file,mask)

% simple routine to compute summary statistics per cluster
%
% FORMAT clusters = limo_get_summary(file,mask)
%
% INPUTS file is a result file like a t-test or ANOVA
%        mask is a N-ary matrix of clusters
%
% OUTPUTS clusters, a structure with summary statistics per cluster
%               - eigenmode is the 'direction' of the effect size (usually,
%                        clusters are right skewed and thus it represents
%                        better the 'average' effect size)
%               - median provided as a comparison point to eigen mode
%               - mean provided as a comparison point to eigen mode
%               - min and max for completeness
%
% Cyril Pernet 2022
% ------------------------------
%  Copyright (C) LIMO Team 2022

clusters = [];

%% check inputs
if nargin == 0
    % no input, ask user to select a file
    [file,filepath] = uigetfile('.mat','select a LIMO stat file');
    if isempty(file)
        return
    else
        file = fullfile(filepath,file);
    end
    
    % no input, check if user want to use current mask
    ismask = evalin( 'base', 'exist(''mask'',''var'') == 1' );
    if ismask
        mask = evalin('base','mask'); 
    else
        if exist('errordlg2','file')
            errordlg2('no mask file found in the workspace, image the statistical results to create one and call this function again')
        else
            errordlg('no mask file found in the workspace, image the statistical results to create one and call this function again')
        end
    end
end

[filepath,filename,ext]=fileparts(file);
if isempty(filepath)
    filepath = pwd;
end
filename = [filename ext];
if ~exist(fullfile(filepath,filename),'file')
    error('file %s not found', filename)
end

if ~exist(fullfile(filepath,'LIMO.mat'),'file')
    error('cannot find a LIMO.mat in the same filder as this file, this is required for this function to work')
else
    LIMO = load(fullfile(filepath,'LIMO.mat'));
    LIMO = LIMO.LIMO;
end

%% get data

stats = load(fullfile(filepath,filename));
stats = stats.(cell2mat(fieldnames(stats)));

if contains(filename,'one_sample','IgnoreCase',true) || contains(filename,'two_samples','IgnoreCase',true) || ...
        contains(filename,'paired_samples','IgnoreCase',true) || contains(filename,'con_','IgnoreCase',true) || ...
        contains(filename,'ess_','IgnoreCase',true)
    if numel(size(stats)) == 3
        stats = squeeze(stats(:,:,4));
    else
        stats = squeeze(stats(:,:,:,4));
    end
elseif contains(filename,'R2.mat') || contains(filename,'semi_partial_coef')
    if numel(size(stats)) == 3
        stats = squeeze(stats(:,:,2));
    else
        stats = squeeze(stats(:,:,:,2));
    end
else
    if numel(size(stats)) == 3
        stats = squeeze(stats(:,:,1));
    else
        stats = squeeze(stats(:,:,:,1));
    end
end

%% deal with clusters
% --------------------
% quickly make if N-ary if binary
if length(unique(mask)) == 2
    mask = limo_findcluster(mask,LIMO.data.neighbouring_matrix,2);
end

num = unique(mask);
num(num==0) = [];
for c = length(num):-1:1
    data                  = stats(mask == num(c));
    clusters(c).eigenmode = sqrt(eig(data'*data)/length(data));
    clusters(c).median    = median(data);
    clusters(c).mean      = mean(data);
    clusters(c).min       = min(data);
    clusters(c).max       = max(data);
    if nargout == 0
       assignin('base','clusters_summary_stats',clusters) 
    end
end


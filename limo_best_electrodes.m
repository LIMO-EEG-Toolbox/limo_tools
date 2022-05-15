function [channel_vector,urchan_vector,freqmap] = limo_best_electrodes(varargin)

% This function finds the channel with the maximum F value in each subject.
% The function works on files R2.mat, Condition_effect.mat, Continuous.mat. 
% The function returns a map of frequency showing how often an electrode is
% selected across subjects. The default list of LIMO.mat can be used for R2.
%
% FORMAT [channel_vector,urchan_vector] = limo_best_electrodes(NameOfListOfFiles,expected_chanlocs)
%
% INPUT if empty user is prompted overwise the name of a file listing LIMO/R2/etc.. files
%       must be provided followed by an a channel location file
%
% OUTPUTS channel_vector the indices of which electrodes had the strongest F values
%         urchan_vector the value read from the urchan field (if present)
%         a frequency map is also presented as graphical output
%
% ------------------------------
%  Copyright (C) LIMO Team 2021

% Cyril Pernet & Guillaume Rousselet

%% file selection

current_dir    = pwd;
channel_vector = [];
urchan_vector  = [];

if nargin == 0
    [name,pathname,FilterIndex]=uigetfile({'*.mat','MAT-files (*.mat)'; '*.txt','Text (*.txt)'}, ...
        'Pick a list of result files (e.g. R2)');
else
   [pathname,name,ext]=fileparts(varargin{1}); name = [name ext];
   if exist(fullfile(pathname,name),'file')
       FilterIndex = 1;
   else
       warning('file %s not found',name);
       return
   end
end


if FilterIndex ~= 0
    
    if strcmp(name(end-3:end),'.txt')
        name = importdata(fullfile(pathname,name));
    elseif strcmp(name(end-3:end),'.mat')
        name = load(fullfile(pathname,name));
        name = getfield(name,cell2mat(fieldnames(name)));
    end
    
    for f=1:size(name,1)
        if ~exist(name{f},'file')
            errordlg(sprintf('%s \n file not found',name{f}));
            return
        end
    end
else
    disp('selection aborded')
    return
end


%% now collect data
Ns = length(name);
channel_vector = NaN(Ns,1);
urchan_vector = NaN(Ns,1);

for i=Ns:-1:1
    tmp = load(name{i});
    if isfield(tmp,'LIMO')
        tmp = load(fullfile(tmp.LIMO.dir,'R2.mat'));
    end
    tmp = getfield(tmp,cell2mat(fieldnames(tmp)));
    
    try
        if numel(size(tmp)) == 4
            data{i} = squeeze(tmp(:,:,:,end-1)); % end-1 because R2 dim is R2,F,p
            % and condition/covariates dim are F/p
            data_size{i} = size(data{i});
            index = find(data{i} == max(data{i}(:)));
            [channel_vector(i),~,~]=ind2sub(data_size{i},index);
        else
            data{i} = squeeze(tmp(:,:,end-1));
            data_size{i} = size(data{i});
            index = find(data{i} == max(data{i}(:)));
            [channel_vector(i),~]=ind2sub(data_size{i},index);
        end
        
        LIMO = load([fileparts(name{i}) filesep 'LIMO.mat']); LIMO =LIMO.LIMO;
        if ~isempty(LIMO.data.chanlocs(channel_vector(i)).urchan)
            urchan_vector(i) = LIMO.data.chanlocs(channel_vector(i)).urchan;
        end
        fprintf('subject %g analysed \n',i);
        
    catch bug
        error('file error with subject %s\n%s',name{i},bug.message);
    end
end

% usually no output just write it down on the disk for latter use
if nargout == 0
    cd(current_dir)
    name = cell2mat(inputdlg('Save electrode vector as','Name'));
    if isempty(name)
        return
    else
        save (name,'channel_vector')
        assignin('base',name,channel_vector)
        save ([name '_urchan'],'urchan_vector')
        assignin('base',[name '_urchan'],urchan_vector)
    end
end


%% do the map
clear data
if sum(isnan(channel_vector)) == 0 && nargin ==0 || ...
        sum(isnan(channel_vector)) == 0 && nargin ==2
    
    if nargin <= 1
        [p,f,filt]=uigetfile('*.mat','load expected chanlocs to check channel positions');
        if filt == 0
            return
        else
            expected_chanlocs = load([f filesep p]);
        end
    else
        expected_chanlocs = varargin{2};
        if ischar(expected_chanlocs)
            expected_chanlocs = load(expected_chanlocs);
        end
    end
    
    structnames = fieldnames(expected_chanlocs);
    if ~isempty(structnames)
        expected_chanlocs = expected_chanlocs.(cell2mat(structnames(structfun(@isstruct,expected_chanlocs))));
    end
    data = zeros(1,length(expected_chanlocs));
    
    for S=1:Ns
        data(channel_vector(S)) = data(channel_vector(S))+1;
    end
    
    % create the frequency map
    figure('Color','w','NumberTitle','off','Name','limo_tools: best electrode frequency map')
    [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
        topoplot(data,expected_chanlocs,'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on','conv','on');
    freqmap = grid_or_val(end:-1:1,:); % reverse row order to get back of the head at the bottom of the image
    freqmap(34,67)=NaN;freqmap(67,34)=NaN;freqmap(34,1)=NaN;freqmap(1,34)=NaN;
    if min(freqmap(:))<0
        freqmap = freqmap + abs(min(freqmap(:)));
    end
    
    % make the figure
    imagesc(freqmap,[0 max(freqmap(:))])
    axis tight;axis square;axis off
    colormap(limo_color_images(freqmap));
    cd(current_dir)
end


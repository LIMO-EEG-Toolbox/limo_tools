function [electrode_vector,urchan_vector] = limo_best_electrodes(varargin)

% This function finds the electrode with the maximum F value in each subject.
% The function works on files R2.mat, Condition_effect.mat, Continuous.mat. 
% The function returns a map of frequency showing how often an electrode is
% selected across subjects. 
%
% FORMAT [electrode_vector,urchan_vector] = limo_best_electrodes(expected_chanlocs)
%
% INPUT if empty user is prompted overwise an expected chanlocs (list of
% fullcap) must be provided -- if the format is wrong or aborded during
% selection the map is not computed ; but the function outputs will
%
% OUTPUTS electrode_vector the indices of which electrodes had the strongest F values
%         urchan_vector the value read from the urchan field (if present)
%         a frequency map is also presented as graphical output
%
% -----------------------------
%  Copyright (C) LIMO Team 2014

% Cyril Pernet 15 July 2010
% GAR, 15 August 2010: updated description, uigetfile prompt, frequency figure
% Cyril Pernet May 2014: revamp + update for time-frequency 

%% file selection

current_dir = pwd;
[name,pathname,FilterIndex]=uigetfile({'*.mat','MAT-files (*.mat)'; '*.txt','Text (*.txt)'}, ...
    'Pick a list of result files (e.g. R2)');

if FilterIndex ~= 0
    
    if strcmp(name(end-3:end),'.txt')
        name = importdata(name);
    elseif strcmp(name(end-3:end),'.mat')
        name = load([pathname name]);
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
electrode_vector = NaN(Ns,1);
urchan_vector = NaN(Ns,1);

for i=1:Ns
    tmp = load(name{i});
    tmp = getfield(tmp,cell2mat(fieldnames(tmp)));
    try
        if numel(size(tmp)) == 4
            data{i} = squeeze(tmp(:,:,:,end-1)); % end-1 because R2 dim is R2,F,p
            % and condition/covariates dim are F/p
            data_size{i} = size(data{i});
            index = find(data{i} == max(data{i}(:)));
            [electrode_vector(i),~,~]=ind2sub(data_size{i},index);
        else
            data{i} = squeeze(tmp(:,:,end-1));
            data_size{i} = size(data{i});
            index = find(data{i} == max(data{i}(:)));
            [electrode_vector(i),~]=ind2sub(data_size{i},index);
        end
        
        try
            load([fileparts(name{i}) filesep 'LIMO.mat'])
            urchan_vector(i) = LIMO.data.chanlocs(electrode_vector(i)).urchan;
        catch nourchan
            fprintf('can''t read data.chanlocs.urchan subject %g \n',i)
        end
        fprintf('subject %g analysed \n',i);
        
    catch ME
        message = sprintf('file error, the map of subject %s is not recognized',name{i});
        error([message])
    end
end

% usually no output just write it down on the disk for latter use
if nargout == 0
    cd(current_dir)
    name = cell2mat(inputdlg('Save electrode vector as','Name'));
    if isempty(name)
        return
    else
        save ([name],'electrode_vector')
        assignin('base',[name],electrode_vector)
        save ([name '_urchan'],'urchan_vector')
        assignin('base',[name '_urchan'],urchan_vector)
    end
end


%% do the map
clear data
if sum(isnan(electrode_vector)) == 0
    
    if nargin == 0
        [p,f,filt]=uigetfile('load expected chanlocs');
        if filt == 0
            return
        else
            load([f filesep p]);
        end
    else
        expected_chanlocs = varargin{1};
    end
    data = zeros(1,length(expected_chanlocs));
    
    for S=1:Ns
        data(electrode_vector(S)) = data(electrode_vector(S))+1;
    end
    
    % create the frequency map
    figure('Color','w','NumberTitle','off','Name','limo_tools: best electrode frequency map')
    [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
        topoplot( data,expected_chanlocs,'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on','conv','on');
    freqmap = grid_or_val(end:-1:1,:); % reverse row order to get back of the head at the bottom of the image
    freqmap(34,67)=NaN;freqmap(67,34)=NaN;freqmap(34,1)=NaN;freqmap(1,34)=NaN;
    freqmap(isnan(freqmap))=max(freqmap(:))+1;
    
    % make the figure
    imagesc(freqmap,[0 max(freqmap(:))])
    axis tight;axis square;axis off
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];cc(end,:)=[1 1 1];colormap(cc);
    cd(current_dir)
end


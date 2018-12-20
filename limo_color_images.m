function cm = limo_color_images(scale,LIMO,varargin)

%% limo color map
% 
% FORNAT cc = limo_color_images_(scale,LIMO,'map','map_name','limits',[min max])
%
% INPUT scale if the data to show (typically stat map masked by NaN)
%       LIMO is the LIMO.mat associated to the data
%       map_name is either 'cubehelix' (default) or any map from limo_eeg\external\color_maps
%       [min max] are the values to set the colormap to
%
% OUPUT cm is your colormap
%
% References:
% Green, D. A., 2011, `A colour scheme for the display of astronomical intensity images', Bulletin of the Astronomical Society of India, 39, 289.
% Pernet & Madan. Data visualization for inference in tomographic brain imaging. https://psyarxiv.com/egc6q.

% defaults
map_name = [];
limits = [];
color_path = [fileparts(which('limo_eeg')) filesep 'external' filesep 'color_maps' filesep];

% possible inputs
for n=1:length(varargin)
    if strcmpi(varargin{n},'map')
        if strcmp(varargin{n+1},'cubehelix') || exist([color_path varargin{n+1}],'file')
            map_name = varargin{n+1};
        end
    elseif strcmpi(varargin{n},'limits') || strcmpi(varargin{n},'limit')
        limits = [];
    end
end

% get the map
if isempty(map_name)
    if min(scale(:)) >= 0
        cm = cubehelixmap('increase',64);
    elseif max(scale(:)) <= 0
        cm = cubehelixmap('decrease',64);
    else
        cm = load([color_path 'diverging_bwr.mat']);
        cm = getfield(cm,cell2mat(fieldnames(cm)));
    end
else
    cm = load([color_path map_name]);
    if isstruct(cm)
        cm = getfield(cm,cell2mat(fieldnames(cm)));
    end
end

% set gray background
if sum(isnan(scale(:))) ~= 0
    cm(1,:)=[.9 .9 .9]; % set NaNs to gray
end

colormap(cm);
if ~isempty(limits)
    caxis(limits);
end

set(gca,'XMinorTick','on','LineWidth',2)
try
    set(gca,'YTick',1:length(LIMO.data.expected_chanlocs));
catch ME
    set(gca,'YTick',1:length(LIMO.data.chanlocs));
end

if strcmp(LIMO.Analysis,'Time')
    xlabel('Time in ms','FontSize',10)
elseif strcmp(LIMO.Analysis,'Frequency')
    xlabel('Frequency in Hz','FontSize',10)
end

if strcmp(LIMO.Type,'Components')
    if size(scale,1) == 1
        label_electrodes = ' ';
        ylabel('optimized component','FontSize',10);
    else
        ylabel('Components','FontSize',10);
        for i=1:size(scale,1)
            label_electrodes{i} = i;
        end
    end
else
    if size(scale,1) == 1
        label_electrodes = ' ';
        ylabel('optimized electrode','FontSize',10);
    else
        ylabel('Electrodes','FontSize',10);
        for i = 1:length(LIMO.data.chanlocs)
            if LIMO.Level == 2
                try
                    label_electrodes{i} = LIMO.data.expected_chanlocs(i).labels;
                catch
                    label_electrodes{i} = i;
                end
            else
                try
                    label_electrodes{i} = LIMO.data.chanlocs(i).labels;
                catch
                    label_electrodes{i} = i;
                end
            end
        end
    end
end
set(gca,'YTickLabel', label_electrodes);
end

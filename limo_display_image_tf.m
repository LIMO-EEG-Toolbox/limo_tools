function varargout = limo_display_image_tf(varargin)
% limo_display_results_tf: interactive GUI to display time*freq results
%
%      INPUT limo_display_results_tf(LIMO,toplot,mask,title)
%       1 - LIMO struct
%       2 - 3D matrix of values to plot, dim elec x freqs x time-bins
%       3 - 3D matrix of significant cells
%       4 - title (from limo_stat_values_tf)
%
% ----------------------------------------------------------------------
%  Copyright (C) LIMO Team 2020

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @limo_display_results_tf_OpeningFcn, ...
    'gui_OutputFcn',  @limo_display_results_tf_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before limo_display_results_tf is made visible.
function limo_display_results_tf_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output     = hObject;
handles.LIMO       = varargin{1};
handles.data3d     = varargin{2};
handles.mask       = varargin{3};
scale              = handles.data3d.*single(handles.mask>0);
scale(scale==0)    = NaN;
handles.cc         = limo_color_images(scale);
handles.scale      = scale;
handles.title      = varargin{4};
handles.plot_sel   = 1;

% get axes right away
if isfield(handles.LIMO.data,'tf_freqs')
    handles.freqs_here = round(handles.LIMO.data.tf_freqs);
else
    handles.freqs_here = round(linspace(handles.LIMO.data.lowf,handles.LIMO.data.highf,size(handles.data3d,2)));
end
if isfield(handles.LIMO.data,'tf_times')
    handles.times_here = round(handles.LIMO.data.tf_times);
else
    handles.times_here = round(linspace(handles.LIMO.data.start,handles.LIMO.data.end,size(handles.data3d,3)));
end

% for each cluster, get start/end/max value
% if unthresholded, uncorrected, tfce or max = mask is made up of ones
handles.n_cluster     = max(handles.mask(:));
handles.cluster_start = NaN(2,handles.n_cluster); % start of each cluster
handles.cluster_end   = NaN(2,handles.n_cluster); % end of each cluster
handles.cluster_maxv  = NaN(1,handles.n_cluster); % max value for each cluster
handles.cluster_maxe  = NaN(1,handles.n_cluster); % channel location of the max value of each cluster
handles.cluster_maxf  = NaN(1,handles.n_cluster); % frame location of the max value of each cluster
handles.cluster_maxt  = NaN(1,handles.n_cluster); % frame location of the max value of each cluster
for c=1:handles.n_cluster
    tmp                               = handles.data3d.*(handles.mask==c);
    sigframes                         = squeeze(sum(tmp,1));
    handles.cluster_start(:,c)        = [find(sum(sigframes,2),1,'first') find(sum(sigframes,1),1,'first')];
    handles.cluster_end(:,c)          = [find(sum(sigframes,2),1,'last')  find(sum(sigframes,1),1,'last')];
    V                                 = max(tmp(:));
    handles.cluster_maxv(c)           = V(1);
    [e,f,t]                           = ind2sub(size(tmp),find(tmp==V(1)));
    handles.cluster_maxe(c)           = e(1);
    handles.cluster_maxf(c)           = f(1);
    handles.cluster_maxt(c)           = t(1);
end
clear varargin scale n_cluster

% Find max values, save idx and value
if (size(handles.data3d,4)) == 3
    dim = inputdlg('which dimension to use','4D input - plotting option');
    handles.data3d = handles.data3d(:,:,:,dim);
    handles.scale  = handles.scale(:,:,:,dim);
end
handles.maxv       = max(handles.scale(:));
handles.maxvi      = find(handles.scale == handles.maxv);

if length(handles.maxvi) ~= 1
    handles.maxvi = handles.maxvi(1);
end
[handles.maxe, handles.maxf, handles.maxt] = ind2sub(size(handles.scale), handles.maxvi);
handles.slider_sel                         = handles.maxt;
plot_data.freqs_here                       = handles.LIMO.data.tf_freqs;
plot_data.times_here                       = handles.LIMO.data.tf_times;

guidata(hObject, plot_data);
guidata(hObject, handles);

% ----------------------------------------------------------
% This sets up the initial plot - only do when invisible
% ----------------------------------------------------------
% show the channel/compoment * freq map at the time values is maximal

if strcmp(get(hObject,'Visible'),'off')
    
    % stat value
    plot_statvalues(handles,squeeze(handles.data3d(handles.maxe,:,handles.maxt)),'Frequency')
    
    % topoplot
    plot_topography(handles,squeeze(handles.data3d(:,handles.maxf,handles.maxt)),handles.maxf,handles.maxt)
    
    % time/freq map
    plot_tfmap(handles,flipud(squeeze(handles.data3d(handles.maxe,:,:))),handles.maxe);
    
    % main display
    plot_main(handles,squeeze(handles.scale(:,:,handles.maxt)),'Frequency',handles.times_here(handles.maxt));
    
    % channel*time
    plot_chantime(handles,squeeze(handles.scale(:,handles.maxf,:)),round(handles.freqs_here(handles.maxf)))
    
    % report all clusters in command window
    if handles.n_cluster > 1
        for c=1:handles.n_cluster
            fprintf('cluster %g starts at %gms %gHz, ends at %gms %gHz, max %g @ %gms %gHz channel %s \n', c, ...
                round(handles.times_here(handles.cluster_start(2,c))),round(handles.freqs_here(handles.cluster_start(1,c))),...
                round(handles.times_here(handles.cluster_end(2,c))),round(handles.freqs_here(handles.cluster_end(1,c))),...
                handles.cluster_maxv(c), handles.times_here(handles.cluster_maxt(c)), handles.freqs_here(handles.cluster_maxf(c)),...
                handles.LIMO.data.chanlocs(handles.cluster_maxe(c)).labels);
        end
    else % no clusters
        fprintf('1st significant frame at %gms %gHz, last signifiant frame at %gms %gHz, max %g @ %gms %gHz channel %s \n', ...
            round(handles.times_here(handles.cluster_start(2,c))),round(handles.freqs_here(handles.cluster_start(1,c))),...
            round(handles.times_here(handles.cluster_end(2,c))),round(handles.freqs_here(handles.cluster_end(1,c))),...
            handles.cluster_maxv(c), handles.times_here(handles.cluster_maxt(c)), handles.freqs_here(handles.cluster_maxf(c)),...
            handles.LIMO.data.chanlocs(handles.cluster_maxe(c)).labels);
    end
end

function varargout = limo_display_results_tf_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

%% interactive part

% ------------------------------------------------------------------
%      UPDATE: swith between Electrode * Frequencies (slide in time)
%              and Electrode * Time (slide in Frequencies)
% ------------------------------------------------------------------
function Main_display_CreateFcn(hObject, eventdata, handles)

function topoplot_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function pop_up_dimensions_CreateFcn(hObject, eventdata, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Electrodes x Frequencies', 'Electrodes x Times'});

% --- Executes on selection change in pop_up_dimensions.
function pop_up_dimensions_Callback(hObject, eventdata, handles)

popup_sel_index = get(hObject,'Value');
switch popup_sel_index
    case 1
        
        % stat value
        plot_statvalues(handles,squeeze(handles.data3d(handles.maxe,:,handles.maxt)),'Frequency')
        
        % topoplot
        plot_topography(handles,squeeze(handles.data3d(:,handles.maxf,handles.maxt)),handles.maxf,handles.maxt)
        
        % time/freq map
        plot_tfmap(handles,flipud(squeeze(handles.data3d(handles.maxe,:,:))),handles.maxe);
        
        % main display
        plot_main(handles,squeeze(handles.scale(:,:,handles.maxt)),'Frequency',handles.times_here(handles.maxt));
        
        % channel*time
        plot_chantime(handles,squeeze(handles.scale(:,handles.maxf,:)),round(handles.freqs_here(handles.maxf)))
        
    case 2
        
        % stat value
        plot_statvalues(handles,squeeze(handles.data3d(handles.maxe,handles.maxf,:)),'Time')
        
        % topoplot
        plot_topography(handles,squeeze(handles.data3d(:,handles.maxf,handles.maxt)),handles.maxf,handles.maxt)
        
        % time/freq map
        plot_tfmap(handles,flipud(squeeze(handles.data3d(handles.maxe,:,:))),handles.maxe);
        
        % main display
        plot_main(handles,squeeze(handles.scale(:,handles.maxf,:)),'Time',handles.freqs_here(handles.maxf));
        
        % channel*freq
        plot_chanfreq(handles,squeeze(handles.scale(:,:,handles.maxt)),round(handles.times_here(handles.maxt)))
end
guidata(hObject, handles);

%  -----------------------------------------------------------
%       SLIDER -- slide in time or in freq
% ------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function slider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_Callback(hObject, eventdata, handles)
plot_data       = guidata(hObject);
popup_sel_index = get(handles.pop_up_dimensions, 'Value');

if popup_sel_index==1
    slider_sel = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.times_here)*slider_sel)); % Scale to get ints 1:5 out
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    elseif slider_sel > size(handles.scale,3)
        slider_sel = size(handles.scale,3);
    end
    D              = squeeze(handles.scale(:,:,slider_sel));
    [maxe,maxf]    = ind2sub(size(D),find(D == max(D(:))));
    if isempty(maxe)
        tmp            = squeeze(handles.data3d(:,:,slider_sel));
        [maxe,maxf]    = ind2sub(size(D),find(tmp == max(tmp(:))));
    end
    
    % stat value
    plot_statvalues(handles,squeeze(handles.data3d(maxe,:,slider_sel)),'Frequency')
    
    % topoplot
    plot_topography(handles,squeeze(handles.data3d(:,maxf,slider_sel)),maxf,slider_sel)
    
    % time/freq map
    plot_tfmap(handles,flipud(squeeze(handles.data3d(maxe,:,:))),maxe);
    
    if sum(isnan(D(:))) ~= numel(D)
        % main display
        plot_main(handles,D,'Frequency',handles.times_here(slider_sel))
        
        % channel*time
        plot_chantime(handles,squeeze(handles.scale(:,maxf,:)),round(handles.freqs_here(maxf)))
    else
        % main display
        plot_main(handles,zeros(size(D,1),size(D,2)),'Frequency',handles.times_here(slider_sel))
        
        % channel*time
        plot_chantime(handles,zeros(size(D,1),size(D,2)),round(handles.freqs_here(maxf)))
        
    end
elseif popup_sel_index==2
    slider_sel2 = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.freqs_here)*slider_sel2)); % Scale slider to correct ints
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    elseif slider_sel > size(handles.scale,2)
        slider_sel = size(handles.scale,2);
    end
    D              = squeeze(handles.scale(:,slider_sel,:));
    [maxe,maxt]    = ind2sub(size(D),find(D == max(D(:))));
    if isempty(maxe)
        tmp            = squeeze(handles.data3d(:,slider_sel,:));
        [maxe,maxf]    = ind2sub(size(D),find(tmp == max(tmp(:))));
    end
    
    % stat value
    plot_statvalues(handles,squeeze(handles.data3d(maxe,slider_sel,:)),'Time')
    
    % topoplot
    plot_topography(handles,squeeze(handles.data3d(:,slider_sel,maxt)),slider_sel,maxt)
    
    % time/freq map
    plot_tfmap(handles,flipud(squeeze(handles.data3d(maxe,:,:))),maxe);
    
    if sum(isnan(D(:))) ~= numel(D)
        % main display
        plot_main(handles,D,'Time',handles.freqs_here(slider_sel))
        
        % channel*freq
        plot_chanfreq(handles,squeeze(handles.scale(:,:,maxt)),round(handles.times_here(maxt)))
    else
        % main display
        plot_main(handles,zeros(size(D,1),size(D,3)),'Time',handles.freqs_here(slider_sel))
        
        % channel*freq
        plot_chanfreq(handles,zeros(size(D,1),size(D,3)),round(handles.times_here(maxt)))
    end
end

handles.slider_sel=slider_sel;
guidata(hObject, handles);


% -------------------------------------------------
%                     MOUSE INPUT
% -------------------------------------------------

% --- Executes on button press in mouse_input.
function mouse_input_Callback(hObject, eventdata, handles)

plot_data       = guidata(hObject);
popup_sel_index = get(handles.pop_up_dimensions, 'Value');
button          = 1;

while button == 1
    [x,y,button] = ginput(1);
    clickedAx   = gca;
    if clickedAx == handles.Main_display
        if x < 1; x=1; end
        if y < 1; y=1; end
    else
        break
    end
    
    if popup_sel_index == 1  % if showing elec x freq
        
        if x < min(plot_data.freqs_here); x=min(plot_data.freqs_here); end
        if x > max(plot_data.freqs_here); x=max(plot_data.freqs_here); end
        freq = round(x); channel = round(y);
        [~,freq_position] = min(abs(plot_data.freqs_here-freq));
        
        % stat value
        if handles.slider_sel > handles.maxt
            handles.slider_sel = handles.maxt;
        end
        plot_statvalues(handles,squeeze(handles.data3d(channel,:,handles.slider_sel)),'Frequency')
        
        % topoplot
        plot_topography(handles,squeeze(handles.data3d(:,freq_position,handles.slider_sel)),freq_position,handles.slider_sel)
        
        % time/freq map
        plot_tfmap(handles,flipud(squeeze(handles.data3d(channel,:,:))),channel);
        
        % reset colormap
        axes(handles.Main_display);
        colormap(gca, handles.cc);
        
        % channel*time
        fprintf('frequency selected: %g Hz\n',freq)
        plot_chantime(handles,squeeze(handles.scale(:,freq_position,:)),freq)
        
        try
            p_values = evalin('base','p_values');
            if ~isnan(p_values(channel,freq,handles.slider_sel))
                fprintf('Stat value: %g, p_value %g \n',handles.data3d(channel,freq,handles.slider_sel),p_values(channel,freq,handles.slider_sel));
            end
        catch pvalerror
            fprintf('couldn''t figure the p value?? %s \n',pvalerror.message)
        end
        
    elseif popup_sel_index == 2  % if showing elec x times on main
        
        if x < min(plot_data.times_here); x=min(plot_data.times_here); end
        if x > max(plot_data.times_here); x=max(plot_data.times_here); end
        time = round(x); channel = round(y);
        [~,time_position] = min(abs(plot_data.times_here-time));
        
        % stat value
        if handles.slider_sel > handles.maxf
            handles.slider_sel = handles.maxf;
        end
        plot_statvalues(handles,squeeze(handles.data3d(channel,handles.slider_sel,:)),'Time')
        
        % topoplot
        plot_topography(handles,squeeze(handles.data3d(:,handles.slider_sel,time_position)),handles.slider_sel,time_position)
        
        % time/freq map
        plot_tfmap(handles,flipud(squeeze(handles.data3d(channel,:,:))),channel);
        
        % reset colormap
        axes(handles.Main_display);
        colormap(gca, handles.cc);
        
        % channel*freq
        fprintf('time selected: %g ms\n',time)
        plot_chanfreq(handles,squeeze(handles.scale(:,:,time_position)),time)
        
        try
            p_values = evalin('base','p_values');
            if ~isnan(p_values(channel,handles.slider_sel,time))
                fprintf('Stat value: %g, p_value %g \n',handles.data3d(channel,handles.slider_sel,time),p_values(channel,handles.slider_sel,time));
            end
        catch pvalerror
            fprintf('couldn''t figure the p value?? %s \n',pvalerror.message)
        end
    end
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%% MAKE A MOVIE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in gifbutton.
function gifbutton_Callback(hObject, eventdata, handles)

vidObj = VideoWriter(handles.title);
open(vidObj);
h=figure('Name',handles.title);
set(gcf,'Color','w','InvertHardCopy','off', 'units','normalized', 'outerposition',[0 0 1 1]);
popup_sel_index = get(handles.pop_up_dimensions, 'Value');
if popup_sel_index == 1  % if showing elec x freq
    for t=1:size(handles.scale,3)
        imagesc(squeeze(handles.scale(:,:,t)));
        colormap(gca, handles.cc); img_prop = get(gca);
        newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.XTick)));
        Xlabels = handles.freqs_here(newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
        if handles.LIMO.Level == 1
            Ylabels  = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
            newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
            Ylabels  = Ylabels(newyticks);
        else
            if isempty(handles.LIMO.design.electrode)
                if isfield(handles.LIMO.data,'chanlocs')
                    Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
                else
                    Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
                end
                newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
                Ylabels  = Ylabels(newyticks);
            else
                ylabel('optimized electrode','fontsize',10);
            end
        end
        if exist('Ylabels','var')
            set(gca,'YTick',newyticks);
            set(gca,'YTickLabel', Ylabels);
        end
        title(['Channels x Freq @ ' num2str(round(handles.times_here(t))) ' ms'],'VerticalAlignment','bottom');
        xlabel('Freq (Hz)','fontsize',10); axis tight
        drawnow
        pause(0.25)
    end
else
    for f=1:size(handles.scale,2)
        imagesc(squeeze(handles.scale(:,f,:)));
        colormap(gca, handles.cc); img_prop = get(gca);
        newxticks = round(linspace(1,length(handles.times_here),length(img_prop.XTick)));
        Xlabels = handles.times_here(newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
        if handles.LIMO.Level == 1
            Ylabels  = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
            newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
            Ylabels  = Ylabels(newyticks);
        else
            if isempty(handles.LIMO.design.electrode)
                if isfield(handles.LIMO.data,'chanlocs')
                    Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
                else
                    Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
                end
                newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
                Ylabels  = Ylabels(newyticks);
            else
                ylabel('optimized electrode','fontsize',10);
            end
        end
        if exist('Ylabels','var')
            set(gca,'YTick',newyticks);
            set(gca,'YTickLabel', Ylabels);
        end
        title(['Channels x Times @ ' num2str(round(handles.freqs_here(f))) ' Hz'],'VerticalAlignment','bottom');
        xlabel('Freq (Hz)','fontsize',10); axis tight
        drawnow
        pause(0.25)
    end
end
warning off
close(vidObj); close(h)
warning on
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%% MENU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)

[file,p] = uigetfile('*.mat','select effect file');
if ~exist(fullfile(p,'LIMO.mat'),'file')
    error('there is no LIMO.mat associated with this file')
end

if ~isequal(file, 0)
    data3d = load(file); data3d = data3d.(cell2mat(fieldnames(data3d)));
    LIMO = load(fullfile(p,'LIMO.mat')); LIMO = LIMO.(cell2mat(fieldnames(LIMO)));
    do_not_update = 0;
    if contains(file,'R2','IgnoreCase',true)
        data3d     = squeeze(data3d(:,:,:,2));
    elseif contains(file,'one_sample','IgnoreCase',true)
        data3d     = squeeze(data3d(:,:,:,4));
    else
        try
            data3d     = squeeze(data3d(:,:,:,1));
        catch nosqueeze
            do_not_update = 1;
            errordlg(sprintf('file not supported \n%s',nosqueeze.message))
        end
    end
    
    % ---------------------
    if do_not_update == 0
        clc; uiresume
        guidata(hObject, handles);
        delete(handles.figure1)
        file(strfind(file,'_')) = ' ';
        limo_display_image_tf(LIMO,data3d,ones(size(data3d)),file)
    end
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)

% saveas(handles.figure1,cell2mat(inputdlg('save figure as: ')))
name = cell2mat(inputdlg('save figure as: '));
fig_id = findall(0,'type','figure');
set(fig_id, 'PaperPositionMode', 'auto');
print(fig_id,'-depsc2',[name '.eps'])
guidata(hObject, handles);

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
clc; uiresume
guidata(hObject, handles);
delete(handles.figure1)

%% routines
function plot_statvalues(handles,courseplot,dim)
axes(handles.tf_course_plot); cla
[M,position] = max(courseplot(:));
if strcmpi(dim,'Time')
    plot(handles.times_here,courseplot,'LineWidth',3);
    title(sprintf('Max stat %g @ %gms',M,handles.times_here(position)),'VerticalAlignment','bottom');
elseif strcmpi(dim,'Frequency')
    plot(handles.freqs_here,courseplot,'LineWidth',3);
    title(sprintf('Max stat %g @ %gHz',M,handles.freqs_here(position)),'VerticalAlignment','bottom');
end
xlabel(dim,'fontsize',10,'VerticalAlignment','top');
ylabel('Stat value','fontsize',10);
box on; grid on; axis tight

function plot_topography(handles,topomap,freq_value,time_value)
axes(handles.topoplot); cla
topoplot(topomap,handles.LIMO.data.chanlocs);
colormap(gca, handles.cc(2:end,:));
mytitle = sprintf('topoplot @%gms & %gHz', handles.times_here(time_value), round(handles.freqs_here(freq_value)));
title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')

function plot_tfmap(handles,tf_map,channel_value)

axes(handles.sub_display2); cla
imagesc(handles.times_here,handles.freqs_here,tf_map);
colormap(gca, limo_color_images(tf_map)); img_prop = get(gca);
newticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
Ylabels  = fliplr(round(handles.freqs_here(newticks))); set(gca,'YTickLabel', split(string(Ylabels)))
title(['Time-Frequency @ channel ' num2str(handles.LIMO.data.chanlocs(channel_value).labels)]);
xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequency','fontsize',10);

function plot_main(handles,scaled_data,dim,dimvalue)

axes(handles.Main_display); cla
if ~isempty(handles.LIMO.design.electrode)
    chan_vect = 1;
else
    chan_vect = 1:length(handles.LIMO.data.chanlocs);
end

if strcmpi(dim,'Frequency')
    imagesc(handles.freqs_here,chan_vect,scaled_data);
    title(sprintf('%s @ %g ms',regexprep(handles.title,'\n+',''), dimvalue),'fontsize',12,'VerticalAlignment','bottom');
    xlabel('Frequency bins (Hz)','VerticalAlignment','top','fontsize',10);
elseif strcmpi(dim,'Time')
    imagesc(handles.times_here,chan_vect,scaled_data);
    title(sprintf('%s @ %gHz',regexprep(handles.title,'\n+',''), dimvalue),'fontsize',12,'VerticalAlignment','bottom');
    xlabel('Time (ms)','VerticalAlignment','top','fontsize',10);
end
if sum(scaled_data(:)) ~= 0
    colormap(gca, handles.cc);
else
    colormap(gca,[0.9 0.9 0.9]);
end
set(gca,'LineWidth',2);
set_Ylabels(handles,get(gca));

function plot_chantime(handles,map,freqvalue)

axes(handles.sub_display1); cla
if ~isempty(handles.LIMO.design.electrode)
    chan_vect = 1;
else
    chan_vect = 1:length(handles.LIMO.data.chanlocs);
end
imagesc(handles.times_here,chan_vect,map);
title(['Channels x Time @ ' num2str(freqvalue) ' Hz'],'VerticalAlignment','bottom');
xlabel('Time (ms)','fontsize',10); axis tight
if sum(map(:)) ~= 0
    colormap(gca, handles.cc);
else
    colormap(gca,[0.9 0.9 0.9]);
end
set_Ylabels(handles,get(gca));

function plot_chanfreq(handles,map,timevalue)

axes(handles.sub_display1); cla
if ~isempty(handles.LIMO.design.electrode)
    chan_vect = 1;
else
    chan_vect = 1:length(handles.LIMO.data.chanlocs);
end

imagesc(handles.freqs_here,chan_vect,map);
title(['Channels x Frequency @ ' num2str(timevalue) 'ms'],'VerticalAlignment','bottom');
xlabel('Frequencies (Hz)','fontsize',10); axis tight
if sum(map(:)) ~= 0
    colormap(gca, handles.cc);
else
    colormap(gca,[0.9 0.9 0.9]);
end
set_Ylabels(handles,get(gca));

function set_Ylabels(handles,img_prop)

if handles.LIMO.Level == 1
    Ylabels   = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
    newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
    Ylabels   = Ylabels(newyticks);
else
    if isempty(handles.LIMO.design.electrode)
        if isfield(handles.LIMO.data,'chanlocs')
            Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
        else
            Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
        end
        newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
        Ylabels  = Ylabels(newyticks);
    else
        ylabel('optimized electrode','fontsize',10);
    end
end

if exist('Ylabels','var')
    set(gca,'YTick',newyticks);
    set(gca,'YTickLabel', Ylabels);
end


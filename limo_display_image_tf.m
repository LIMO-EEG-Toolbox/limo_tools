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
    axes(handles.tf_course_plot); 
    D = squeeze(handles.data3d(handles.maxe,:,handles.maxt));
    plot(handles.freqs_here,D,'LineWidth',3); M = max(D(:));
    ylabel('Stat value','fontsize',10); xlabel('Frequencies','fontsize',10,'VerticalAlignment','top');
    title(sprintf('Max stat %g @ %gHz',M,handles.freqs_here(handles.maxf)),'VerticalAlignment','bottom'); grid on; axis tight

    % topoplot
    axes(handles.topoplot);
    topoplot(squeeze(handles.scale(:,handles.maxf,handles.maxt)),handles.LIMO.data.chanlocs);
    colormap(gca, handles.cc(2:end,:));
    mytitle = sprintf('topoplot @%gms & %gHz', handles.times_here(handles.maxt), round(handles.freqs_here(handles.maxf)));
    title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')

    % main display
    axes(handles.Main_display);   
    imagesc(squeeze(handles.scale(:,:,handles.maxt)));
    colormap(gca, handles.cc); img_prop = get(gca); set(gca,'LineWidth',2);
    title(sprintf('%s @ %g ms',regexprep(handles.title,'\n+',''), round(handles.times_here(handles.maxt))),'fontsize',12,'VerticalAlignment','bottom');
    newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.XTick)));
    Xlabels = handles.freqs_here(newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    xlabel('Frequency bins (Hz)','VerticalAlignment','top','fontsize',10);
    if handles.LIMO.Level == 1
        Ylabels  = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
        newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        Ylabels  = Ylabels(newyticks);
    else
        if isempty(handles.LIMO.design.electrode)
            if isfield(handles.LIMO.data,'chanlocs')
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
            else
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
            end
            newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
            Ylabels  = Ylabels(newyticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    
    % channel * time
    axes(handles.sub_display1); cla;
    imagesc(squeeze(handles.scale(:,handles.maxf,:)));
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
            newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
            Ylabels  = Ylabels(newticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    title(['Channels x Times @ ' num2str(round(handles.freqs_here(handles.maxf))) ' Hz'],'VerticalAlignment','bottom');
    xlabel('Time (ms)','fontsize',10); axis tight

    % time * frequency 
    axes(handles.sub_display2); cla;
    imagesc(squeeze(handles.scale(handles.maxe,:,:)));
    colormap(gca, handles.cc); img_prop = get(gca);
    set(gca,'XTick',newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    newticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
    Ylabels = fliplr(round(handles.freqs_here(newticks))); set(gca,'YTickLabel', split(string(Ylabels)))
    title(['Frequency x time @ channel ' num2str(handles.LIMO.data.chanlocs(handles.maxe).labels)]);
    xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequencies','fontsize',10);
    
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
function pop_up_dimensions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', {'Electrodes x Frequencies', 'Electrodes x Times'});

% --- Executes on selection change in pop_up_dimensions.
function pop_up_dimensions_Callback(hObject, eventdata, handles)

popup_sel_index = get(hObject,'Value');  
axes(handles.Main_display); cla;
timep = numel(handles.times_here);
freqp = numel(handles.freqs_here);

switch popup_sel_index
    case 1
        
    % stat value
    axes(handles.tf_course_plot); 
    D = squeeze(handles.data3d(handles.maxe,:,handles.maxt));
    plot(handles.freqs_here,D,'LineWidth',3); M = max(D(:));
    ylabel('Stat value','fontsize',10); xlabel('Frequencies','fontsize',10,'VerticalAlignment','top');
    title(sprintf('Max stat %g @ %gHz',M,handles.freqs_here(handles.maxf)),'VerticalAlignment','bottom'); grid on; axis tight

    % topoplot
    axes(handles.topoplot);
    topoplot(squeeze(handles.scale(:,handles.maxf,handles.maxt)),handles.LIMO.data.chanlocs);
    colormap(gca, handles.cc(2:end,:));
    mytitle = sprintf('topoplot @%gms & %gHz', handles.times_here(handles.maxt), round(handles.freqs_here(handles.maxf)));
    title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')

    % main display
    axes(handles.Main_display);   
    imagesc(squeeze(handles.scale(:,:,handles.maxt)));
    colormap(gca, handles.cc); img_prop = get(gca); set(gca,'LineWidth',2);
    title(sprintf('%s @ %g ms',regexprep(handles.title,'\n+',''), round(handles.times_here(handles.maxt))),'fontsize',12,'VerticalAlignment','bottom');
    newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.XTick)));
    Xlabels = handles.freqs_here(newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    xlabel('Frequency bins (Hz)','VerticalAlignment','top','fontsize',10);
    if handles.LIMO.Level == 1
        Ylabels  = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
        newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        Ylabels  = Ylabels(newyticks);
    else
        if isempty(handles.LIMO.design.electrode)
            if isfield(handles.LIMO.data,'chanlocs')
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
            else
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
            end
            newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
            Ylabels  = Ylabels(newyticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    
    % channel * time
    axes(handles.sub_display1); cla;
    imagesc(squeeze(handles.scale(:,handles.maxf,:)));
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
    title(['Channels x Times @ ' num2str(round(handles.freqs_here(handles.maxf))) ' Hz'],'VerticalAlignment','bottom');
    xlabel('Time (ms)','fontsize',10); axis tight

    % time * frequency 
    axes(handles.sub_display2); cla;
    imagesc(squeeze(handles.scale(handles.maxe,:,:)));
    colormap(gca, handles.cc); img_prop = get(gca);
    set(gca,'XTick',newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    newticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
    Ylabels = fliplr(round(handles.freqs_here(newticks))); set(gca,'YTickLabel', split(string(Ylabels)))
    title(['Frequency x time @ channel ' num2str(handles.LIMO.data.chanlocs(handles.maxe).labels)]);
    xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequencies','fontsize',10);
        
    case 2
        
    % stat value
    axes(handles.tf_course_plot); 
    D = squeeze(handles.data3d(handles.maxe,handles.maxf,:));
    plot(handles.times_here,D,'LineWidth',3); M = max(D(:));
    ylabel('Stat value','fontsize',10); xlabel('Time','fontsize',10,'VerticalAlignment','top');
    title(sprintf('Max stat %g @ %gms',M,handles.times_here(handles.maxt)),'VerticalAlignment','bottom'); grid on; axis tight

    % topoplot
    axes(handles.topoplot);
    topoplot(squeeze(handles.scale(:,handles.maxf,handles.maxt)),handles.LIMO.data.chanlocs);
    colormap(gca, handles.cc(2:end,:));
    mytitle = sprintf('topoplot @%gms & %gHz', handles.times_here(handles.maxt), round(handles.freqs_here(handles.maxf)));
    title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')

    % main display
    axes(handles.Main_display);   
    imagesc(squeeze(handles.scale(:,handles.maxf,:)));
    colormap(gca, handles.cc); img_prop = get(gca); set(gca,'LineWidth',2);
    title(sprintf('%s @ %gHz',regexprep(handles.title,'\n+',''), round(handles.freqs_here(handles.maxf))),'fontsize',12,'VerticalAlignment','bottom');
    newxticks = round(linspace(1,length(handles.times_here),length(img_prop.XTick)));
    Xlabels = handles.times_here(newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    xlabel('Time (ms)','VerticalAlignment','top','fontsize',10);
    if handles.LIMO.Level == 1
        Ylabels  = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
        newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        Ylabels  = Ylabels(newyticks);
    else
        if isempty(handles.LIMO.design.electrode)
            if isfield(handles.LIMO.data,'chanlocs')
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
            else
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
            end
            newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
            Ylabels  = Ylabels(newyticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    
    % channel * freq
    axes(handles.sub_display1); cla;
    imagesc(squeeze(handles.scale(:,:,handles.maxt)));
    colormap(gca, handles.cc); img_prop = get(gca);
    newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.XTick)));
    Xlabels = round(handles.freqs_here(newxticks)); set(gca,'XTickLabel', split(string(Xlabels)))
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
            newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
            Ylabels  = Ylabels(newticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    title(['Channels x Frequencies @ ' num2str(round(handles.times_here(handles.maxt))) 'ms'],'VerticalAlignment','bottom');
    xlabel('Frequencies (Hz)','fontsize',10); axis tight

    % time * frequency 
    axes(handles.sub_display2); cla;
    imagesc(squeeze(handles.scale(handles.maxe,:,:)));
    colormap(gca, handles.cc); img_prop = get(gca);
    newxticks = round(linspace(1,length(handles.times_here),length(img_prop.XTick)));
    Xlabels = round(handles.times_here(newxticks)); set(gca,'XTickLabel', split(string(Xlabels)))
    newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
    Ylabels = fliplr(round(handles.freqs_here(newxticks))); set(gca,'YTickLabel', split(string(Ylabels)))
    title(['Frequency x time @ channel ' num2str(handles.LIMO.data.chanlocs(handles.maxe).labels)]);
    xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequencies','fontsize',10);
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
timep           = numel(handles.times_here);
freqp           = numel(handles.freqs_here);

if popup_sel_index==1
    slider_sel = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.times_here)*slider_sel)); % Scale to get ints 1:5 out
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    end
    D = squeeze(handles.scale(:,:,slider_sel));
    [maxe,maxf]=ind2sub(size(D),find(D == max(D(:))));
    plot_data.maxe = maxe;
    plot_data.maxf = maxf;
    plot_data.maxt = slider_sel;
    
    % stat value
    axes(handles.tf_course_plot);
    d = squeeze(handles.data3d(maxe,:,slider_sel)); axis tight
    plot(handles.freqs_here,d,'LineWidth',3); M = max(d(:));
    ylabel('Stat value','fontsize',10); xlabel('Frequencies','fontsize',10,'VerticalAlignment','top');
    title(sprintf('Max stat %g @ %gHz',M,handles.freqs_here(maxf)),'VerticalAlignment','bottom'); grid on; axis tight
    
    % topoplot at max freq
    axes(handles.topoplot); cla
    topoplot(D(:,maxf),handles.LIMO.data.chanlocs);
    colormap(gca, handles.cc(2:end,:));
    mytitle = sprintf('topoplot @%gms & %gHz', handles.times_here(slider_sel), round(handles.freqs_here(maxf)));
    title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')
    
    % show electrode * freq
    axes(handles.Main_display);
    imagesc(D); colormap(gca, handles.cc); img_prop = get(gca); set(gca,'LineWidth',2);
    title(sprintf('%s @ %g ms',regexprep(handles.title,'\n+',''), round(handles.times_here(slider_sel))),'fontsize',12,'VerticalAlignment','bottom');
    newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.XTick)));
    Xlabels = handles.freqs_here(newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    xlabel('Frequency bins (Hz)','VerticalAlignment','top','fontsize',10);
    if handles.LIMO.Level == 1
        Ylabels  = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
        newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        Ylabels  = Ylabels(newyticks);
    else
        if isempty(handles.LIMO.design.electrode)
            if isfield(handles.LIMO.data,'chanlocs')
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
            else
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
            end
            newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
            Ylabels  = Ylabels(newyticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    
    % channel * time
    axes(handles.sub_display1); cla;
    imagesc(squeeze(handles.scale(:,maxf,:)));
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
            newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)));
            Ylabels  = Ylabels(newticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    title(['Channels x Times @ ' num2str(round(handles.freqs_here(maxf))) ' Hz'],'VerticalAlignment','bottom');
    xlabel('Time (ms)','fontsize',10); axis tight
    
    % time * frequency
    axes(handles.sub_display2); cla;
    imagesc(squeeze(handles.scale(maxe,:,:)));
    colormap(gca, handles.cc); img_prop = get(gca);
    set(gca,'XTick',newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    newticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
    Ylabels = fliplr(round(handles.freqs_here(newticks))); set(gca,'YTickLabel', split(string(Ylabels)))
    title(['Frequency x time @ channel ' num2str(handles.LIMO.data.chanlocs(maxe).labels)]);
    xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequencies','fontsize',10);
    
elseif popup_sel_index==2
    slider_sel2 = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.freqs_here)*slider_sel2)); % Scale slider to correct ints
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    end
    D = squeeze(handles.scale(:,slider_sel,:));
    [maxe,maxt]=ind2sub(size(D),find(D == max(D(:))));
    plot_data.maxe = maxe;
    plot_data.maxf = slider_sel;
    plot_data.maxt = maxt;

    % stat value
    axes(handles.tf_course_plot);
    d = squeeze(handles.data3d(maxe,slider_sel,:)); axis tight
    plot(handles.times_here,d,'LineWidth',3); M = max(d(:));
    ylabel('Stat value','fontsize',10); xlabel('Time','fontsize',10,'VerticalAlignment','top');
    title(sprintf('Max stat %g @ %gms',M,handles.times_here(maxt)),'VerticalAlignment','bottom'); grid on; axis tight
    
    % topoplot at max freq
    axes(handles.topoplot); cla
    topoplot(D(:,maxt),handles.LIMO.data.chanlocs);
    colormap(gca, handles.cc(2:end,:));
    mytitle = sprintf('topoplot @%gms & %gHz', round(handles.times_here(maxt)),handles.freqs_here(slider_sel));
    title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')
    
    % show electrode * time
    axes(handles.Main_display);
    imagesc(D); colormap(gca, handles.cc); img_prop = get(gca); set(gca,'LineWidth',2);
    title(sprintf('%s @ %g Hz',regexprep(handles.title,'\n+',''), round(handles.freqs_here(slider_sel))),'fontsize',12,'VerticalAlignment','bottom');
    newxticks = round(linspace(1,length(handles.times_here),length(img_prop.XTick)));
    Xlabels = handles.times_here(newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
    xlabel('Time (ms)','VerticalAlignment','top','fontsize',10);
    if handles.LIMO.Level == 1
        Ylabels  = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
        newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
        Ylabels  = Ylabels(newyticks);
    else
        if isempty(handles.LIMO.design.electrode)
            if isfield(handles.LIMO.data,'chanlocs')
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.chanlocs, 'UniformOutput', false);
            else
                Ylabels = arrayfun(@(x)(x.labels), handles.LIMO.data.expected_chanlocs, 'UniformOutput', false);
            end
            newyticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
            Ylabels  = Ylabels(newyticks);
        else
            ylabel('optimized electrode','fontsize',10);
        end
    end
    if exist('Ylabels','var')
        set(gca,'YTick',newyticks);
        set(gca,'YTickLabel', Ylabels);
    end
    
    % channel * freq
    axes(handles.sub_display1); cla;
    imagesc(squeeze(handles.scale(:,:,maxt)));
    colormap(gca, handles.cc); img_prop = get(gca);
    newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.XTick)));
    Xlabels = round(handles.freqs_here(newxticks)); set(gca,'XTickLabel', split(string(Xlabels)))
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
    title(['Channels x Frequencies @ ' num2str(round(handles.times_here(maxt))) 'ms'],'VerticalAlignment','bottom');
    xlabel('Frequencies (Hz)','fontsize',10); axis tight
    
    % time * frequency
    axes(handles.sub_display2); cla;
    imagesc(squeeze(handles.scale(maxe,:,:)));
    colormap(gca, handles.cc); img_prop = get(gca);
    newticks = round(linspace(1,length(handles.times_here),length(img_prop.YTick)));
    Xlabels = round(handles.times_here(newxticks)); set(gca,'XTickLabel', split(string(Xlabels)))
    newticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
    Ylabels = fliplr(round(handles.freqs_here(newticks))); set(gca,'YTickLabel', split(string(Ylabels)))
    title(['Frequency x time @ channel ' num2str(handles.LIMO.data.chanlocs(maxe).labels)]);
    xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequencies','fontsize',10);
end

handles.slider_sel=slider_sel;
guidata(hObject, plot_data);
guidata(hObject, handles);


% -------------------------------------------------
%                     MOUSE INPUT 
% -------------------------------------------------

% --- Executes on button press in mouse_input.
function mouse_input_Callback(hObject, eventdata, handles)

plot_data       = guidata(hObject);
popup_sel_index = get(handles.pop_up_dimensions, 'Value');
[x,y,button]    = ginput(1);

while button == 1
    clickedAx = gca;
    if clickedAx ==handles.Main_display
        if x < 1; x=1; end
        if y < 1; y=1; end
        popup_sel_index_str = {'Frequency', 'Time', 'Electrode'};
        timep = numel(handles.times_here);
        freqp = numel(handles.freqs_here);
    else
        break
    end
    
    if popup_sel_index == 1  % if showing elec x freq
        
        if x > numel(plot_data.freqs_here)
            x=numel(plot_data.freqs_here); 
        end
        freq = floor(x); channel = floor(y);
        
        
        % stat value
        axes(handles.tf_course_plot);
        if handles.slider_sel > handles.maxt
            handles.slider_sel = handles.maxt;
        end
        d = squeeze(handles.data3d(channel,:,handles.slider_sel)); axis tight
        plot(handles.freqs_here,d,'LineWidth',3); M = max(d(:));
        ylabel('Stat value','fontsize',10); xlabel('Frequencies','fontsize',10,'VerticalAlignment','top');
        title(sprintf('Max stat %g @ %gHz',M,handles.freqs_here(freq)),'VerticalAlignment','bottom'); grid on; axis tight
        
        % topoplot at max freq
        axes(handles.topoplot); cla
        topoplot(squeeze(handles.data3d(:,freq,handles.slider_sel)),handles.LIMO.data.chanlocs);
        colormap(gca, handles.cc(2:end,:));
        mytitle = sprintf('topoplot @%gms & %gHz', handles.times_here(handles.slider_sel), round(handles.freqs_here(freq)));
        title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')
        axes(handles.Main_display);
        colormap(gca, handles.cc);

        % channel * time
        axes(handles.sub_display1); cla;
        imagesc(squeeze(handles.scale(:,freq,:)));
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
        title(['Channels x Times @ ' num2str(round(handles.freqs_here(freq))) ' Hz'],'VerticalAlignment','bottom');
        xlabel('Time (ms)','fontsize',10); axis tight
        
        % time * frequency
        axes(handles.sub_display2); cla;
        imagesc(squeeze(handles.scale(channel,:,:)));
        colormap(gca, handles.cc); img_prop = get(gca);
        set(gca,'XTick',newxticks); set(gca,'XTickLabel', split(string(Xlabels)))
        newticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
        Ylabels = fliplr(round(handles.freqs_here(newticks))); set(gca,'YTickLabel', split(string(Ylabels)))
        title(['Frequency x time @ channel ' num2str(handles.LIMO.data.chanlocs(channel).labels)]);
        xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequencies','fontsize',10);
        
        try
            p_values = evalin('base','p_values');
            if ~isnan(p_values(channel,freq,handles.slider_sel))
                fprintf('Stat value: %g, p_value %g \n',handles.data3d(channel,freq,handles.slider_sel),p_values(channel,freq,handles.slider_sel));
            end
        catch pvalerror
            fprintf('couldn''t figure the p value?? %s \n',pvalerror.message)
        end
        
    elseif popup_sel_index == 2  % if showing elec x times on main
        
        if x > numel(plot_data.times_here); x=numel(plot_data.times_here); end
        time = floor(x); channel = floor(y);
        
        % stat value
        axes(handles.tf_course_plot);
        if handles.slider_sel > handles.maxf
            handles.slider_sel = handles.maxf;
        end
        d = squeeze(handles.data3d(channel,handles.slider_sel,:)); axis tight
        plot(handles.times_here,d,'LineWidth',3); M = max(d(:));
        ylabel('Stat value','fontsize',10); xlabel('Time','fontsize',10,'VerticalAlignment','top');
        title(sprintf('Max stat %g @ %gms',M,handles.times_here(handles.slider_sel)),'VerticalAlignment','bottom'); grid on; axis tight
        
        % topoplot at max freq
        axes(handles.topoplot); cla
        topoplot(squeeze(handles.data3d(:,handles.slider_sel,time)),handles.LIMO.data.chanlocs);
        colormap(gca, handles.cc(2:end,:));
        mytitle = sprintf('topoplot @%gms & %gHz', round(handles.times_here(time)),round(handles.freqs_here(handles.slider_sel)));
        title(mytitle, 'Units', 'normalized', 'Position', [1.2, 0.5],'Rotation',-90,'FontWeight','bold','VerticalAlignment','top')
        
        axes(handles.Main_display);
        colormap(gca, handles.cc);
        
        % channel * freq
        axes(handles.sub_display1); cla;
        imagesc(squeeze(handles.scale(:,:,time)));
        colormap(gca, handles.cc); img_prop = get(gca);
        newxticks = round(linspace(1,length(handles.freqs_here),length(img_prop.XTick)));
        Xlabels = round(handles.freqs_here(newxticks)); set(gca,'XTickLabel', split(string(Xlabels)))
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
        title(['Channels x Frequencies @ ' num2str(round(handles.times_here(time))) 'ms'],'VerticalAlignment','bottom');
        xlabel('Frequencies (Hz)','fontsize',10); axis tight
        
        % time * frequency
        axes(handles.sub_display2); cla;
        imagesc(squeeze(handles.scale(channel,:,:)));
        colormap(gca, handles.cc); img_prop = get(gca);
        newticks = round(linspace(1,length(handles.times_here),length(img_prop.YTick)));
        Xlabels = round(handles.times_here(newxticks)); set(gca,'XTickLabel', split(string(Xlabels)))
        newticks = round(linspace(1,length(handles.freqs_here),length(img_prop.YTick)));
        Ylabels = fliplr(round(handles.freqs_here(newticks))); set(gca,'YTickLabel', split(string(Ylabels)))
        title(['Frequency x time @ channel ' num2str(handles.LIMO.data.chanlocs(channel).labels)]);
        xlabel('Time (ms)','fontsize',10,'VerticalAlignment','top'); ylabel('Frequencies','fontsize',10);
        
        try
            p_values = evalin('base','p_values');
            if ~isnan(p_values(channel,handles.slider_sel,time))
                fprintf('Stat value: %g, p_value %g \n',handles.data3d(channel,handles.slider_sel,time),p_values(channel,handles.slider_sel,time));
            end
        catch pvalerror
            fprintf('couldn''t figure the p value?? %s \n',pvalerror.message)
        end
    end
    
    [x,y,button]=ginput(1);
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%% MAKE A MOVIE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in gifbutton.
function gifbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gifbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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

%% set axes and labels 
% -------------------------------------------------------------------------
function set_imgaxes(LIMO,scale)

img_prop = get(gca);
set(gca,'LineWidth',2)

% ----- X --------
if strcmp(LIMO.Analysis,'Time')
    xlabel('Time in ms','FontSize',10)
elseif strcmp(LIMO.Analysis,'Frequency')
    xlabel('Frequency in Hz','FontSize',10)
end
Xlabels = LIMO.data.start:LIMO.data.end;
newticks = round(linspace(1,length(Xlabels),length(img_prop.XTick)));
Xlabels = Xlabels(newticks);
set(gca,'XTick',newticks);
set(gca,'XTickLabel', split(string(Xlabels)))
 
% ----- Y --------
if strcmp(LIMO.Type,'Components')
    if size(scale,1) == 1
        ylabel('Optimized component','FontSize',10);
    else
        ylabel('Components','FontSize',10);
    end
else
    if size(scale,1) == 1
        ylabel('Optimized channel','FontSize',10);
    else
        ylabel('Channels','FontSize',10);
    end
end

if isfield(LIMO.data, 'chanlocs')
    Ylabels = arrayfun(@(x)(x.labels), LIMO.data.chanlocs, 'UniformOutput', false);
else
    Ylabels = arrayfun(@(x)(x.labels), LIMO.data.expected_chanlocs, 'UniformOutput', false);
end
newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
Ylabels = Ylabels(newticks);
set(gca,'YTick',newticks);
set(gca,'YTickLabel', Ylabels);

% ----- Colormap limits --------
try
    maxval = max(abs(max(scale(:))),abs(min(scale(:))));
    if max(scale(:)) < 0
        caxis([-maxval 0])
    elseif min(scale(:)) > 0 
        caxis([0 maxval])
    else
        caxis([-maxval maxval])
    end
catch caxiserror
    fprintf('axis issue: %s\n',caxiserror.message)
end


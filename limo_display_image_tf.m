function varargout = limo_display_image_tf(varargin)

% limo_display_results_tf: interactive GUI to display time*freq results
%
%      INPUT limo_display_results_tf(LIMO,toplot,mask,title)
%       1 - LIMO struct
%       2 - 3D matrix of values to plot, dim elec x freqs x time-bins
%       3 - 3D matrix of significant cells
%       4 - title (from limo_stat_values_tf)
%
% A. Stewart v1 mar14 axs - basic GUI set up + 3D tf plots from 3 point-of-view 
% C. Pernet added mask,title, - fixed various bugs - made it show topoplot and time courses +
% fixed the code to update at each 'click' 9-04-2014
% ---------------------------------------------------------------
%  Copyright (C) LIMO Team 2014

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
scale              = handles.data3d.*handles.mask; 
scale(scale==0)    = NaN;
handles.scale      = scale;
handles.title      = varargin{4};
handles.freqs_here = linspace(handles.LIMO.data.lowf,handles.LIMO.data.hightf,size(handles.data3d,2));
handles.times_here = linspace(handles.LIMO.data.start,handles.LIMO.data.end,size(handles.data3d,3));
handles.plot_sel   = 1;
clear varargin scale

% Find max values, save idx and value
if (size(handles.data3d,4)) == 3
    dim = inputdlg('which dimension to use','4D input - plotting option');
    handles.data3d = handles.data3d(:,:,:,dim);
end
handles.maxv         = max(handles.data3d(:));
handles.maxvi        = find(handles.data3d == handles.maxv);
if length(handles.maxvi) ~= 1
    handles.maxvi = handles.maxvi(1);
end
[handles.maxe, handles.maxf, handles.maxt] = ind2sub(size(handles.data3d), handles.maxvi);
handles.clims        = [0 handles.maxv];  % Set the global default scale of the colour bar to be 0:max value
handles.slider_sel   = handles.maxt; % 0.5
handles.ef           = squeeze(handles.data3d(:,:,handles.maxt));
plot_data.freqs_here = handles.LIMO.data.tf_freqs; % (handles.LIMO.data.trim_low_f:handles.LIMO.data.trim_high_f);
plot_data.times_here = handles.LIMO.data.tf_times; % (handles.LIMO.data.trim1:handles.LIMO.data.trim2);
guidata(hObject, plot_data);
guidata(hObject, handles);


% ----------------------------------------------------------
% This sets up the initial plot - only do when we are invisible
% By Default we plot Electrode * Frequencies (slide in time)
% so window can get raised using limo_display_results_tf.
% ----------------------------------------------------------

if strcmp(get(hObject,'Visible'),'off')
    ef = handles.ef;
    freqp = numel(handles.freqs_here);
    axes(handles.Main_display);
    imagesc(squeeze(handles.scale(:,:,handles.maxt)),handles.clims);
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    mytitle = sprintf('%s plotted at %g ms',handles.title, round(handles.times_here(handles.maxt)));
    title(mytitle,'fontsize',12);
    set(gca, 'XTick',[1 2 3 4 5 6],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
    xlabel('Frequency bin (Hz)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    try
        set(gca,'YTick',1:length(handles.LIMO.data.expected_chanlocs));
    catch ME
        set(gca,'YTick',1:length(handles.LIMO.data.chanlocs));
    end
    
    if handles.LIMO.Level == 1
        for i = 1 : length(handles.LIMO.data.chanlocs)
            try
                label_electrodes{i} = handles.LIMO.data.expected_chanlocs(i).labels;
            catch ME
                label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
            end
        end
    else
        if isempty(handles.LIMO.design.electrode)
            for i = 1 : length(handles.LIMO.data.chanlocs)
                label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
            end
        else
            if length(handles.LIMO.design.electrode) == 1
                label_electrodes = handles.LIMO.design.electrode;
            else
                label_electrodes = ' ';
                ylabel('optimized electrode','FontSize',14);
            end
        end
    end
    set(gca,'YTickLabel', label_electrodes);
    
    % topoplot
    axes(handles.topoplot);
    topoplot(ef(:,handles.maxf),handles.LIMO.data.chanlocs);
    title(['topoplot @ ' num2str(round(handles.freqs_here(handles.maxf))) ' Hz'],'FontSize',12)
    
    % time/freq courses
    axes(handles.tf_course_plot);
    plot(handles.freqs_here,ef(handles.maxe,:),'LineWidth',3); grid on; axis tight
    mytitle = sprintf('power spectrum @ \n electrode %s (%g)', handles.LIMO.data.chanlocs(handles.maxe).labels,handles.LIMO.data.chanlocs(handles.maxe).urchan);
    title(mytitle,'FontSize',12);
    
    % sub-displays
    axes(handles.sub_display1); cla;
    timep = numel(handles.times_here);
    imagesc(squeeze(handles.scale(:,handles.maxf,:)));
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    title(['Electrode x times @ ' num2str(round(handles.freqs_here(handles.maxf))) ' Hz'],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    
    axes(handles.sub_display2); cla;
    imagesc(flipud(squeeze(handles.scale(handles.maxe,:,:))));
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    set(gca, 'YTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    set(gca, 'YTickLabel',{round(handles.freqs_here(freqp)),round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(1))},'fontsize',12);
    title(['Frequency x time @ electrode ' num2str(handles.LIMO.data.chanlocs(handles.maxe).labels)],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Frequencies','fontsize',12);
    
end

function varargout = limo_display_results_tf_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% UIWAIT makes limo_display_results_tf wait for user response (see UIRESUME)
% uiwait(handles.figure1);

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
        
        % show electrode * freq
        ef = handles.data3d(:,:,handles.maxt);
        freqp = numel(handles.freqs_here);
        axes(handles.Main_display);
        imagesc(handles.scale(:,:,handles.maxt),handles.clims);
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        mytitle = sprintf('%s plotted at %g ms',handles.title, round(handles.times_here(handles.maxt)));
        title(mytitle,'fontsize',12);
        set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
        xlabel('Frequency bin (Hz)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
        try
            set(gca,'YTick',1:length(handles.LIMO.data.expected_chanlocs));
        catch ME
            set(gca,'YTick',1:length(handles.LIMO.data.chanlocs));
        end
        if handles.LIMO.Level == 1
            for i = 1 : length(handles.LIMO.data.chanlocs)
                try
                    label_electrodes{i} = handles.LIMO.data.expected_chanlocs(i).labels;
                catch ME
                    label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
                end
            end
        else
            if isempty(handles.LIMO.design.electrode)
                for i = 1 : length(handles.LIMO.data.chanlocs)
                    label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
                end
            else
                if length(handles.LIMO.design.electrode) == 1
                    label_electrodes = handles.LIMO.design.electrode;
                else
                    label_electrodes = ' ';
                    ylabel('optimized electrode','FontSize',14);
                end
            end
        end
        set(gca,'YTickLabel', label_electrodes);
        
        % topoplot
        axes(handles.topoplot); cla;
        topoplot(ef(:,handles.maxf),handles.LIMO.data.chanlocs);
        title(['topoplot @ ' num2str(round(handles.freqs_here(handles.maxf))) ' Hz'],'FontSize',12)
        
        % time/freq courses
        axes(handles.tf_course_plot); cla;
        plot(handles.freqs_here,ef(handles.maxe,:),'LineWidth',3); grid on; axis tight
        try
            mytitle = sprintf('power spectrum @ \n electrode %s (%g)', handles.LIMO.data.chanlocs(handles.maxe).labels,handles.LIMO.data.chanlocs(maxe).urchan);
        catch urchan_issue
            mytitle = sprintf('power spectrum @ \n electrode %s', handles.LIMO.data.chanlocs(handles.maxe).labels);
        end
        title(mytitle,'FontSize',12);
        
        % sub-displays
        axes(handles.sub_display1); cla;
        timep = numel(handles.times_here);
        imagesc(squeeze(handles.scale(:,handles.maxf,:)));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        title(['Electrode x times @ ' num2str(round(handles.freqs_here(handles.maxf))) ' Hz'],'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
        
        axes(handles.sub_display2); cla;
        imagesc(flipud(squeeze(handles.scale(handles.maxe,:,:))));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        set(gca, 'YTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
        set(gca, 'YTickLabel',{round(handles.freqs_here(freqp)),round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(1))},'fontsize',12);
        title(['Frequency x time @ electrode ' num2str(handles.LIMO.data.chanlocs(handles.maxe).labels)],'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Frequencies','fontsize',12);
        
        handles.ef = ef;
        
    case 2
        
        % show electrode * time
        et = squeeze(handles.data3d(:,handles.maxf,:));
        axes(handles.Main_display); cla
        imagesc(squeeze(handles.scale(:,handles.maxf,:)),handles.clims);
        mytitle = sprintf('%s plotted at %g Hz',handles.title, round(handles.freqs_here(handles.maxf)));
        title(mytitle,'fontsize',12);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
        try
            set(gca,'YTick',1:length(handles.LIMO.data.expected_chanlocs));
        catch ME
            set(gca,'YTick',1:length(handles.LIMO.data.chanlocs));
        end
        if handles.LIMO.Level == 1
            for i = 1 : length(handles.LIMO.data.chanlocs)
                try
                    label_electrodes{i} = handles.LIMO.data.expected_chanlocs(i).labels;
                catch ME
                    label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
                end
            end
        else
            if isempty(handles.LIMO.design.electrode)
                for i = 1 : length(handles.LIMO.data.chanlocs)
                    label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
                end
            else
                if length(handles.LIMO.design.electrode) == 1
                    label_electrodes = handles.LIMO.design.electrode;
                else
                    label_electrodes = ' ';
                    ylabel('optimized electrode','FontSize',14);
                end
            end
        end
        set(gca,'YTickLabel', label_electrodes);
        
        % topoplot
        axes(handles.topoplot); cla;
        topoplot(et(:,handles.maxt),handles.LIMO.data.chanlocs);
        title(['topoplot @ ' num2str(round(handles.times_here(handles.maxt))) ' ms'],'FontSize',12)
        
        % time/freq courses
        axes(handles.tf_course_plot); cla;
        plot(handles.times_here,et(handles.maxe,:),'LineWidth',3); grid on; axis tight
        try
            mytitle = sprintf('Amplitude @ \n electrode %s (%g)', handles.LIMO.data.chanlocs(handles.maxe).labels,handles.LIMO.data.chanlocs(maxe).urchan);
        catch urchan_issue
            mytitle = sprintf('Amplitude @ \n electrode %s', handles.LIMO.data.chanlocs(handles.maxe).labels);
        end
        title(mytitle,'FontSize',12);
        
        % sub-displays
        axes(handles.sub_display1); cla;
        imagesc(squeeze(handles.scale(:,:,handles.maxt)));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        title(['Electrode x Frequencies @' num2str(round(handles.times_here(handles.maxt))) ' ms'],'fontsize',12);
        set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
        xlabel('Frequencies (Hz)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
        
        axes(handles.sub_display2); cla;
        imagesc(flipud(squeeze(handles.scale(handles.maxe,:,:))));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        title(['Frequency x times @ electrode ' num2str(handles.LIMO.data.chanlocs(handles.maxe).labels)],'fontsize',12);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Frequencies','fontsize',12);
        
        handles.et = et;
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
plot_data = guidata(hObject);
popup_sel_index=get(handles.pop_up_dimensions, 'Value');
timep = numel(handles.times_here);
freqp = numel(handles.freqs_here);

if popup_sel_index==1;
    slider_sel = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.times_here)*slider_sel)); % Scale to get ints 1:5 out
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    end
    
    % show electrode * freq
    ef = handles.data3d(:,:,slider_sel);
    axes(handles.Main_display);
    imagesc(squeeze(handles.scale(:,:,slider_sel)),handles.clims);
    mytitle = sprintf('%s plotted at %g ms',handles.title, round(handles.times_here(slider_sel)));
    title(mytitle,'fontsize',12);
    set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    try
        set(gca,'YTick',1:length(handles.LIMO.data.expected_chanlocs));
    catch ME
        set(gca,'YTick',1:length(handles.LIMO.data.chanlocs));
    end
    if handles.LIMO.Level == 1
        for i = 1 : length(handles.LIMO.data.chanlocs)
            try
                label_electrodes{i} = handles.LIMO.data.expected_chanlocs(i).labels;
            catch ME
                label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
            end
        end
    else
        if isempty(handles.LIMO.design.electrode)
            for i = 1 : length(handles.LIMO.data.chanlocs)
                label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
            end
        else
            if length(handles.LIMO.design.electrode) == 1
                label_electrodes = handles.LIMO.design.electrode;
            else
                label_electrodes = ' ';
                ylabel('optimized electrode','FontSize',14);
            end
        end
    end
    set(gca,'YTickLabel', label_electrodes);
    [maxe,maxf]=ind2sub(size(ef),find(ef == max(ef(:))));
    
    % topoplot at max freq
    axes(handles.topoplot); cla
    topoplot(ef(:,maxf),handles.LIMO.data.chanlocs);
    title(['topoplot @ ' num2str(round(handles.freqs_here(maxf))) ' Hz'],'FontSize',12)
    
    % time/freq courses at max elec
    axes(handles.tf_course_plot); cla
    plot(handles.freqs_here,ef(maxe,:),'LineWidth',3); grid on; axis tight
    try
        mytitle = sprintf('power spectrum @ \n electrode %s (%g)', handles.LIMO.data.chanlocs(maxe).labels,handles.LIMO.data.chanlocs(maxe).urchan);
    catch urchan_issue
        mytitle = sprintf('power spectrum @ \n electrode %s', handles.LIMO.data.chanlocs(maxe).labels);
    end
    title(mytitle,'FontSize',12);
    
    % sub-display
    axes(handles.sub_display1); cla;
    timep = numel(handles.times_here);
    imagesc(squeeze(handles.scale(:,maxf,:)));
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    title(['Electrode x times @ ' num2str(round(handles.freqs_here(maxf))) ' Hz'],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    
    axes(handles.sub_display2); cla;
    imagesc(flipud(squeeze(handles.scale(maxe,:,:))));
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    set(gca, 'YTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    set(gca, 'YTickLabel',{round(handles.freqs_here(freqp)),round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(1))},'fontsize',12);
    title(['Frequency x time @ electrode ' num2str(handles.LIMO.data.chanlocs(maxe).labels)],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Frequencies','fontsize',12);
    
    handles.ef = ef;
end

if popup_sel_index==2;
    slider_sel2 = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.freqs_here)*slider_sel2)); % Scale slider to correct ints
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    end
    
    % show electrode * time
    et = squeeze(handles.data3d(:,slider_sel,:));
    axes(handles.Main_display); cla;
    imagesc(squeeze(handles.scale(:,slider_sel,:)),handles.clims);
    title(['Values at freq of ',num2str(round(handles.freqs_here(slider_sel))),' Hz'],'fontsize',12);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    try
        set(gca,'YTick',1:length(handles.LIMO.data.expected_chanlocs));
    catch ME
        set(gca,'YTick',1:length(handles.LIMO.data.chanlocs));
    end
    if handles.LIMO.Level == 1
        for i = 1 : length(handles.LIMO.data.chanlocs)
            try
                label_electrodes{i} = handles.LIMO.data.expected_chanlocs(i).labels;
            catch ME
                label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
            end
        end
    else
        if isempty(handles.LIMO.design.electrode)
            for i = 1 : length(handles.LIMO.data.chanlocs)
                label_electrodes{i} = handles.LIMO.data.chanlocs(i).labels;
            end
        else
            if length(handles.LIMO.design.electrode) == 1
                label_electrodes = handles.LIMO.design.electrode;
            else
                label_electrodes = ' ';
                ylabel('optimized electrode','FontSize',14);
            end
        end
    end
    set(gca,'YTickLabel', label_electrodes);
    [maxe,maxt]=ind2sub(size(et),find(et == max(et(:))));
    
    % topoplot at max time
    axes(handles.topoplot); cla;
    topoplot(et(:,maxt),handles.LIMO.data.chanlocs);
    title(['topoplot @ ' num2str(round(handles.times_here(maxt))) ' ms'],'FontSize',12)
    
    % time/freq courses at max elec
    axes(handles.tf_course_plot); cla;
    plot(handles.times_here,et(maxe,:),'LineWidth',3); grid on; axis tight
    mytitle = sprintf('Amplitude @ \n electrode %s (%g)', handles.LIMO.data.chanlocs(maxe).labels,handles.LIMO.data.chanlocs(maxe).urchan);
    title(mytitle,'FontSize',12);
    
    % sub-displays
    axes(handles.sub_display1); cla;
    imagesc(squeeze(handles.scale(:,:,maxt)));
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    title(['Electrode x Frequencies @' num2str(round(handles.times_here(maxt))) ' ms'],'fontsize',12);
    set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
    xlabel('Frequencies (Hz)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    
    axes(handles.sub_display2); cla;
    imagesc(flipud(squeeze(handles.scale(maxe,:,:))));
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    title(['Frequency x times @ electrode ' num2str(handles.LIMO.data.chanlocs(maxe).labels)],'fontsize',12);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Frequencies','fontsize',12);
    
    handles.et = et;
end

handles.slider_sel=slider_sel;
guidata(hObject, handles);


% -------------------------------------------------
%                     MOUSE INPUT 
% -------------------------------------------------

% --- Executes on button press in mouse_input.
function mouse_input_Callback(hObject, eventdata, handles)

plot_data = guidata(hObject);
popup_sel_index=get(handles.pop_up_dimensions, 'Value');
[x,y,button]=ginput(1);
while button == 1
    fprintf('The current mouse location is: %g %g \n',x,y);
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
        
        if x > numel(plot_data.freqs_here); x=numel(plot_data.freqs_here); end
        sel_str = ['Selection is at electrode ',num2str(floor(y)),' (',plot_data.LIMO.data.chanlocs(1,floor(y)).labels,') and at a freq of ',num2str(plot_data.freqs_here(floor(x))),' Hz.'];
        x = floor(x); y = floor(y);
               
        % topoplot at selected freq
        axes(handles.topoplot); cla;
        topoplot(handles.ef(:,x),handles.LIMO.data.chanlocs);
        title(['topoplot @ ' num2str(round(handles.freqs_here(x))) ' Hz'],'FontSize',12)
        
        % time/freq courses at selected elec
        axes(handles.tf_course_plot); 
        plot(handles.freqs_here,handles.ef(y,:),'LineWidth',3); grid on; axis tight
        try
            mytitle = sprintf('power spectrum @ \n electrode %s (%g)', handles.LIMO.data.chanlocs(y).labels,handles.LIMO.data.chanlocs(maxe).urchan);
        catch urchan_issue
            mytitle = sprintf('power spectrum @ \n electrode %s', handles.LIMO.data.chanlocs(y).labels);
        end
        title(mytitle,'FontSize',12);
        
        
        % sub-display
        axes(handles.sub_display1); cla;
        timep = numel(handles.times_here);
        imagesc(squeeze(handles.scale(:,x,:)));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        title(['Electrode x times @ ' num2str(round(handles.freqs_here(x))) ' Hz'],'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
        
        axes(handles.sub_display2); cla;
        imagesc(flipud(squeeze(handles.scale(y,:,:))));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        set(gca, 'YTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
        set(gca, 'YTickLabel',{round(handles.freqs_here(freqp)),round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(1))},'fontsize',12);
        title(['Frequency x time @ electrode ' num2str(handles.LIMO.data.chanlocs(y).labels)],'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Frequencies','fontsize',12);
        
        
    elseif popup_sel_index == 2  % if showing elec x times on main
        
        if x > numel(plot_data.times_here); x=numel(plot_data.times_here); end
        sel_str = ['Selection is at electrode ',num2str(floor(y)),' (',plot_data.LIMO.data.chanlocs(1,floor(y)).labels,') and at a time of ',num2str(plot_data.times_here(floor(x))),' ms.'];
        x = floor(x); y = floor(y);
        
        % topoplot at max time
        axes(handles.topoplot); cla;
        topoplot(handles.et(:,x),handles.LIMO.data.chanlocs);
        title(['topoplot @ ' num2str(round(handles.times_here(x))) ' ms'],'FontSize',12)
        
        % time/freq courses at max elec
        axes(handles.tf_course_plot); cla;
        plot(handles.times_here,handles.et(y,:),'LineWidth',3); grid on; axis tight
        try
            mytitle = sprintf('Amplitude @ \n electrode %s (%g)', handles.LIMO.data.chanlocs(y).labels,handles.LIMO.data.chanlocs(maxe).urchan);
        catch urchan_issue
            mytitle = sprintf('Amplitude @ \n electrode %s', handles.LIMO.data.chanlocs(y).labels);
        end
        title(mytitle,'FontSize',12);
        
        % sub-displays
        axes(handles.sub_display1); cla;
        imagesc(squeeze(handles.scale(:,:,x)));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        title(['Electrode x Frequencies @' num2str(round(handles.times_here(x))) ' ms'],'fontsize',12);
        set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
        xlabel('Frequencies (Hz)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
        
        axes(handles.sub_display2); cla;
        imagesc(flipud(squeeze(handles.scale(y,:,:))));
        cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
        title(['Frequency x times @ electrode ' num2str(handles.LIMO.data.chanlocs(y).labels)],'fontsize',12);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Frequencies','fontsize',12);
        
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




%%%%%%%%%%%%%%%%%%%% MENU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)

[file,p] = uigetfile('*.mat','select effect file');
if ~isequal(file, 0)
    cd(p); load(file); load LIMO;
    do_not_update = 0;
    if strncmp(file,'Covariate_effect',16)
        data3d     = squeeze(Covariate_effect(:,:,:,1));
    elseif strncmp(file,'Condition_effect',16)
        data3d     = squeeze(Condition_effect(:,:,:,1));
    elseif strncmp(file,'R2',2)
        data3d     = squeeze(R2(:,:,:,2));
    elseif strncmp(file,'one_sample',10)
        data3d     = squeeze(one_sample(:,:,:,4));
    else
        do_not_update = 1;
        errordlg('file not supported')
    end
    
    % ---------------------
    if do_not_update == 0
        clc; uiresume
        guidata(hObject, handles);
        delete(handles.figure1)
        limo_display_results_tf(LIMO,data3d,ones(size(data3d)),file)
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

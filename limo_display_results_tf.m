function varargout = limo_display_results_tf(varargin)
% limo_display_results_tf MATLAB code for limo_display_results_tf.fig
%      
%      limo_display_results_tf shows an interactive plot to show 3D plot of
%      electrodes x frequency x time-point data, in a 3D matrix passed to
%      it.
%
%      INPUT
%       1 - LIMO struct
%       2 - 3D matrix of values to plot, orientated in elec x freqs x
%       time-bins
%

%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help limo_display_results_tf

% v1 mar14 axs - 3D tf plots from 3 point-of-view working

% Begin initialization code - DO NOT EDIT
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
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to limo_display_results_tf (see VARARGIN)

% Choose default command line output for limo_display_results_tf

handles.output = hObject;
handles.LIMO = varargin{1};
handles.data3d = varargin{2};

handles.freqs_here = handles.LIMO.data.tf_freqs(handles.LIMO.data.trim_low_f:handles.LIMO.data.trim_high_f);
handles.times_here = handles.LIMO.data.tf_times(handles.LIMO.data.trim_low_t:handles.LIMO.data.trim_high_t);

handles.plot_sel=1;

% Find max values, save idx and value
if (size(handles.data3d,4)) == 3
    handles.data3d = handles.data3d(:,:,:,1);
    disp('4D passed to display tf, taking first dim')
end
   




handles.maxv = max(handles.data3d(:));
handles.maxvi = find(handles.data3d == handles.maxv);

[handles.maxe, handles.maxf, handles.maxt] = ind2sub(size(handles.data3d), handles.maxvi);

handles.clims = [0 handles.maxv];  % Set the global default scale of the colour bar to be 0:max value


handles.slider_sel = 0.5;

plot_data.freqs_here = handles.LIMO.data.tf_freqs(handles.LIMO.data.trim_low_f:handles.LIMO.data.trim_high_f);
plot_data.times_here = handles.LIMO.data.tf_times(handles.LIMO.data.trim_low_t:handles.LIMO.data.trim_high_t);


guidata(hObject, plot_data);



% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using limo_display_results_tf.
if strcmp(get(hObject,'Visible'),'off')
    disp('Electrodes x Frequencies across times')
    
    
    ef = handles.data3d(:,:,handles.maxt);
    disp(handles.maxt)
    
    freqp = numel(handles.freqs_here);
    axes(handles.axes1);
    imagesc(ef,handles.clims);
    title(['Values at time of ',num2str(round(handles.times_here(handles.maxt))),' ms'],'fontsize',12);
    set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
    xlabel('Frequency bin (Hz)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    
    axes(handles.tplot);
    % Topoplot here 
    
    clear et ft
    axes(handles.axes2); cla;
    et = squeeze(handles.data3d(:,handles.maxf,:));
    imagesc(et);
    title(['Electrode x times'],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    
    axes(handles.axes3); cla;
    ft = squeeze(handles.data3d(handles.maxe,:,:));
    imagesc(ft);
    title(['Frequency x times'],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Frequencies','fontsize',12);
    
    
    
end

% UIWAIT makes limo_display_results_tf wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = limo_display_results_tf_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;


timep = numel(handles.times_here);
freqp = numel(handles.freqs_here);

popup_sel_index = get(handles.popupmenu1, 'Value');




switch popup_sel_index
    case 1
        disp('Electrodes x frequencies across times')
        handles.plot_sel=1;
        disp(handles.plot_sel)
        
        ef = handles.data3d(:,:,handles.maxt,1);
        imagesc(ef,handles.clims);
        title(['Values at time of ',num2str(round(handles.times_here(handles.maxt))),' ms'],'fontsize',12);
        set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
        
        
        xlabel('Frequency bin (Hz)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
        
        
    case 2
        disp('Electrodes x Times across frequencies')
        handles.plot_sel=2;
        disp(handles.plot_sel)
        
        
        
        et = squeeze(handles.data3d(:,handles.maxf,:));
        imagesc(et,handles.clims);
        title(['Values at freq of ',num2str(round(handles.freqs_here(handles.maxf))),' Hz'],'fontsize',12);
        set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
        set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
        xlabel('Time (ms)','fontsize',12);
        ylabel('Electrodes','fontsize',12);
end







% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
%     ['Close ' get(handles.figure1,'Name') '...'],...
%     'Yes','No','Yes');
% if strcmp(selection,'No')
%     return;
% end



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'Electrodes x frequencies across times', 'Electrodes x Times across frequencies'});


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_data = guidata(hObject);


%disp(handles.plot_sel)
popup_sel_index=get(handles.popupmenu1, 'Value');

timep = numel(handles.times_here);
freqp = numel(handles.freqs_here);


if popup_sel_index==1;
    
    slider_sel = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.times_here)*slider_sel)); % Scale to get ints 1:5 out
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    end
    
    %disp(slider_sel);
    
    
    
    
    
    
    ef = handles.data3d(:,:,slider_sel,1);
    
    axes(handles.axes1);
    imagesc(ef,handles.clims);
    title(['Values at time of ',num2str(round(handles.times_here(slider_sel))),' ms'],'fontsize',12);
    set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
    xlabel('Frequency (Hz)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    set(gca, 'YTickLabel',plot_data.LIMO.data.chanlocs(1,:));
    
    %tstr = ['Topoplot at ',num2str(round(handles.times_here(slider_sel))),' ms'];
    axes(handles.tplot);
    cla;
    % topoplot here
    
    
    
end

if popup_sel_index==2;
    slider_sel2 = get(hObject,'Value');
    slider_sel = int32(ceil(numel(handles.freqs_here)*slider_sel2)); % Scale slider to correct ints
    if slider_sel == 0
        slider_sel = 1; % Don't want the 0th entry, set to 1 instead
    end
    
    %disp(slider_sel);
    
    
    et = squeeze(handles.data3d(:,slider_sel,:,1));
    axes(handles.axes1); cla;
    imagesc(et,handles.clims);
    title(['Values at freq of ',num2str(round(handles.freqs_here(slider_sel))),' Hz'],'fontsize',12);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    
    xlabel('Time (ms)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
end

handles.slider_sel=slider_sel;


% Update handles structure
guidata(hObject, handles);



% set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep]);
% set(gca, 'XTickLabel',{handles.times_here(1), handles.times_here(round(timep/4)), handles.times_here(round(timep/2)), handles.times_here(round(3*timep/4)),handles.times_here(timep)});



% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in gifbutton.
function gifbutton_Callback(hObject, eventdata, handles)
% hObject    handle to gifbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in ginput.
function ginput_Callback(hObject, eventdata, handles)
% hObject    handle to ginput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plot_data = guidata(hObject);
popup_sel_index=get(handles.popupmenu1, 'Value');



[x,y,button]=ginput(1);
mousestr = ['The current mouse location is:',num2str(x),'x and ',num2str(y)];
disp(mousestr);
popup_sel_index_str = {'Frequency', 'Time', 'Electrode'};
% Ensure x and y are within range
if x < 1; x=1; end
if y < 1; y=1; end

timep = numel(handles.times_here);
freqp = numel(handles.freqs_here);


if popup_sel_index == 1  % if showing elec x freq
    
    
    
    if x > numel(plot_data.freqs_here); x=numel(plot_data.freqs_here); end
    %mousestr2 = ['Selection is at electrode ',LIMO.data.chanlocs(1,y).labels,' and at a freq of ',plot_data.freqs_here(floor(x)),' Hz.'];
    sel_str = ['Selection is at electrode ',num2str(floor(y)),' (',plot_data.LIMO.data.chanlocs(1,floor(y)).labels,') and at a freq of ',num2str(plot_data.freqs_here(floor(x))),' Hz.'];
    
    clear et ft
    axes(handles.axes2); cla;
    et = squeeze(handles.data3d(:,floor(x),:));
    imagesc(et);
    title(['Electrode x times'],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    
    
    axes(handles.axes3); cla;
    ft = squeeze(handles.data3d(floor(y),:,:));
    imagesc(flipud(ft));
    title(['Frequency x times, at electrode ',plot_data.LIMO.data.chanlocs(1,floor(y)).labels],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Frequencies','fontsize',12);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    ftyticks={round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))};
    ftyticks=fliplr(ftyticks);
    set(gca, 'YTickLabel',ftyticks,'fontsize',12);
    
    elec_here = plot_data.LIMO.data.chanlocs(1,floor(y)).labels;
    
    axes(handles.erpplot);
    
    
    
elseif popup_sel_index == 2  % if showing elec x times on main
    
    if x > numel(plot_data.times_here); x=numel(plot_data.times_here); end
    sel_str = ['Selection is at electrode ',num2str(floor(y)),' (',plot_data.LIMO.data.chanlocs(1,floor(y)).labels,') and at a time of ',num2str(plot_data.times_here(floor(x))),' ms.'];
    
    axes(handles.tplot); cla;
    %pop_topoplot(plot_data.LIMO.data.chanlocs,1, plot_data.times_here(floor(x)),'selected time',[1 1] ,0,'colorbar','off');
    topoplot(plot_data.times_here(floor(x)),plot_data.LIMO.data.chanlocs);
    
    axes(handles.axes2); cla;
    ef = squeeze(handles.data3d(:,:,floor(x)));
    imagesc(ef);
    title(['Electrodes x frequencies, at ',num2str(round(plot_data.times_here(floor(x)))),'ms'],'fontsize',12);
    xlabel('Freq (Hz)','fontsize',12);
    ylabel('Electrodes','fontsize',12);
    set(gca, 'XTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))},'fontsize',12);
    
    
    axes(handles.axes3); cla;
    ft = squeeze(handles.data3d(floor(y),:,:));
    
    imagesc(flipud(ft));
    
    title(['Frequency x times, at electrode ',plot_data.LIMO.data.chanlocs(1,floor(y)).labels],'fontsize',12);
    xlabel('Time (ms)','fontsize',12);
    ylabel('Frequencies','fontsize',12);
    set(gca, 'XTick',[1 timep/4 timep/2 3*timep/4 timep],'fontsize',12);
    set(gca, 'XTickLabel',{round(handles.times_here(1)), round(handles.times_here(round(timep/4))), round(handles.times_here(round(timep/2))), round(handles.times_here(round(3*timep/4))),round(handles.times_here(timep))},'fontsize',12);
    set(gca, 'YTick',[1 freqp/4 freqp/2 3*freqp/4 freqp],'fontsize',12);
    
    ftyticks={round(handles.freqs_here(1)), round(handles.freqs_here(round(freqp/4))), round(handles.freqs_here(round(freqp/2))), round(handles.freqs_here(round(3*freqp/4))),round(handles.freqs_here(freqp))};
    ftyticks=fliplr(ftyticks);
    set(gca, 'YTickLabel',ftyticks,'fontsize',12);
    
    
    
end

disp(sel_str)

set(handles.sel_text,'String',sel_str);








% --- Executes during object creation, after setting all properties.
function sel_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function tplot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate tplot


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

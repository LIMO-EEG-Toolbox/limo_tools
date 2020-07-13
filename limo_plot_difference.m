function Data = limo_plot_difference(varargin)

% allows to plot the Difference between two data set with alpha % Hightest
% Density Intervals (ie Bayesian bootstrap CI) 
% WARNING this is computed channel/time frame wise, ie this is  not simultaneous CI
%
% FORMATS
% Diff = limo_plot_Difference
% Diff = limo_plot_Difference(data1,data2,'type','paired/independent')
% Diff = limo_plot_Difference(data1,data2,'type','paired/independent','percent', 20, 'alpha', 0.05, 'fig', 'on')
%
% INPUTS
% data1/2 matrices of data up to 4D, if 4D because coming from robust averaging, 
%         i.e. the last dimension is the estimator with CI, only the estimator is 
%         used ie level 2 if dim 4
% type    'paired' or 'independent'
% percent the amount of trimming 0% is a mean, 20% is the default trimmed mean, 
%         50% is the median (in case 50 is used, the Harrell-Davis estimator
%         of the median is used)
% alpha   is the 1-alpha level of the HDI
% fig     'on' (default) or 'off' indicates to produce a figure or not
%
% OUPUT
% Data    a structure with
%         .Diff the 3D matrix of Difference with HDI
%         .limo the LIMO structure for channel/time-freq info (optional) 
%
% Cyril Pernet - all options in Septembre 2019
% ---------------------------------------------
%  Copyright (C) LIMO Team 2019

%% check inputs

percent     = 20/100; % defines the amount of trimming done
alpha_level = 5/100;  % 1-alpha CI
figure_flag = 1;      % make a figure
wrapdata    = 0;      % for 4D wrap into 3D

if nargin < 3 
    % select 1st dataset 
    % ------------------
    [file,locpath,idx]=uigetfile('.mat','Select 1st dataset');
    if idx == 0
        return
    end
    cd(locpath); 
    data1 = load(file); 
    data1 = data1.(cell2mat(fieldnames(data1)));
    if isfield(data1,'limo')
        limo1 = data1.limo;
    end
    if isfield(data1,'data')
        data1 = data1.data;
    end
    if isnumeric(data1)
        fprintf('%s loaded \n',file);
    else
       error('couldn''t load the data') 
    end
        
    % select 2nd dataset 
    % ------------------
    [file,locpath,idx]=uigetfile('.mat','Select 1st dataset');
    if idx == 0
        return
    end
    cd(locpath);
    data2 = load(file); 
    data2 = data2.(cell2mat(fieldnames(data2)));
    if isfield(data2,'limo')
        limo2 = data2.limo;
    end
    if isfield(data2,'data')
        data2 = data2.data;
    end
    if isnumeric(data2)
        fprintf('%s loaded \n',file);
    else
       error('couldn''t load the data') 
    end
    
    % type 
    % ----
    type = questdlg('are the data','analysis option','paired','independent','paired');
    if isempty(type)
        warning('selection aborded'); return
    end
    
    % percent 
    % ----
    percent = questdlg('which summary','analysis option','mean','20% trimmed mean','median','20% trimmed mean');
    if isempty(percent)
        warning('selection aborded'); return
    elseif strcmpi(percent,'mean')     
        percent = 0;
    elseif strcmpi(percent,'20% trimmed mean')     
        percent = 20;
    elseif strcmpi(percent,'median')     
        percent = 50;       
    end
    
    % alpha_level
    % -----------
    v = inputdlg('set alpha_level level (%) for confidence intervals');
    if isempty(v)
        return
    else
        alpha_level = eval(cell2mat(v));
        if alpha_level > 1
            alpha_level = alpha_level / 100;
        end
    end    
    
elseif nargin>= 4 && nargin <=6
    data1       = varargin{1};
    data2       = varargin{2};
    type        = varargin{3};
    percent     = varargin{4};
    if nargin >= 5
        alpha_level = varargin{5}; 
        if alpha_level > 1
            alpha_level = alpha_level / 100;
        end
    end
    if nargin == 6
        figure_flag = varargin{5};
    end
elseif nargin > 5
    error('too many arguments')
end

if alpha_level == 0
    disp('stange value for alpha - adjusted to 5%');
    alpha_level = 5/100;
end

%% check possible dimensions issues
squish = size(data1)==1;
if squish(end-1)==1
    if size(data1,1) == 1
        tmp = squeeze(data1); clear data1
        if ndims(tmp) == 3
            data1(1,:,:,:) = squeeze(tmp); clear tmp;
        else
            data1(1,:,:) = squeeze(tmp); clear tmp;
        end
    else
        data1 = squeeze(data1);
    end
end

squish = size(data2)==1;
if squish(end-1)==1
    if size(data2,1) == 1
        tmp = squeeze(data2); clear data2
        if ndims(tmp) == 3
            data2(1,:,:,:) = squeeze(tmp); clear tmp;
        else
            data2(1,:,:) = squeeze(tmp); clear tmp;
        end
    else
        data2 = squeeze(data2);
    end
end

if strcmpi(type,'paired')
    if numel(data1) ~= numel(data2)
        error('for paired data, each matrix must be of the same size')
    end
else
    if ndims(data1) ~= ndims(data2)
        error('the number of dimensions between datasets don''t match');
    end
end

% if 4D, wrap dim 2 and 3
if ndims(data1) == 4 && ndims(data2) == 4
    wrapdata      = 1;
    [C1,F1,T1,S1] = size(data1);
    data1         = limo_tf_4d_reshape(data1,[C1,F1*T1,S1]);
    [C2,F2,T2,S2] = size(data2);
    data2         = limo_tf_4d_reshape(data2,[C2,F2*T2,S2]);
end

if size(data1,3) < 6
    error('there are only %g observations to bootstrap, computation aborded (min=6)',size(data1,3))
end
if size(data2,3) < 6
    error('there are only %g observations to bootstrap, computation aborded (min=6)',size(data1,3))
end
Diff = NaN(size(data1,1),size(data1,2),3);

%% compute the mean Difference and CI

if strcmpi(type,'Paired')
    
        % the Difference 
        D = data1-data2;
        for channel=size(data1,1):-1:1
            fprintf('bootstraping data for CI estimation channel %g \n',channel)
            if percent == 0
                [est1(channel,:),CI1(channel,:,:)] = limo_central_estimator(squeeze(data1(channel,:,:)),'Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:)] = limo_central_estimator(squeeze(data2(channel,:,:)),'Mean',1-alpha_level);
                [Diff(channel,:,2),CID]            = limo_central_estimator(squeeze(D(channel,:,:)),'Mean',1-alpha_level);
                Diff(channel,:,1)                  = CID(1,:);
                Diff(channel,:,3)                  = CID(2,:);
            elseif percent == 20
                [est1(channel,:),CI1(channel,:,:)] = limo_central_estimator(squeeze(data1(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:)] = limo_central_estimator(squeeze(data2(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [Diff(channel,:,2),CID]            = limo_central_estimator(squeeze(D(channel,:,:)),'Trimmed Mean',1-alpha_level);
                Diff(channel,:,1)                  = CID(1,:); 
                Diff(channel,:,3)                  = CID(2,:);
            elseif percent == 50
                [est1(channel,:),CI1(channel,:,:)] = limo_central_estimator(squeeze(data1(channel,:,:)),'HD',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:)] = limo_central_estimator(squeeze(data2(channel,:,:)),'HD',1-alpha_level);
                [Diff(channel,:,2),CID]            = limo_central_estimator(squeeze(D(channel,:,:)),'HD',1-alpha_level);
                Diff(channel,:,1)                  = CID(1,:); 
                Diff(channel,:,3)                  = CID(2,:);
            end
        end
        
        % --------------------------------------------------------
        
elseif strcmpi(type,'Independent')
        
        est1 = NaN(size(data1,1),size(data1,2));   est2 = est1;
        CI1  = NaN(size(data1,1),2,size(data1,2)); CI2 = CI1;
        for channel=size(data1,1):-1:1
            fprintf('bootstraping data for CI estimation channel %g \n',channel)
            
            if percent == 0
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'Mean',1-alpha_level);
                Diff(channel,:,2)                      = est1(channel,:)-est2(channel,:);
            elseif percent == 20
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'Trimmed Mean',1-alpha_level);
                Diff(channel,:,2)                      = est1(channel,:)-est2(channel,:);
            elseif percent == 50
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'HD',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'HD',1-alpha_level);
                Diff(channel,:,2)                      = est1(channel,:)-est2(channel,:);
            end
            
            sorted_data   = sort(sort(bb1)-sort(bb2),2);
            upper_centile = floor((1-alpha_level)*size(sorted_data,2)); % upper bound
            nCIs          = size(sorted_data,2) - upper_centile;
            for frame = 1:size(sorted_data,1)
                tmp       = sorted_data(frame,:);
                ci        = 1:nCIs; 
                ciWidth   = tmp(ci+upper_centile) - tmp(ci); % all centile distances
                [~,index] = find(ciWidth == min(ciWidth));   % densest centile
                if length(index) > 1; index = index(1); end  % many similar values
                Diff(channel,frame,1) = tmp(index);
                Diff(channel,frame,3) = tmp(index+upper_centile);
            end
        end
end

% if 4D, unwrap dim 2 and 3
if wrapdata == 1
    est1 = reshape(est1,[C1,F1,T1]);
    CI1  = reshape(CI1,[C1,2,F1,T1]);
    est2 = reshape(est2,[C2,F2,T2]);
    CI2  = reshape(CI2,[C1,2,F1,T1]);
    Diff = limo_tf_4d_reshape(Diff, [C1,F1,T1,3]);
end

%% time/freq info
if exist('limo1','var'); LIMO = limo1; end
if exist('limo2','var'); LIMO = limo2; end
if ~exist('LIMO','var')
    [file,locpath,ind]=uigetfile({'LIMO.mat'},'Select any LIMO with right info');
    if ind == 0
        return
    else
        if strcmpi(file,'LIMO.mat')
            LIMO = load(fullfile(locpath,file));
            LIMO = LIMO.LIMO;
        else
            warning('selection aborded'); return
        end
    end
end

%% save
name      = cell2mat(inputdlg('save as [?]','name option'));
newname   = sprintf('%s',name);
newpath   = uigetdir('destination folder?');
Data.Diff = Diff;
if exist('LIMO','var')
    Data.limo = LIMO;
end
save (fullfile(newpath,newname),'Data');

%% Plot
% -----
if strcmp(figure_flag,'on') || figure_flag == 1
    
    if ndims(Diff) == 4
        whichdim = questdlg('which domain to plot?','averaging dimension','time','frequency','time');
        if isempty(whichdim)
            return
        elseif strcmpi(whichdim,'Frequency')
            est1 = squeeze(mean(est1,3));
            CI1  = squeeze(mean(CI1 ,4));
            est2 = squeeze(mean(est2,3));
            CI2  = squeeze(mean(CI2 ,4));
            Diff = squeeze(mean(Diff,3));
        elseif strcmpi(whichdim,'Time')
            est1 = squeeze(mean(est1,2));
            CI1  = squeeze(mean(CI1 ,3));
            est2 = squeeze(mean(est2,2));
            CI2  = squeeze(mean(CI2 ,3));
            Diff = squeeze(mean(Diff,2));
        else
            return
        end
    end
    
    % channel info
    if size(Diff,1) > 1
        channel = inputdlg('which channel to plot','Plotting option');
        if isempty(channel)  % cancel
            return
        elseif strcmp(channel,'') % ok empty
            [~,channel]=max(max(squeeze(Diff(:,:,2)),[],2));
        else
            channel = eval(cell2mat(channel));
            if length(channel) > 1
                error('1 channel only can be plotted')
            elseif channel > size(Diff,1)
                error('channel number invalid')
            end
        end
    else
        channel = 1;
    end
    
    if strcmpi(LIMO.Analysis,'Time')
        if isfield(data.LIMO.data,'timevect')
            vect = LIMO.data.timevect;
        else
            vect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;  % in msec
        end
    elseif strcmpi(LIMO.Analysis,'Frequency')
        if isfield(LIMO.data,'freqlist')
            vect = LIMO.data.freqlist;
        else
            vect = linspace(LIMO.data.start,LIMO.data.end,size(Diff,2));
        end
    elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        if strcmpi(whichdim,'Time')
            if isfield(LIMO.data,'tf_times')
                vect = LIMO.data.tf_times;
            else
                vect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;  % in msec
            end
        elseif strcmpi(whichdim,'Frequency')
            if isfield(LIMO.data,'tf_freqs')
                vect = LIMO.data.tf_freqs;
            else
                vect = linspace(LIMO.data.low_f,LIMO.data.high_f,size(Diff,2));
            end
        end
    end
    
    if  length(vect) ~= length(squeeze(est1(channel,:)))
        if length(vect) ~= length(squeeze(est1(channel,:)))
                fprintf('error in computing %s frames \n',LIMO.Analysis)
        end
        dlg = sprintf('enter %s interval by hand e.g. [0:0.5:40]',LIMO.Analysis);
        v = inputdlg(dlg);
        if isempty(v)
            return
        else
            try
                vect = eval(cell2mat(v));
                if length(vect) ~= size(Diff,2)
                    fprintf('%s interval invalid format \n',LIMO.Analysis)
                    vect = 1:size(Diff,2);
                end
            catch ME
                fprintf('%s interval invalid format \n',ME.message)
                vect = 1:size(Diff,2);
            end
        end
    end
    
    % figure
    figure;set(gcf,'Color','w'); subplot(3,2,[1 2 3 4]); hold on
    plot(vect,squeeze(est1(channel,:)),'LineWidth',3);
    fillhandle = patch([vect fliplr(vect)],[squeeze(CI1(channel,1,:))' fliplr(squeeze(CI1(channel,2,:))')], [0 0 1]);
    set(fillhandle,'EdgeColor',[0 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    plot(vect,squeeze(est2(channel,:)),'r','LineWidth',3);
    fillhandle = patch([vect fliplr(vect)],[squeeze(CI2(channel,1,:))' fliplr(squeeze(CI2(channel,2,:))')], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    grid on; axis tight; ylabel('Amplitude','FontSize',12); set(gca,'FontSize',12,'layer','top'); box on
    if strcmpi(percent,'mean') || percent == 0
        title(['Means and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    elseif strcmpi(percent,'20% trimmed mean') || percent == 20
        title(['Trimmed Means and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    elseif strcmpi(percent,'median') || percent == 50
        title(['Medians and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    end

    subplot(3,2,[5 6]); hold on
    plot(vect,squeeze(Diff(channel,:,2)),'LineWidth',3);
    fillhandle = patch([vect fliplr(vect)], [squeeze(Diff(channel,:,1)),fliplr(squeeze(Diff(channel,:,3)))], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    grid on; axis tight; xlabel('Time ','FontSize',12)
    ylabel('Amplitude Difference','FontSize',12); set(gca,'FontSize',12,'layer','top'); box on
    if strcmpi(percent,'mean') || percent == 0
        title(['Mean Difference and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    elseif strcmpi(percent,'20% trimmed mean') || percent == 20
        title(['Trimmed Mean Difference and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    elseif strcmpi(percent,'median') || percent == 50
        title(['Median Difference and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    end
end

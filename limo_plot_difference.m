function Data = limo_plot_difference(varargin)

% allows to plot the Difference between two data set with alpha % Hightest
% Density Intervals (ie Bayesian bootstrap CI) 
% WARNING this is computed channel/time frame wise, ie this is not simultaneous CI
%
% FORMATS
% Diff = limo_plot_difference
% Diff = limo_plot_difference(data1,data2,'type','paired/independent')
% Diff = limo_plot_difference(data1,data2,'type','paired/independent',options)
%             'LIMO',LIMOfilename,'percent', 20, 'alpha', 0.05, 'fig', 'on')
%
% INPUTS 
% data1/2 matrices of data (or name of those matrices) 
%         up to 5D when coming from limo_central_tendency_and_ci.m 
% type    'paired' or 'independent'
% LIMO    name of the LIMO file to use 
% options are declared as key value pairs
%         'percent'  is 'mean', 'trimmed mean', 'median' 
%                    or the amount of trimming 0% is a mean, 20% is the default 
%                    trimmed mean, 50% is the median  
%         'alpha'    is the 1-alpha level of the Highest Density Interval
%         'fig'      'on' (default) or 'off' indicates to produce a figure or not
%         'channel'  to indicate the channel to plot
%         'restrict' is 'Time' or 'Frequency' for ERSP analyses
%                    --> only for display, differences and HDI are saved on
%                        drive for both dimensions
%
% OUTPUT
% Data    a structure with
%         .Diff the 3D matrix of Difference with HDI
%         .limo the LIMO structure for channel/time-freq info (optional) 
%
% Cyril Pernet 
% ---------------------------------------------
%  Copyright (C) LIMO Team 2021

%% check inputs

percent     = 20;     % defines the percentage of trimming done
alpha_level = 5/100;  % 1-alpha CI
figure_flag = 1;      % make a figure
channel     = [];     % ask user
restrict    = [];     % for ERSP
wrapdata    = 0;      % for 4D wrap into 3D

if nargin < 4
    
    % select 1st dataset 
    % ------------------
    [file,locpath,idx]=uigetfile('.mat','Select 1st dataset');
    if idx == 0; return; end
    data1 = load(fullfile(locpath,file)); 
    data1 = data1.(cell2mat(fieldnames(data1)));
        
    % select 2nd dataset 
    % ------------------
    [file,locpath,idx]=uigetfile('.mat','Select 1st dataset');
    if idx == 0; return; end
    data2 = load(fullfile(locpath,file)); 
    data2 = data2.(cell2mat(fieldnames(data2)));
       
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
    
    % LIMO
    % ------
    if exist(fullfile(pwd,'LIMO.mat'),'file')
        LIMO = load(fullfile(pwd,'LIMO.mat'));
        LIMO = LIMO.LIMO;
    else
        [file,locpath,idx]=uigetfile('.mat','Select 1st dataset');
        if idx == 0; return; end
        LIMO = load(fullfile(locpath,file)); 
        LIMO = LIMO.LIMO;
    end
           
else % check all inputs
    
    if ischar(varargin{1})
        data1 = load(varargin{1}); 
        data1 = data1.(cell2mat(fieldnames(data1)));
    else
        data1 = varargin{1};
    end

    if ischar(varargin{2})
        data2 = load(varargin{2}); 
        data2 = data2.(cell2mat(fieldnames(data2)));
    else
        data2 = varargin{2};
    end
    
    for n=3:nargin
        if strcmpi(varargin{n},'type')
            type = varargin{n+1};
            if ~strcmpi(type,{'paired','independent'})
                error('unknown value %s for key ''type''',type)
            end
        elseif strcmpi(varargin{n},'percent')
            percent = varargin{n+1};
            if strcmpi(percent,'mean')
                percent = 0;
            elseif strcmpi(percent,'20% trimmed mean')
                percent = 20;
            elseif strcmpi(percent,'median')
                percent = 50;
            end
            if percent < 0
                percent = percent * 100;
            end
        elseif strcmpi(varargin{n},'LIMO')
            LIMO = load(varargin{n+1});
            LIMO = LIMO.LIMO;
        elseif strcmpi(varargin{n},'name')
            name = varargin{n+1};
        elseif strcmpi(varargin{n},'channel')
            channel = varargin{n+1};
        elseif strcmpi(varargin{n},'restrict')
            restrict = varargin{n+1};
        elseif ischar(varargin{n}) % needed to use contains
            if contains(varargin{n},'fig','IgnoreCase',true)
                figure_flag = varargin{n+1};
            elseif contains(varargin{n},'alpha','IgnoreCase',true)
                alpha_level = varargin{n+1};
                if alpha_level > 1
                    alpha_level = alpha_level/100;
                end
            end
        end
    end
end

if ~exist('type','var')
    type = questdlg('are the data','analysis option','paired','independent','paired');
    if isempty(type)
        warning('selection aborded'); return
    end
end

if isfield(data1,'limo'); limo1 = data1.limo; end
if isfield(data1,'data'); data1 = data1.data; end
if isnumeric(data1) 
    fprintf('1st data file loaded \n');
else
    error('couldn''t load the data \n')
end

if isfield(data2,'limo'); limo2 = data2.limo; end
if isfield(data2,'data'); data2 = data2.data; end
if isnumeric(data2)
    fprintf('2nd data file loaded \n');
else
    error('couldn''t load the data \n')
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
            elseif percent == 50
                [est1(channel,:),CI1(channel,:,:)] = limo_central_estimator(squeeze(data1(channel,:,:)),'HD',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:)] = limo_central_estimator(squeeze(data2(channel,:,:)),'HD',1-alpha_level);
                [Diff(channel,:,2),CID]            = limo_central_estimator(squeeze(D(channel,:,:)),'HD',1-alpha_level);
                Diff(channel,:,1)                  = CID(1,:); 
                Diff(channel,:,3)                  = CID(2,:);
            else
                [est1(channel,:),CI1(channel,:,:)] = limo_central_estimator(squeeze(data1(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:)] = limo_central_estimator(squeeze(data2(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [Diff(channel,:,2),CID]            = limo_central_estimator(squeeze(D(channel,:,:)),'Trimmed Mean',1-alpha_level);
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
if ~exist('LIMO','var')
    
    if exist('limo1','var')
        LIMO = limo1;
    elseif exist('limo2','var')
        LIMO = limo2;
    else
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
end

%% save
Data.Diff = Diff;
if exist('LIMO','var')
    Data.limo = LIMO;
end

if ~exist('name','var')
    name      = cell2mat(inputdlg('save as [?]','name option'));
    newname   = sprintf('%s',name);
    if isempty(newname)
        warning  on ; warning('no name - data are not saved')
    else
        newpath   = uigetdir('destination folder?');
        save(fullfile(newpath,newname),'Data');
    end
else
    if isempty(fileparts(name))
        save(fullfile(pwd,name),'Data');
    else
        save(name,'Data');
    end
end

%% Plot
% -----
if strcmp(figure_flag,'on') || figure_flag == 1
    
    if ndims(Diff) == 4
        if isempty(restrict)
            restrict = questdlg('which domain to plot?','showing means only','time','frequency','time');
            if isempty(restrict)
                return
            end
        end
    end
    
    % channel info
    if size(Diff,1) > 1
        if isempty(channel)
            channel = inputdlg('which channel to plot','Plotting option');
            
            if isempty(channel)  % cancel
                return
            elseif strcmp(channel,'') % ok empty
                if ndims(Diff) == 4
                    if strcmpi(restrict ,'time')
                        [~,channel]=max(max(max(squeeze(Diff(:,:,:,2)),[],3)'));
                    else
                        [~,channel]=max(max(squeeze((max(squeeze(Diff(:,:,:,2)),[],2))),[],2));
                    end
                else
                    [~,channel]=max(max(squeeze(Diff(:,:,2)),[],2));
                end
            else
                channel = eval(cell2mat(channel));
            end
        end
        
        if length(channel) > 1
            error('1 channel only can be plotted')
        elseif channel > size(Diff,1)
            error('channel number invalid')
        end
    else
        channel = 1; % should already be one given the loop above
    end
    
    if strcmpi(LIMO.Analysis,'Time')
        if isfield(LIMO.data,'timevect')
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
        if strcmpi(restrict,'Time')
            if isfield(LIMO.data,'tf_times')
                vect = LIMO.data.tf_times;
            else
                vect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;  % in msec
            end
        elseif strcmpi(restrict,'Frequency')
            if isfield(LIMO.data,'tf_freqs')
                vect = LIMO.data.tf_freqs;
            else
                vect = linspace(LIMO.data.low_f,LIMO.data.high_f,size(Diff,2));
            end
        end
    else
        v   = inputdlg('time/freq info missing: enter x-axis interval by hand e.g. [0:0.5:40]');
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
                warning(ME.identifier,'xaxis interval invalid:%s - using default',ME.message)
                vect = 1:size(Diff,2);
            end
        end
    end
    
    % figure
    figure;set(gcf,'Color','w'); 
    subplot(3,2,[1 2 3 4]); 
    if ndims(Diff) == 4
        if strcmpi(restrict ,'time')
            cm = limo_color_images(size(est1,2)*2);
            for f=1:size(est1,2)
                top_plot(vect,squeeze(est1(channel,f,:)),squeeze(est2(channel,f,:)),[],[],cm([f f+size(est1,2)],:),'ss'); 
            end
        else
            cm = limo_color_images(size(est1,3)*2);
            for t=1:size(est1,3)
                top_plot(vect,squeeze(est1(channel,:,t)),squeeze(est2(channel,:,t)),[],[],cm([t t+size(est1,3)],:),'ss');
            end
        end
    else
        top_plot(vect,squeeze(est1(channel,:)),squeeze(est2(channel,:)),...
            squeeze(CI1(channel,:,:)),squeeze(CI2(channel,:,:)),[0 0 1; 1 0 0],'all')
    end
    grid on; axis tight; ylabel('Amplitude','FontSize',12); 
    set(gca,'FontSize',12,'layer','top'); box on
    if strcmpi(percent,'mean') || percent == 0
        Sumstat = 'Means';
    elseif strcmpi(percent,'20% trimmed mean') || percent == 20
        Sumstat = 'Trimmed Means';
    elseif strcmpi(percent,'median') || percent == 50
        Sumstat = 'Medians';
    end
    
    if ndims(Diff) ~= 4
        title([Sumstat ' and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    else
        if strcmpi(restrict,'time')
            title([Sumstat ' at each frequency'],'FontSize',14); drawnow;
        else
            title([Sumstat ' at each time frame'],'FontSize',14); drawnow;
        end
    end
   
    subplot(3,2,[5 6]); 
    if ndims(Diff) == 4
        if strcmpi(restrict ,'time')
            cm = limo_color_images(size(est1,2)); 
            for f=1:size(est1,2)
                bottom_plot(vect,squeeze(Diff(channel,f,:,:)),cm(f,:),'ss')
            end
        else
            cm = limo_color_images(size(est1,3)); 
            for t=1:size(est1,3)
                bottom_plot(vect,squeeze(Diff(channel,:,t,:)),cm(t,:),'ss')
            end
        end
    else
        bottom_plot(vect,squeeze(Diff(channel,:,:)),[0 1 0.2],'all')
    end
    grid on; axis tight; xlabel(LIMO.Analysis,'FontSize',12)
    ylabel('Amplitude Difference','FontSize',12); set(gca,'FontSize',12,'layer','top'); box on
    
    if ndims(Diff) ~= 4
        title([Sumstat(1:end-1) ' difference and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',14); drawnow;
    else
        if strcmpi(restrict,'time')
            title([Sumstat ' differences at each frequency'],'FontSize',14); drawnow;
        else
            title([Sumstat ' differences at each time frame'],'FontSize',14); drawnow;
        end
    end
end
end

function top_plot(vect,est1,est2,CI1,CI2,cm,option)
plot(vect,est1,'LineWidth',3,'Color',cm(1,:)); hold on
if strcmpi(option,'all')
    fillhandle = patch([vect fliplr(vect)],[CI1(1,:) fliplr(CI1(2,:))], cm(1,:));
    set(fillhandle,'EdgeColor',cm(1,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);
end
plot(vect,est2,'LineWidth',3,'Color',cm(2,:));
if strcmpi(option,'all')
    fillhandle = patch([vect fliplr(vect)],[CI2(1,:) fliplr(CI2(2,:))], cm(2,:));
    set(fillhandle,'EdgeColor',cm(2,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);
end
end

function bottom_plot(vect,Diff,cm,option)
plot(vect,Diff(:,2)','LineWidth',3,'Color',cm);hold on
if strcmpi(option,'all')
    fillhandle = patch([vect fliplr(vect)], [Diff(:,1)',flipud(Diff(:,3))'], cm);
    set(fillhandle,'EdgeColor',cm,'FaceAlpha',0.2,'EdgeAlpha',0.8);
end
end

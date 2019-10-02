function Data = limo_plot_difference(varargin)

% allows to plot the difference between two data set with alpha % Hightest
% Density Intervals (ie Bayesian bootstrap CI) 
% WARNING this is computed electrode/time frame wise, ie this is  not simultaneous CI
%
% FORMATS
% diff = limo_plot_difference
% diff = limo_plot_difference(data1,data2,'type','paired/independent')
% diff = limo_plot_difference(data1,data2,'type','paired/independent','percent', 20, 'alpha', 0.05, 'fig', 'on')
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
%         .diff the 3D matrix of difference with HDI
%         .limo the LIMO structure for channel/time-freq info (optional) 
%
% Cyril Pernet - all options in Septembre 2019
% ---------------------------------------------
%  Copyright (C) LIMO Team 2019

%% check inputs

percent     = 20/100; % defines the amount of trimming done
nboot       = 1000;   % bootstrap for CI
alpha_level = 5/100;  % 1-alpha CI
figure_flag = 1;      % make a figure

if nargin < 3 
    % select 1st dataset 
    % ------------------
    [file,locpath,filter]=uigetfile('.mat','Select 1st dataset');
    if filter == 0
        return
    end
    cd(locpath); 
    D = load(file); D = D.Data;
    if isstruct(D)
        name  = fieldnames(D);
        tmp   = sprintf('D.%s',cell2mat(name(1)));
        data1 = eval(tmp);
    else
        data1 = D;    
    end
    data1 = squeeze(data1);
    
    if numel(size(data1)) ~=3
        error('this function only works with 3D matrices');
    else
        fprintf('%s loaded \n',file);
        disp('dimension'); size(data1)
    end
    
    % select 2nd dataset 
    % ------------------
    [file,locpath]=uigetfile('.mat','Select 1st dataset');
    cd(locpath); 
    D = load(file); D = D.Data;
    if isstruct(D)
        name = fieldnames(D);
        tmp = sprintf('D.%s',cell2mat(name(1)));
        data2 = eval(tmp);
    else
        data2 = D;
    end
    
    data2 = squeeze(data2);
    if numel(size(data2)) ~=3
        error('this function only works with 3D matrices');
    else
        fprintf('%s loaded \n',file);
        disp('dimension'); size(data2)
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
if strcmpi(type,'paired')
    if numel(size(data1)) ~= numel(size(data2))
        error('for paired data, each matrix must be of the same size')
    elseif numel(size(data1)) == 3
        if sum(size(data1) == size(data2)) ~= 3
            error('for paired data, each matrix must be of the same size')
        end
    elseif numel(size(data1)) == 4
        if sum(size(data1) == size(data2)) ~= 4
            error('for paired data, each matrix must be of the same size')
        end
    end
else
    if size(data1,1) ~= size(data2,1) || size(data1,2) ~= size(data2,2)
        error('the number of elements in dim 1 or 2 don''t match');
    end
end

diff = NaN(size(data1,1),size(data1,2),3);

%% compute the mean difference and CI

if strcmpi(type,'Paired')
    
        if numel(size(data1)) == 4
            data1 = squeeze(data1(:,:,:,2));
            data2 = squeeze(data2(:,:,:,2));
        end
 
        % the difference 
        D = data1-data2;
        
        for channel=1:size(data1,1)
            fprintf('bootstraping data for CI estimation electrode %g \n',channel)
            if strcmpi(percent,'mean') || percent == 0
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'Mean',1-alpha_level);
                [diff(channel,:,2),CID] = limo_central_estimator(squeeze(D(channel,:,:)),'Mean',1-alpha_level);
                diff(channel,:,1) = CID(1,:); diff(channel,:,3) = CID(2,:);
            elseif strcmpi(percent,'20% trimmed mean') || percent == 20
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [diff(channel,:,2),CID] = limo_central_estimator(squeeze(D(channel,:,:)),'Trimmed Mean',1-alpha_level);
                diff(channel,:,1) = CID(1,:); diff(channel,:,3) = CID(2,:);
            elseif strcmpi(percent,'median') || percent == 50
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'HD',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'HD',1-alpha_level);
                [diff(channel,:,2),CID] = limo_central_estimator(squeeze(D(channel,:,:)),'HD',1-alpha_level);
                diff(channel,:,1) = CID(1,:); diff(channel,:,3) = CID(2,:);
            end
            
        end
        
        % --------------------------------------------------------
        
elseif strcmpi(type,'Independent')
        
        if numel(size(data1)) == 4
            data1 = squeeze(data1(:,:,:,2));
            data2 = squeeze(data2(:,:,:,2));
        end
        
        est1 =NaN(size(data1,1),size(data1,2)); est2 = est1;
        CI1 =NaN(size(data1,1),2,size(data1,2)); CI2 = CI1;
        for channel=1:size(data1,1)
            fprintf('bootstraping data for CI estimation electrode %g \n',channel)
            
            if strcmpi(percent,'mean') || percent == 0
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'Mean',1-alpha_level);
                diff(channel,:,2) = est1(channel,:)-est2(channel,:);
            elseif strcmpi(percent,'20% trimmed mean') || percent == 20
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'Trimmed Mean',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'Trimmed Mean',1-alpha_level);
                diff(channel,:,2) = est1(channel,:)-est2(channel,:);
            elseif strcmpi(percent,'median') || percent == 50
                [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'HD',1-alpha_level);
                [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'HD',1-alpha_level);
                diff(channel,:,2) = est1(channel,:)-est2(channel,:);
            end
            
            sorted_data   = sort(sort(bb1)-sort(bb2),2);
            upper_centile = floor((1-alpha_level)*size(sorted_data,2)); % upper bound
            nCIs          = size(sorted_data,2) - upper_centile;
            for frame = 1:size(sorted_data,1)
                tmp = sorted_data(frame,:);
                ci = 1:nCIs; ciWidth = tmp(ci+upper_centile) - tmp(ci); % all centile distances
                [~,index]=find(ciWidth == min(ciWidth)); % densest centile
                if length(index) > 1; index = index(1); end % many similar values
                diff(channel,frame,1) = tmp(index);
                diff(channel,frame,3) = tmp(index+upper_centile);
            end
        end
end

%% Plot
% -----
if strcmp(figure_flag,'on') || figure_flag == 1
    % electrode info
    if size(diff,1) > 1
        channel = inputdlg('which electrode to plot','Plotting option');
        if isempty(channel)  % cancel
            return
        elseif strcmp(channel,'') % ok empty
            [~,channel]=max(max(squeeze(diff(:,:,2)),[],2));
        else
            channel = eval(cell2mat(channel));
            if length(channel) > 1
                error('1 electrode only can be plotted')
            elseif channel > size(diff,1)
                error('electrode number invalid')
            end
        end
    else
        channel = 1;
    end
    
    % time/freq info
    [file,locpath,ind]=uigetfile({'LIMO.mat'},'Select any LIMO with right info');
    if strcmp(file,'LIMO.mat')
        load(fullfile(locpath,file));
        if strcmpi(LIMO.Analysis,'Time')
            vect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
        elseif strcmpi(LIMO.Analysis,'Frequency')
            vect = LIMO.data.freqlist;
            if size(vect,2) == 1; vect = vect'; end
            if size(vect,2) ~= size(toplot,2)
                vect = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
            end
        end
    else
        warning('selection aborded'); return
    end
    
    if  single(ind == 0) || length(vect) ~= length(squeeze(est1(channel,:)))
        if length(vect) ~= squeeze(est1(channel,:))
                fprintf('error in computing %s frames \n',LIMO.Analysis)
        end
        dlg = sprintf('enter %s interval by hand e.g. [0:0.5:40]',LIMO.Analysis);
        v = inputdlg(dlg);
        if isempty(v)
            return
        else
            try
                vect = eval(cell2mat(v));
                if length(vect) ~= size(diff,2)
                    fprintf('%s interval invalid format \n',LIMO.Analysis)
                    vect = 1:size(diff,2);
                end
            catch ME
                fprintf('%s interval invalid format \n',LIMO.Analysis)
                vect = 1:size(diff,2);
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
        title(['Means and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',16); drawnow;
    elseif strcmpi(percent,'20% trimmed mean') || percent == 20
        title(['Trimmed Means and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',16); drawnow;
    elseif strcmpi(percent,'median') || percent == 50
        title(['Medians and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',16); drawnow;
    end

    subplot(3,2,[5 6]); hold on
    plot(vect,squeeze(diff(channel,:,2)),'LineWidth',3);
    fillhandle = patch([vect fliplr(vect)], [squeeze(diff(channel,:,1)),fliplr(squeeze(diff(channel,:,3)))], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    grid on; axis tight; xlabel('Time ','FontSize',12)
    ylabel('Amplitude Difference','FontSize',12); set(gca,'FontSize',12,'layer','top'); box on
    if strcmpi(percent,'mean') || percent == 0
        title(['Mean difference and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',16); drawnow;
    elseif strcmpi(percent,'20% trimmed mean') || percent == 20
        title(['Trimmed Mean difference and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',16); drawnow;
    elseif strcmpi(percent,'median') || percent == 50
        title(['Median difference and ' num2str(100-alpha_level*100) '%HDI'],'FontSize',16); drawnow;
    end

end

%% save
name    = cell2mat(inputdlg('save as [?]','name option'));
newname = sprintf('%s',name);
newpath = uigetdir('destination folder?');

Data.diff = diff;
if exist('LIMO','var')
    Data.limo = LIMO;
end
save (fullfile(newpath,newname),'Data');

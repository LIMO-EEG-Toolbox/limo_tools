function diff = limo_plot_difference(varargin)

% allows to plot the difference between two data set with alpha_level% confidence
% interval - WARNING this is computed electrode/time frame wise, ie this is 
% not simultaneous CI (alpha_level level is not controlled over all differences)
%
% FORMAT
% diff = limo_plot_difference(data1,data2,type)
% diff = limo_plot_difference(data1,data2,type,alpha_level)
% diff = limo_plot_difference(data1,data2,type,alpha_level,flag)
%
% INPUT
% data1/2 matrices of data up to 4D, if 4D because coming from robust
% averaging, i.e. the last dimension is the estimator with CI, only the
% estimator is used ie level 2 if dim 4
% type [1/2] 1 for paired data and 2 for independent data
% alpha_level is the alpha_level level
% flag [0/1] indicate to produce a figure or not
%
% OUPUT
% diff is a 3D matrix of difference with CI
%
% Cyril Pernet 01/11/2011
% -----------------------------
% Copyright (C) LIMO Team 2010

%% check inputs

percent = 20/100; % defines the amount of trimming done
nboot = 1000;
alpha_level = 5/100;
figure_flag = 1;
if nargin < 3 
    % data1
    [file,locpath,filter]=uigetfile('.mat','Select 1st dataset');
    if filter == 0
        return
    end
    cd(locpath); 
    D = load(file); D = D.Data;
    if isstruct(D)
        name = fieldnames(D);
        tmp = sprintf('D.%s',cell2mat(name(1)));
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
    % data2
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
    answer = questdlg('are the data','analysis option','paired','independent','paired');
    if strcmp(answer,'paired')
        type = 1;
    else
        type = 2;
    end
    % alpha_level
    v = inputdlg('set alpha_level level (%) for confidence intervals');
    if isempty(v)
        return
    else
        alpha_level = eval(cell2mat(v));
        if alpha_level > 1
            alpha_level = alpha_level / 100;
        end
    end    
elseif nargin>= 3 && nargin <=5
    data1       = varargin{1};
    data2       = varargin{2};
    type        = varargin{3};
    if nargin >= 4
        alpha_level = varargin{4}; 
        if alpha_level > 1
            alpha_level = alpha_level / 100;
        end
    end
    if nargin == 5
        figure_flag = varargin{5};
    end
elseif nargin > 5
    error('too many arguments')
end

if alpha_level == 0
    disp('stange value for alpha - adjusted to 5%');
    alpha_level = 5/100;
end

% check possible dimensions issues
if type == 1
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

switch type
    
    case {1} % paired data
        
        if numel(size(data1)) == 4
            data1 = squeeze(data1(:,:,:,2));
            data2 = squeeze(data2(:,:,:,2));
        end
 
        % the difference 
        D = data1-data2;
        
        for channel=1:size(data1,1)
            fprintf('bootstraping data for CI estimation electrode %g \n',channel)
            [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'trimmed mean',1-alpha_level);
            [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'trimmed mean',1-alpha_level);
            
            [diff(channel,:,2),CID] = limo_central_estimator(squeeze(D(channel,:,:)),'trimmed mean',1-alpha_level);
            diff(channel,:,1) = CID(1,:);
            diff(channel,:,3) = CID(2,:);
        end
        
        % --------------------------------------------------------
        
    case {2} % independent data
        
        if numel(size(data1)) == 4
            data1 = squeeze(data1(:,:,:,2));
            data2 = squeeze(data2(:,:,:,2));
        end
        
        est1 =NaN(size(data1,1),size(data1,2)); est2 = est1;
        CI1 =NaN(size(data1,1),2,size(data1,2)); CI2 = CI1;
        for channel=1:size(data1,1)
            fprintf('bootstraping data for CI estimation electrode %g \n',channel)
            
            [est1(channel,:),CI1(channel,:,:),bb1] = limo_central_estimator(squeeze(data1(channel,:,:)),'trimmed mean',1-alpha_level);
            [est2(channel,:),CI2(channel,:,:),bb2] = limo_central_estimator(squeeze(data2(channel,:,:)),'trimmed mean',1-alpha_level);
            diff(channel,:,2) = est1(channel,:)-est2(channel,:);

            sorted_data = sort(sort(bb1)-sort(bb2),2);
            upper_centile = floor((1-alpha_level)*size(sorted_data,2)); % upper bound
            nCIs = size(sorted_data,2) - upper_centile;
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

% save
name = cell2mat(inputdlg('save as [?]','name option'));
newname = sprintf('%s',name);
newpath = uigetdir('destination folder?');
cd(newpath); save (newname,'diff');


% plot
% -----
if figure_flag == 1
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
    % timing info
    [file,locpath,ind]=uigetfile('.mat','Select any LIMO with right timing info');
    if ind ==1
        cd(locpath); load LIMO;
        timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
    end
    
    if  single(ind == 0) || length(timevect) ~= length(squeeze(est1(channel,:)))
        if length(timevect) ~= squeeze(est1(channel,:))
            disp('error in computing time frames')
        end
        v = inputdlg('enter time interval by hand e.g. [0:0.5:200]');
        if isempty(v)
            return
        else
            try
                timevect = eval(cell2mat(v));
                if length(timevect) ~= size(diff,2)
                    disp('time interval invalid format');
                    timevect = 1:size(diff,2);
                end
            catch ME
                disp('time interval invalid format');
                timevect = 1:size(diff,2);
            end
        end
    end
    
    
    % figure
    figure;set(gcf,'Color','w'); subplot(3,2,[1 2 3 4]); hold on
    plot(timevect,squeeze(est1(channel,:)),'LineWidth',3);
    fillhandle = patch([timevect fliplr(timevect)],[squeeze(CI1(channel,1,:))' fliplr(squeeze(CI1(channel,2,:))')], [0 0 1]);
    set(fillhandle,'EdgeColor',[0 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    plot(timevect,squeeze(est2(channel,:)),'r','LineWidth',3);
    fillhandle = patch([timevect fliplr(timevect)],[squeeze(CI2(channel,1,:))' fliplr(squeeze(CI2(channel,2,:))')], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    grid on; axis tight; ylabel('Amplitude','FontSize',12); set(gca,'FontSize',12,'layer','top'); box on
    title(['Trimmed Means and ' num2str(100-alpha_level*100) '%CI'],'FontSize',16); drawnow;

    subplot(3,2,[5 6]); hold on
    plot(timevect,squeeze(diff(channel,:,2)),'LineWidth',3);
    fillhandle = patch([timevect fliplr(timevect)], [squeeze(diff(channel,:,1)),fliplr(squeeze(diff(channel,:,3)))], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    grid on; axis tight; xlabel('Time ','FontSize',12)
    ylabel('Amplitude Difference','FontSize',12); set(gca,'FontSize',12,'layer','top'); box on
    title(['Trimmed Mean difference and ' num2str(100-alpha_level*100) '%CI'],'FontSize',16); drawnow;

end


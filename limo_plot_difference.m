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
    [file,locpath]=uigetfile('.mat','Select 1st dataset');
    cd(locpath); D = load(file);
    if isstruct(D)
        name = fieldnames(D);
        tmp = sprintf('D.%s',cell2mat(name));
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
    cd(locpath); D = load(file);
    if isstruct(D)
        name = fieldnames(D);
        tmp = sprintf('D.%s',cell2mat(name));
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

% univariate bounds
diff = NaN(size(data1,1),size(data1,2),3);
low = round(nboot.*alpha_level./2);
high = nboot - low;

%% compute the mean difference and CI

switch type
    
    case {1} % paired data
        
        if numel(size(data1)) == 4
            data1 = squeeze(data1(:,:,:,2));
            data2 = squeeze(data2(:,:,:,2));
        end
 
        na=size(data1,3);
        ga=floor(percent*na);
        data1=sort(data1,3);
        data1=data1(:,:,(ga+1):(na-ga));

        na=size(data2,3);
        ga=floor(percent*na);
        data2=sort(data2,3);
        data2=data2(:,:,(ga+1):(na-ga));

        % the difference 
        D = data1-data2;
        
        % compute bootstrap differences
        % ------------------------------
        boot_table = zeros(size(data1,3),nboot); B=1;
        while B~=nboot+1
            tmp = randi(size(data1,3),size(data1,3),1);
            if length(unique(tmp)) ~= 1
                boot_table(:,B) = tmp;
                B=B+1;
            end
        end
        
        disp('bootstraping data for CI estimation ...')
        loop = nboot / 100; index1 = 1; index2 = 100; 
        for b=1:loop % rather than doing b loops we work by blocks of 100
            D_boot = D(:,:,boot_table(:,index1:index2));
            avg_boot(:,:,index1:index2) = squeeze(nanmean(reshape(D_boot,size(data1,1),size(data1,2),size(data1,3),100),3));
            index1 = index1+100; index2 = index2+100;
        end
        
        avg_boot = sort(avg_boot,3);
        % --------------------
        diff(:,:,1) = avg_boot(:,:,low);
        diff(:,:,2) = nanmean(D,3);
        diff(:,:,3) = avg_boot(:,:,high);
        
        TM1 = NaN(size(data1,1),size(data1,2),nboot);
        TM2 = NaN(size(data2,1),size(data2,2),nboot);
        for electrode=1:size(data1,1)
            % means over resampled subjects
            TM1(electrode,:,:) = squeeze(nanmean(reshape(data1(electrode,:,boot_table),size(data1,2),size(data1,3),nboot),2));
            TM2(electrode,:,:) = squeeze(nanmean(reshape(data2(electrode,:,boot_table),size(data2,2),size(data2,3),nboot),2));
        end
        TM1 = sort(TM1,3);
        TM2 = sort(TM2,3);
        
        % --------------------------------------------------------
        
    case {2} % independent data
        
        na=size(data1,3);
        ga=floor(percent*na);
        data1=sort(data1,3);
        data1=data1(:,:,(ga+1):(na-ga));

        boot_table1 = zeros(size(data1,3),nboot); B=1;
        while B~=nboot+1
            tmp = randi(size(data1,3),size(data1,3),1);
            if length(unique(tmp)) ~= 1
                boot_table1(:,B) = tmp;
                B=B+1;
            end
        end
        
        na=size(data2,3);
        ga=floor(percent*na);
        data2=sort(data2,3);
        data2=data2(:,:,(ga+1):(na-ga));

        boot_table2 = zeros(size(data2,3),nboot); B=1;
        while B~=nboot+1
            tmp = randi(size(data2,3),size(data2,3),1);
            if length(unique(tmp)) ~= 1
                boot_table2(:,B) = tmp;
                B=B+1;
            end
        end
        
        TM1 = NaN(size(data1,1),size(data1,2),nboot);
        TM2 = NaN(size(data2,1),size(data2,2),nboot);
        for electrode=1:size(data1,1)
            fprintf('bootstraping data for CI estimation electrode %g \n',electrode)
            % means over resampled subjects
            means_gp1 = squeeze(nanmean(reshape(data1(electrode,:,boot_table1),size(data1,2),size(data1,3),nboot),2));
            means_gp2 = squeeze(nanmean(reshape(data2(electrode,:,boot_table2),size(data2,2),size(data2,3),nboot),2));
            D = sort(means_gp1-means_gp2,2);
            
            TM1(electrode,:,:) = means_gp1;
            TM2(electrode,:,:) = means_gp2;
            % --------------------
            diff(electrode,:,1) = D(:,low);
            diff(electrode,:,2) = nanmean(squeeze(data1(electrode,:,:)),2)-nanmean(squeeze(data2(electrode,:,:)),2);
            diff(electrode,:,3) = D(:,high);
        end
        TM1 = sort(TM1,3);
        TM2 = sort(TM2,3);
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
        electrode = inputdlg('which electrode to plot','Plotting option');
        if isempty(electrode)
            return
        else
            electrode = eval(cell2mat(electrode));
            if length(electrode) > 1
                error('1 electrode only can be plotted')
            elseif electrode > size(diff,1)
                error('electrode number invalid')
            end
        end
    else
        electrode = 1;
    end
    % timing info
    [file,locpath,ind]=uigetfile('.mat','Select any LIMO with right timing info');
    if ind == 0
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
    else
        cd(locpath); load LIMO;
        timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
    end
    
    % figure
    figure;set(gcf,'Color','w'); hold on
    plot(timevect,nanmean(squeeze(data1(electrode,:,:)),2),'LineWidth',3);
    fillhandle = patch([timevect fliplr(timevect)], [squeeze(TM1(electrode,:,low)),fliplr(squeeze(TM1(electrode,:,high)))], [0 0 1]);
    set(fillhandle,'EdgeColor',[0 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    plot(timevect,nanmean(squeeze(data2(electrode,:,:)),2),'r','LineWidth',3);
    fillhandle = patch([timevect fliplr(timevect)], [squeeze(TM2(electrode,:,low)),fliplr(squeeze(TM2(electrode,:,high)))], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    grid on; axis tight; box on; 
    xlabel('Time ','FontSize',16)
    ylabel('Amplitude','FontSize',16)
    title(['Trimmed Means and ' num2str(100-alpha_level*100) '%CI'],'FontSize',18); drawnow;
    set(gca,'FontSize',14,'layer','top');

    figure;set(gcf,'Color','w'); hold on
    plot(timevect,squeeze(diff(electrode,:,2)),'LineWidth',3);
    fillhandle = patch([timevect fliplr(timevect)], [squeeze(diff(electrode,:,1)),fliplr(squeeze(diff(electrode,:,3)))], [1 0 0]);
    set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
    grid on; axis tight; box on; 
    xlabel('Time ','FontSize',16)
    ylabel('Amplitude','FontSize',16)
    title(['Trimmed Mean difference and ' num2str(100-alpha_level*100) '%CI'],'FontSize',18); drawnow;
    set(gca,'FontSize',14,'layer','top');

end


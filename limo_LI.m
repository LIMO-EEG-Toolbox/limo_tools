function LI_stats = limo_LI(varargin)

% routine to compute the LI for thresholded maps on 'file'
% this assumes that TFCE was computed - at each thresold one computes
% left-right / left+right then integrate the results ; H0 is computing
% permuting left/right channels ; the bias is computed using the null
% bootrapped/TFCE data (if present)
%
% FORMAT LI_stats = limo_LI(file)
%        LI_stats = limo_LI(file,'channelpairs',channelindices,'summary','single','alpha',0.05,'figure','on')
%
% INPUT is any known LIMO stat files
%       'channelpairs' key matches the channelindices values, a n*2 matrix for chanloc to pair
%       'summary' can be 'single' (default) or 'split' reflecting the
%                 summary statistic for the lateralization curve, a single mean or split negative  
%                 and positive values of the map (for instance for a t-test between 2 conditions)   
%       'alpha' key matches the alpha level confience interval value
%       'figure' is 'on' by default or switch it 'off'
%
%  'channelpairs' is optional but it is strongly advised to provide it
%  making sure the correct channels are paired -- this is done easily using
%  limo_pair_channels.m
%
% OUTPUT LI_stats is a structure with the following fields:
%     LI          -- the lateralization curve based on spatial-temporal clustering with tfce
%     thresholds  -- the thresholds used in the tfce computation
%     mean        -- the LI mean (split between - and + for t-tests)
%     H0          -- the null LI means under permutation of channels
%     H0CI        -- the alphav conficence intervals testing mean LI
%     bias        -- the LI means under null data
%     biasCI      -- the alphav conficence intervals of null data (should be symmetric)
%
% Note if a t-test or con map is used, this likely contains negative and
% positive values ; the output is thus split and alpha is divided by 2
%
% see also limo_pair_channels.m limo_lateralization.m
%
% Cyril Pernet v1 25-03-2021
% ------------------------------
%  Copyright (C) LIMO Team 2021

%% input check
[filepath,FileName]=fileparts(varargin{1});
if isempty(filepath)
    filepath = pwd;
end

% defaults
alphav     = 0.05;
summary    = 'single';
fig_option = 'on';

% options
if nargin >1
    for n=1:nargin
        if strcmpi(varargin{n},'alpha')
            alphav = varargin{n+1};
            if alphav > 1
                alphav = alphav / 100;
            end
        elseif contains(varargin{n},'channelpair')
            channels = varargin{n+1};
        elseif strcmpi(varargin{n},'summary')
            if any(contains(varargin{n+1},{'single','split'}))
                summary = varargin{n+1};
            else
               error('the option ''channelpairs'' must be with the value ''single'' or ''split''') 
            end
        elseif strcmpi(varargin{n},'figure')
            fig_option = varargin{n+1};
        end
    end
end

% check metadata
if ~exist(fullfile(filepath,'LIMO.mat'))
    warning('no LIMO file found')
    [Name,Path,filter] = uigetfile('LIMO.mat','select LIMO.mat file');
    if isempty(filter) || filter == 0
        warning('selection aborded')
        return
    else
        LIMO = load(fullfile(Path,Name));
        LIMO = LIMO.(cell2mat(fieldnames(LIMO)));
    end
else
    LIMO = load(fullfile(filepath,'LIMO.mat'));
    if isfield(LIMO,'LIMO')
        LIMO = LIMO.LIMO;
    else
        error('can''t load LIMO file')
    end
end

% check file
M = [];
matfile = load(fullfile(LIMO.dir,FileName));
if strcmpi(LIMO.Analysis,'Time-Frequency')
    if strcmp(FileName,'R2.mat')
        M         = squeeze(matfile.R2(:,:,:,2));
    elseif strncmp(FileName,'Condition_effect',16)
        M         = squeeze(matfile.Condition_effect(:,:,:,1));
    elseif strncmp(FileName,'Covariate_effect',16)
        M         = squeeze(matfile.Covariate_effect(:,:,:,1));
    elseif strncmp(FileName,'Interaction_effect',18)
        M         = squeeze(matfile.Interaction_effect(:,:,:,1));
    elseif strncmp(FileName,'semi_partial_coef',17)
        M         = squeeze(matfile.semi_partial_coef(:,:,:,2));
    elseif strncmp(FileName,'con_',4)
        M         = squeeze(matfile.con(:,:,:,4));
    elseif contains(FileName,'ttest')
        matfile   = matfile.(cell2mat(fieldnames(matfile)));
        M         = matfile(:,:,:,4);
    elseif strncmp(FileName,'ess_',4)
        M         = squeeze(matfile.ess(:,:,:,end-1));
    elseif contains(FileName,'Rep_ANOVA')
        M         = matfile.(cell2mat(fieldnames(matfile)))(:,:,:,1);
    end
    
else  % same with one less dimention
    if strcmp(FileName,'R2.mat')
        M         = squeeze(matfile.R2(:,:,2));
    elseif strncmp(FileName,'Condition_effect',16)
        M         = squeeze(matfile.Condition_effect(:,:,1));
    elseif strncmp(FileName,'Covariate_effect',16)
        M         = squeeze(matfile.Covariate_effect(:,:,1));
    elseif strncmp(FileName,'Interaction_effect',18)
        M         = squeeze(matfile.Interaction_effect(:,:,1));
    elseif strncmp(FileName,'semi_partial_coef',17)
        M         = squeeze(matfile.semi_partial_coef(:,:,2));
    elseif strncmp(FileName,'con_',4)
        M         = squeeze(matfile.con(:,:,4));
    elseif contains(FileName,'ttest')
        matfile   = matfile.(cell2mat(fieldnames(matfile)));
        M         = matfile(:,:,4);
    elseif strncmp(FileName,'ess_',4)
        M         = squeeze(matfile.ess(:,:,end-1));
    elseif contains(FileName,'Rep_ANOVA')
        M         = matfile.(cell2mat(fieldnames(matfile)))(:,:,1);
    end
end

if isempty(M)
    if exist('errordlg2','file')
        errordlg2('File not supported'); return
    else
        errordlg('File not supported'); return
    end
else
    MCC_FileName  = fullfile(LIMO.dir,['H0' filesep 'H0_' FileName]);
end

%% compute

% do what limo_random_robust would do but get all thresholded tfce maps
% and use this to compare left vs. right, then integrate
[~,all_maps] = limo_tfce_handling(fullfile(LIMO.dir,FileName),'checkfile','no');
if iscell(all_maps)
    th_maps = all_maps{1};
    H0_maps = all_maps{2};
else
    th_maps = all_maps;
end
clear all_maps

% LI curve is the LI per threshold
% A tfce map contains the weighted clusters over space and frames, we just do
% the mean tfce values (these are the ones surviving a given threshold) for
% left and right hemispheres and compute LI. Just like for a tfce score,
% the final result inegrate over thresholds - here we use a 20% trimmed
% mean to avoid to large/small values, giving a LI over time.

if ~exist('channel','var')
    channels = limo_pair_channels(LIMO.data.chanlocs);
end
LI         = getLI(LIMO,channels,th_maps);
thresholds = linspace(min(M(:)),max(M(:)),size(th_maps,numel(size(th_maps))));
if any(thresholds<0) && strcmpi(summary,'split')
    meanLI(1) = nanmean(LI(thresholds<0));
    meanLI(2) = nanmean(LI(thresholds>0));
    H0_LI = NaN(2,1000);
else
    meanLI = mean(LI);
    H0_LI = NaN(1,1000);
end
LI_stats.LI         = LI;
LI_stats.thresholds = thresholds;
LI_stats.mean       = meanLI;

% LI curve null is how is no topographical organization
% Use permutation to test random left vs right pairs
fprintf('Estimating the null topography ...\n')
for p =1:1000
    LI = NaN(1,size(th_maps,numel(size(th_maps))));
    parfor th = 1:size(th_maps,numel(size(th_maps)))
        % invert left/right for an arbritrary number of channels at random
        RSelect                = randperm(27);
        RShitf                 = randperm(27);
        RSelect                = RSelect(1:RShitf(1));  % take a random subset
        Rchannels              = channels;
        Rchannels(RSelect,:)   = [channels(RSelect,2) channels(RSelect,1)]; % invert left/right for this subset
        left                   = squeeze(th_maps(Rchannels(:,1),:,th));
        right                  = squeeze(th_maps(Rchannels(:,2),:,th));
        left                   = trimmean(left(left~=0),40);
        if isnan(left); left   = 0; end
        right                  = trimmean(right(right~=0),40);
        if isnan(right); right = 0; end
        LI(th)                 = ((left-right)./(left+right)).*100;
    end
    
    if any(thresholds<0) && strcmpi(summary,'split')
        H0_LI(1,p) = nanmean(LI(thresholds<0));
        H0_LI(2,p) = nanmean(LI(thresholds>0));
    else
        H0_LI(p) = nanmean(LI);
    end
end

LI_stats.H0   = sort(H0_LI,2);
if size(LI_stats.H0,1) == 1
    tmp                = LI_stats.H0;
    if ~isempty(isnan(tmp))
        tmp(isnan(tmp))    = [];
    end
    low                = round(alphav*length(tmp));
    high               = length(tmp) - low;
    LI_stats.H0CI      = [tmp(low) tmp(high)];
else
    tmp                = LI_stats.H0(1,:);
    if ~isempty(isnan(tmp))
        tmp(isnan(tmp))    = [];
    end
    low                = round(alphav/2*length(tmp));
    high               = floor(length(tmp) - low);
    LI_stats.H0CI(:,1) = [tmp(low) tmp(high)];
    tmp                = LI_stats.H0(2,:);
    if ~isempty(isnan(tmp))
        tmp(isnan(tmp))    = [];
    end
    low                = round(alphav/2*length(tmp));
    high               = floor(length(tmp) - low);
    LI_stats.H0CI(:,2) = [tmp(low) tmp(high)];
end

%% figure
if strcmpi(fig_option,'on')
    hfig = figure;
    subplot(3,6,1:12); plot(LI_stats.thresholds,LI_stats.LI,'LineWidth',3); grid on; box on
    xlabel('TFCE thresholds'); ylabel('lateralization index'); ax = get(gca);
    if length(LI_stats.mean) == 1
        title(sprintf('Lateralization curve - mean LI %g ',LI_stats.mean))
        if isfield(LI_stats,'H0CI')
            hold on; plot(LI_stats.thresholds,repmat(LI_stats.H0CI(1),1,length(LI_stats.thresholds)),'k--','LineWidth',2)
            plot(LI_stats.thresholds,repmat(LI_stats.H0CI(2),1,length(LI_stats.thresholds)),'k--','LineWidth',2)
        end
    else
        title(sprintf('Lateralization curve \n negative map mean LI %g positve map mean LI %g ',LI_stats.mean(1),LI_stats.mean(2)))
        if isfield(LI_stats,'H0CI')
            hold on; plot(LI_stats.thresholds,repmat(LI_stats.H0CI(1,:)',1,length(LI_stats.thresholds)),'k--','LineWidth',2)
            plot(LI_stats.thresholds,repmat(LI_stats.H0CI(2,:)',1,length(LI_stats.thresholds)),'k--','LineWidth',2)
        end
    end
    axis([LI_stats.thresholds(1)-0.1 LI_stats.thresholds(end)+0.1 ax.YAxis.TickValues(1)-1 ax.YAxis.TickValues(end)+1]);
    opt = {'maplimits','absmax','electrodes','off','verbose','off','colormap', ...
        limo_color_images(trimmean(th_maps(:,:,:),40,3))};
    frames = round(linspace(1,length(LI_stats.thresholds),6));
    for f=1:length(frames)
        subplot(3,6,12+f);
        topoplot(sum(squeeze(th_maps(:,:,frames(f))),2),LIMO.data.chanlocs,opt{:});
    end
end

%% addtional bias analysis

% there is no effect, i.e. under H0, the tfce maps left vs right
% should be symetric - with chance at 50%
if exist('H0_maps','var')
    fprintf('Testing bias using null data ...\n')
    H0_data = load(MCC_FileName);
    H0_data = H0_data.(cell2mat(fieldnames(H0_data)));
    if contains(FileName,'R2') || contains(FileName,'semi_partial')
        bootM = squeeze(H0_data(:,:,2,:)); % get all F values under H0
    else
        bootM = squeeze(H0_data(:,:,1,:));
    end
    
    clear thresholds
    parfor m=1:length(H0_maps)
        H0_curve{m} = getLI(LIMO,channels,H0_maps{m});
        H0M = bootM(:,:,m);
        thresholds{m} = linspace(min(H0M(:)),max(H0M(:)),size(H0_maps{m},numel(size(H0_maps{m}))));
    end
    
    for m=1:length(H0_maps)
        if any(thresholds{m}<0)
            H0_LI(1,m) = nanmean(H0_curve{m}(thresholds{m}<0));
            H0_LI(2,m) = nanmean(H0_curve{m}(thresholds{m}>0));
        else
            H0_LI(m) = nanmean(H0_curve{m});
        end
    end
    
    LI_stats.bias            = sort(H0_LI,2);
    if size(LI_stats.bias,1) == 1
        tmp                  = LI_stats.bias;
        if ~isempty(isnan(tmp))
            tmp(isnan(tmp))    = [];
        end
        low                  = round(alphav*length(tmp));
        high                 = length(tmp) - low;
        LI_stats.biasCI      = [tmp(low) tmp(high)];
    else
        tmp                  = LI_stats.bias(1,:);
        if ~isempty(isnan(tmp))
            tmp(isnan(tmp))    = [];
        end
        low                  = round(alphav/2*length(tmp));
        high                 = length(tmp) - low;
        LI_stats.biasCI(:,1) = [tmp(low) tmp(high)];
        tmp                  = LI_stats.bias(2,:);
        if ~isempty(isnan(tmp))
            tmp(isnan(tmp))    = [];
        end
        low                  = round(alphav/2*length(tmp));
        high                 = length(tmp) - low;
        LI_stats.biasCI(:,2) = [tmp(low) tmp(high)];
    end
end

%% routine to compute LI

function LI = getLI(LIMO,channels,th_maps)

if strcmpi(LIMO.Analysis,'Time-Frequency')
    disp('to do')
else
    for th = size(th_maps,numel(size(th_maps))):-1:1
        left   = squeeze(th_maps(channels(:,1),:,th));
        right  = squeeze(th_maps(channels(:,2),:,th));
        left   = trimmean(left(left~=0),40);
        if isnan(left); left = 0; end
        right  = trimmean(right(right~=0),40);
        if isnan(right); right = 0; end
        LI(th) = ((left-right)./(left+right)).*100;
    end
end



function res = limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,varargin)

% This function displays various results
% The arguments specify cases for the
% different kind of figures, thresholds etc ..
%
% FORMAT:
% limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,options)
%
% INPUTS:
%   Type      = type of images/plot to do
%               1 - 2D images with a intensity plotted as function of time (x) and electrodes (y)
%               2 - topographic plot a la eeglab
%               3 - plot the ERP data (original or modeled)
%   Filename  = Name of the file to image
%   PathName  = Path of the file to image
%   p         = threshold p value e.g. 0.05
%   MCC       = Multiple Comparison technique
%               1=None, 2= Cluster, 3=TFCE, 4=T max
%   LIMO      = LIMO structure
%   flag      = interactivity  (1) or not (0)
%
% OPTIONAL INPUTS  (Usage: {''key'', value, ... })
% 'channels' : Provide the index of the channel to be used.
% 'regressor': Provide the index of the regressor to be used.
% 'plot3type': Type of plots to show when 'Type' is 3. Select between {'Original', 'Modeled', 'Adjusted'}
% 'sumstats' : Course plot summary statistics 'Mean' or 'Trimmed'
% 'restrict' : for time-frequency data, plot restrict plot to 'Time' or 'Frequency'
% 'dimvalue' : for time-frequency data, what value to resctrict on (e.g. restrict to 'Time' with dimvalue 5Hz)
%
% Although the function is mainly intented to be used via the GUI, some figures
% can be generated automatically, for instance limo_display_results(1,'R2.mat',pwd,0.05,5,LIMO,0);
% would load the R2.mat file from the current directory, and plot all
% electrodes/time frames F values thresholded using tfce at alpha 0.05
% topoplot and ERP like figures can't be automated since they require user
% input
%
% Cyril Pernet, Guillaume Rousselet, Carl Gaspar,
% Nicolas Chauveau, Andrew Stewart, Ramon Martinez-Cancino, Arnaud Delorme
%
% see also limo_stat_values limo_display_image limo_display_image_tf topoplot
% ----------------------------------------------------------------------
%  Copyright (C) LIMO Team 2024

if ~ischar(Type)
    options = { 'type', Type, 'filename', FileName, 'pathname', PathName, 'p', p, 'MCC', MCC, 'LIMO', LIMO, varargin{:} };
else
    options = {Type FileName PathName p MCC LIMO flag varargin{:} };
end

try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else
        g = [];
    end
catch
    limo_errordlg('limo_display_results() error: calling convention {''key'', value, ... } error');
    return
end

try g.channels;  catch, g.channels  = [];  end % No default values
try g.regressor; catch, g.regressor = [];  end % No default values
try g.plot3type; catch, g.plot3type = [];  end % No default values
try g.sumstats;  catch, g.sumstats  = [];  end % No default values
try g.restrict;  catch, g.restrict  = [];  end % No default values
try g.dimvalue;  catch, g.dimvalue  = [];  end % No default values
try g.fig;       catch, g.fig       = [];  end % Existing figure

if isequal(g.regressor, 0); g.regressor = []; end
if ~isempty(g.plot3type)
    extra = {'Original','Modelled','Adjusted'};
    if isnumeric(g.plot3type)
        extra = extra{g.plot3type};
    else
        extra(contains(extra,g.plot3type,'IgnoreCase',true));
    end
end
res = '';
toplot = load(fullfile(PathName,FileName));
toplot = toplot.(cell2mat(fieldnames(toplot)));

if nargin <= 6
    flag = 1;
end

choice = 'use theoretical p values'; % threshold based on what is computed since H0 is used for clustering
% see limo_stat_values - discontinuated empirical threshold (misleading)

% Load LIMO structure if a path was provided
% Load LIMO structure if a path was provided
if ischar(LIMO)
    load(LIMO, 'LIMO');
end

[~,FileNameTmp,ext] = fileparts(FileName);
if MCC == 2 || MCC == 4 % cluster and MAX correction
    LIMO.design.bootstrap = 1;

    % deal with bootstrap
    if ~exist([PathName filesep 'H0' filesep 'H0_' FileNameTmp ext],'file')
        if LIMO.Level == 1
            if strncmp(FileNameTmp,'con',3) || strncmp(FileNameTmp,'ess',3)
                limo_warndlg(sprintf('This contrast cannot be bootstrapped now, \nbootstrap the model and recompute the contrast'))
            else
                if strcmp(limo_questdlg('Level 1: are you sure to compute all bootstraps for that subject?','bootstrap turned on','Yes','No','No'),'Yes')
                    LIMO.design.bootstrap = 800;
                    if handles.tfce == 1
                        LIMO.design.tfce  = 1;
                    end
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                    limo_eeg(4);
                end
            end
        else % LIMO.Level == 2
            res = limo_questdlg('This option requires to compute bootstraps (this may take time)','Bootstraping data','Cancel','Continue','Continue');
            if ~strcmp(res,'Continue')
                return;
            end
            if ~isfield(LIMO.design, 'bootstrap') || LIMO.design.bootstrap == 1
                fprintf('Bootstrap repetition set to 1000')
                LIMO.design.bootstrap = 1000;
            end
            if contains(FileNameTmp,'one_sample')
                limo_random_robust(1,fullfile(LIMO.dir,'Yr.mat'),...
                    str2double(FileNameTmp(max(strfind(FileNameTmp,'_'))+1:end)),LIMO);
            elseif contains(FileNameTmp,'two_samples')
                limo_random_robust(2,fullfile(LIMO.dir,'Y1r.mat'),...
                    fullfile(LIMO.dir,'Y1r.mat'), str2double(FileNameTmp(max(strfind(FileNameTmp,'_'))+1:end)),LIMO);
            elseif contains(FileNameTmp,'paired_samples')
                underScoresPos = strfind(FileNameTmp,'_');
                param1         = str2num(FileNameTmp(underScoresPos(end-1)+1:underScoresPos(end)-1));
                param2         = str2num(FileNameTmp(underScoresPos(end)+1:end));
                limo_random_robust(3,fullfile(LIMO.dir,'Y1r.mat'),...
                    fullfile(LIMO.dir,'Y1r.mat'), [param1 param2],LIMO);
            elseif contains(FileNameTmp,'Covariate_effect') && contains(LIMO.design.name,'Regression')
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
                limo_eeg(4,LIMO.dir);
            elseif contains(FileNameTmp,'ANOVA') && ~strncmpi(FileNameTmp,'Rep_ANOVA',9)
                limo_random_robust(5,fullfile(LIMO.dir,'Yr.mat'), LIMO.data.Cat,LIMO.data.Cont,LIMO,'go','yes');
            elseif contains(FileNameTmp,'Rep_ANOVA')
                if strncmp(FileNameTmp,'con',3)
                    if exist([PathName filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'],'file')
                        limo_contrast([PathName filesep 'Yr.mat'], ...
                            [PathName filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'], LIMO, 0,3);
                    else
                        limo_errordlg('there is no GLM bootstrap file for this contrast file')
                    end
                elseif strncmp(FileNameTmp,'ess',3)
                    if exist([PathName filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'],'file')
                        limo_contrast([PathName filesep 'Yr.mat'], ...
                            [PathName filesep 'H0' filesep 'H0_' filesep 'H0_Betas.mat'], LIMO, 1,3);
                    else
                        limo_errordlg('there is no bootstrap file for this contrast file')
                    end
                else
                    disp('Bootstraping Repeated Measure ANOVA')
                    limo_random_robust(6,fullfile(PathName,'Yr.mat'),LIMO.data.Cat, ...
                        LIMO.design.repeated_measure, LIMO, 'go','yes')
                end
            end
        end
    end
elseif MCC == 3
    LIMO.design.tfce      = 1;
    currentfile = fullfile(PathName, FileName);
    if ~exist([PathName filesep 'H0' filesep 'tfce_H0_' FileNameTmp ext],'file')
        limo_tfce_handling(currentfile,'checkfile','yes')
    end
end

if ~isfield(LIMO,'Level')
    if Type == 3 % likely a summary stat file
        LIMO.Level = 2; % even if for a subject, calls limo_add_plots
    end
end

% -------------------------------------------------------------------------
% -------------------      LEVEL 1     ------------------------------------
% -------------------  SINGLE SUBJECT  ------------------------------------
% -------------------------------------------------------------------------
if LIMO.Level == 1

    switch Type

        case{1}

            %--------------------------
            % imagesc of the results
            %--------------------------

            if  strcmpi(LIMO.design.type_of_analysis,'Mass-univariate')

                % univariate results from 1st level analysis
                % ------------------------------------------

                % if previously plotted recover data from the cache
                data_cached = 0;
                if isfield(LIMO,'cache')
                    try
                        if strcmpi(LIMO.cache.fig.name, FileName) && ...
                                LIMO.cache.fig.MCC == MCC && ...
                                LIMO.cache.fig.threshold == p

                            disp('using cached data');
                            mask = LIMO.cache.fig.mask;
                            if isempty(mask)
                                data_cached = 0;
                            elseif sum(mask(:)) == 0
                                limo_errordlg('  no values under threshold  ','no significant effect','modal');
                                return
                            else
                                M           = LIMO.cache.fig.pval;
                                mytitle     = LIMO.cache.fig.title;
                                toplot      = LIMO.cache.fig.stats;
                                data_cached = 1;
                                assignin('base','p_values',M)
                                assignin('base','mask',mask)
                            end
                        end
                    catch no_cache
                        fprintf('could not load cached data %s',no_cache.message)
                        data_cached = 0;
                    end
                end

                % ------------------
                % compute the plot
                % ------------------
                if data_cached == 0

                    [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO);

                    if isempty(mask)
                        disp('no values computed'); return
                    elseif sum(mask(:)) == 0
                        limo_errordlg('  no values under threshold  ','no significant effect','modal');
                        LIMO.cache.fig.name       = FileName;
                        LIMO.cache.fig.MCC        = MCC;
                        LIMO.cache.fig.stats      = [];
                        LIMO.cache.fig.threshold  = p;
                        LIMO.cache.fig.pval       = M;
                        LIMO.cache.fig.mask       = mask;
                        LIMO.cache.fig.title      = mytitle;
                        % do an exception for designs with just the constant
                        if strcmpi(FileName,'R2.mat') && size(LIMO.design.X,2)==1
                            mask = ones(size(mask)); LIMO.cache.fig.mask = mask;
                            mytitle = 'R^2 Coef unthresholded'; LIMO.cache.fig.title = mytitle;
                            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
                        else
                            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
                            return
                        end
                    else
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                    end

                    if contains(FileName,'R2','IgnoreCase',true)
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,1)); % plot R2 values instead of F
                        else
                            toplot = squeeze(toplot(:,:,1));
                        end
                        assignin('base','R2_values',toplot)

                    elseif contains(FileName,'Condition_effect','IgnoreCase',true) || ...
                            contains(FileName,'Covariate_effect','IgnoreCase',true) || ...
                            contains(FileName,'Interaction_effect','IgnoreCase',true) || ...
                            contains(FileName,'semi_partial_coef','IgnoreCase',true)

                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,1)); % plot F values
                        else
                            toplot = squeeze(toplot(:,:,1));
                        end

                        if contains(FileName,'semi_partial_coef','IgnoreCase',true)
                            assignin('base','semi_partial_coef',toplot)
                        else
                            assignin('base','F_values',toplot)
                        end

                    elseif strcmpi(FileName(1:4),'con_')
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,4)); % plot T values
                        else
                            toplot = squeeze(toplot(:,:,4));
                        end
                        assignin('base','T_values',toplot)

                    elseif strcmpi(FileName(1:4),'ess_')
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,end-1)); % plot F values
                        else
                            toplot = squeeze(toplot(:,:,end-1));
                        end
                        assignin('base','F_values',toplot)

                    else
                        limo_errordlg('file not supported');
                        return
                    end
                end

                % replace plotting value with user regressor selection
                if ~isempty(g.regressor) && ~isequal(g.regressor, 0)
                    if ~exist('freq_index', 'var'), freq_index = []; end
                    toplot = limo_get_model_data(LIMO, g.regressor, extra, p, freq_index);
                end

                % -------------------------------------------------------------------------
                %              Actual plot takes place here
                % -------------------------------------------------------------------------
                if ~isempty(toplot)

                    % cache the results for next time
                    if data_cached == 0 && ~all(mask(:)==1)
                        LIMO.cache.fig.name       = FileName;
                        LIMO.cache.fig.MCC        = MCC;
                        LIMO.cache.fig.stats      = toplot;
                        LIMO.cache.fig.threshold  = p;
                        LIMO.cache.fig.pval       = M;
                        LIMO.cache.fig.mask       = mask;
                        LIMO.cache.fig.title      = mytitle;
                        if exist(LIMO.dir,"dir")
                            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
                        end
                    end

                    if ndims(toplot)==3
                        limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
                    else
                        limo_display_image(LIMO,toplot,mask,mytitle,flag);
                    end
                end

            else

                % mutivariate results from 1st level analysis
                % ------------------------------------------
                if strncmp(FileName,'R2',2) || strncmp(FileName,'Condition_effect',16) || strncmp(FileName,'Covariate_effect',16)   % MANOVA PLOTTING
                    if strncmp(FileName,'R2',2)

                        R2_EV     = load(fullfile(LIMO.dir,'R2_EV.mat'));
                        R2_EV     = R2_EV.R2_EV;
                        EV        = R2_EV(1:size(R2_EV,1),:); % no point plotting 0, just pick 5 1st Eigen values
                        R2_EV_var = load(fullfile(LIMO.dir,'R2_EV_var.mat'));
                        R2_EV_var = R2_EV_var.R2_EV_var;
                        test      =  sum(R2_EV_var(1,:) > 95) / size(R2_EV_var,2); % If more than 50% of the time-frames have a
                        % first eigenvalue with a proportion higher than 90%, the results of Roy's test are displayed,
                        if test > .50
                            choice = 'Roy';
                        else
                            choice = 'Pillai';
                        end
                        clear R2_EV;

                        F_values(:,1) = squeeze(toplot(:,2));
                        F_values(:,2) = squeeze(toplot(:,4));
                        [M, mask, mytitle] = limo_mstat_values(Type,FileName,p,MCC,LIMO,choice);
                        if isempty(mask)
                            return
                        elseif sum(mask(:)) == 0
                            limo_errordlg('  no values under threshold  ','no significant effect','modal');
                            return
                        else
                            toplot = squeeze(toplot(:,1)); % plot R2 values instead of F
                            assignin('base','F_values',F_values)
                            assignin('base','p_values',M)
                            assignin('base','mask',mask)
                            clear R2
                        end

                    else

                        if strcmpi(FileName(end-6:end),'_EV.mat')
                            FileName = [FileName(1:end-7) '.mat'];
                            toplot   = load(fullfile(PathName,FileName));
                            toplot   = toplot.(cell2mat(fieldnames(toplot)));
                        end
                        name   = sprintf('%s_%g_EV',FileName(1:end-4),str2double(FileName(max(strfind(FileName,'_')):end-4)));
                        EV     = load(fullfile(LIMO.dir,name));
                        EV     = EV.(cell2mat(fieldnames(EV)));
                        EV     = EV(1:size(Condition_effect_EV,1),:); % no point plotting 0, just pick 5 1st Eigen values
                        name   = sprintf('%s_%g_EV_var',FileName(1:end-4),str2double(FileName(max(strfind(FileName,'_')):end-4)));
                        EV_var = load(fullfile(LIMO.dir,name));
                        EV_var = EV_var.(cell2mat(fieldnames(EV_var)));
                        EV_var = EV_var(1:size(EV_var,1),:);
                        test =  sum(EV_var(1,:) > 95) / size(EV_var,2); % If more than 50% of the time-frames have a
                        %first eigenvalue with a proportion higher than 90%, the results of Roy's test are displayed,
                        if test > .50
                            choice = 'Roy';
                        else
                            choice = 'Pillai';
                        end

                        F_values(:,1) = squeeze(toplot(:,1));
                        F_values(:,2) = squeeze(toplot(:,3));
                        [M, mask, mytitle] = limo_mstat_values(Type,FileName,p,MCC,LIMO,choice);
                        if isempty(mask)
                            return
                        elseif sum(mask(:)) == 0
                            limo_errordlg('  no values under threshold  ','no significant effect','modal');
                            return
                        else
                            if strcmpi(choice,'Roy')
                                toplot = F_values(:,1);
                            else
                                toplot = F_values(:,2);
                            end
                            assignin('base','F_values',F_values)
                            assignin('base','p_values',M)
                            assignin('base','mask',mask)
                            clear R2
                        end
                    end

                    figure; set(gcf,'Color','w');
                    % imagesc eigen values
                    subplot(3,3,[4 5 7 8]);
                    timevect = linspace(LIMO.data.start,LIMO.data.end,size(EV,2));
                    scale = EV; scale(scale==0)=NaN;
                    imagesc(timevect,1:size(EV,1),scale);
                    color_images_(scale,LIMO);  colorbar
                    ylabel('Eigen Values','Fontsize',14)
                    set(gca,'YTickLabel',{'1','2','3','4','5'});
                    title('non-zero Eigen values','Fontsize',14)

                    % imagesc effect values
                    subplot(3,3,[1 2]);
                    scale = toplot'.*mask; scale(scale==0)=NaN;
                    imagesc(timevect,1,scale);
                    caxis([min(scale(:)), max(scale(:))]);
                    color_images_(scale,LIMO); xlabel(' ')
                    title(mytitle,'Fontsize',18); colorbar
                    ylabel(' '); set(gca,'YTickLabel',{''});

                    % ERP plot1 - Roy -
                    subplot(3,3,6);
                    plot(timevect, F_values(:,1),'LineWidth',3); grid on; axis tight
                    mytitle2 = sprintf('F values - Roy');
                    title(mytitle2,'FontSize',14)

                    % ERP plot2 - Pillai -
                    subplot(3,3,9);
                    plot(timevect, F_values(:,2),'LineWidth',3); grid on; axis tight
                    mytitle2 = sprintf('F value - Pillai');
                    title(mytitle2,'FontSize',14)
                end % end of MANOVA PLOTTING

                if strncmp(FileName,'Discriminant_coeff',18) || strncmp(FileName,'Discriminant_scores',19)
                    Discriminant_coeff      = load(fullfile(LIMO.dir,'Discriminant_coeff'));
                    Discriminant_coeff      = Discriminant_coeff.Discriminant_coeff;
                    Discriminant_scores     = load(fullfile(LIMO.dir,'Discriminant_scores'));
                    Discriminant_scores     = Discriminant_scores.Discriminant_scores;
                    Condition_effect_EV_var = load(fullfile(LIMO.dir,'Condition_effect_1_EV_var.mat'));
                    Condition_effect_EV_var = Condition_effect_EV_var.Condition_effect_EV_var;

                    time = linspace(LIMO.data.start,LIMO.data.end, size(Discriminant_coeff,2));
                    input_title = sprintf('which time-frame to plot (in ms)?: ');
                    timepoint = inputdlg(input_title,'Plotting option');
                    t = dsearchn(time', str2double(timepoint{1}));
                    groupcolors = 'rgbcwmryk';
                    groupsymbols = 'xo*+.sdv<>';
                    [class,~] = find(LIMO.design.X(:,1:LIMO.design.nb_conditions)');
                    k = LIMO.design.nb_conditions;

                    if k>2
                        figure;set(gcf,'Color','w');
                        subplot(2,2,[1 2]); % 2D plot of two discriminant functions
                        gscatter(squeeze(Discriminant_scores(1,t,:)), squeeze(Discriminant_scores(2,t,:)), class, groupcolors(1:k), groupsymbols(1:k));
                        grid on; axis tight;
                        xlabel(['Z1, var: ' num2str(round(Condition_effect_EV_var(1,t)),2) '%'],'Fontsize',14);
                        ylabel(['Z2, var: ' num2str(round(Condition_effect_EV_var(2,t)),2) '%'],'Fontsize',14);
                        title(['Results of the discriminant analysis at ' num2str(time(t)) 'ms'], 'Fontsize', 18);
                        z1 = subplot(2,2,3); % First discriminant coeff
                        cc = limo_color_images(Discriminant_coeff(:,t,1)); % get a color map commensurate to that
                        topoplot(Discriminant_coeff(:,t,1),LIMO.data.chanlocs, 'electrodes','off','style','map','whitebk', 'on','colormap',cc);colorbar;
                        title('Z1','Fontsize',14); colormap(z1, 'hot');
                        z2 = subplot(2,2,4); % Second discriminant coeff
                        cc = limo_color_images(Discriminant_coeff(:,t,2)); % get a color map commensurate to that
                        topoplot(Discriminant_coeff(:,t,2),LIMO.data.chanlocs, 'electrodes','off','style','map','whitebk', 'on','colormap',cc);colorbar;
                        title('Z2','Fontsize',14); colormap(z2, 'hot');
                    elseif k==2
                        figure;set(gcf,'Color','w');
                        subplot(2,2,[1 2]); % 1D plot of two discriminant functions
                        data = squeeze(Discriminant_scores(1,t,:));
                        class1 = data(class == 1);
                        class2 = data(class == 2);
                        histogram(class1, 'BinWidth', 0.1);
                        hold on
                        histogram(class2,'BinWidth',0.1);
                        hold off
                        legend show
                        grid on; axis tight;
                        xlabel(['Z1, var: ' num2str(round(Condition_effect_EV_var(1,t)),2) '%'],'Fontsize',14);
                        title(['Results of the discriminant analysis at ' num2str(time(t)) 'ms'], 'Fontsize', 18);
                        z1 = subplot(2,2,[3,4]); % First discriminant coeff
                        cc = limo_color_images(Discriminant_coeff(:,t,1)); % get a color map commensurate to that
                        topoplot(Discriminant_coeff(:,t,1),LIMO.data.chanlocs, 'electrodes','off','style','map','whitebk', 'on','colormap',cc);colorbar;
                        title('Z1','Fontsize',14); colormap(z1, 'hot');
                    end
                    limo_display_image(LIMO,abs(Discriminant_coeff(:,:,1)),abs(Discriminant_coeff(:,:,1)),'Discriminant coefficients Z1',flag)

                    %                     figure;set(gcf,'Color','w');
                    %                     for t=1:size(Discriminant_coeff,2)
                    %                     topoplot(Discriminant_coeff(:,t,1),LIMO.data.chanlocs, 'electrodes','numbers','style','map');
                    %                     title(['Discriminant values first discriminant at timepoint ' num2str(t) ' corresponding to ' num2str(time(t)) ' ms']);
                    %                     pause(.01)
                    %                     end;
                end

                if strncmp(FileName,'Linear_Classification',21)
                    Linear_Classification = load(fullfile(LIMO.dir,'Linear_Classification'));
                    Linear_Classification = Linear_Classification.Linear_Classification;
                    [~, mask, mytitle] = limo_mstat_values(Type,FileName,p,MCC,LIMO,choice);
                    timevect = linspace(LIMO.data.start,LIMO.data.end,size(Linear_Classification,1));
                    figure;set(gcf,'Color','w');
                    subplot(3,1,[1 2]); % lineplot
                    plot(timevect,Linear_Classification(:,2),'LineWidth',3);title(mytitle, 'Fontsize', 18);
                    ylabel('decoding accuracies', 'Fontsize', 14);grid on; axis tight; hold on;
                    plot(timevect, Linear_Classification(:,2) + 2*Linear_Classification(:,3), 'k-','LineWidth',1); hold on;
                    plot(timevect, Linear_Classification(:,2) - 2*Linear_Classification(:,3), 'k-','LineWidth',1)
                    line([0,0],[0,1], 'color', 'black')
                    subplot(3, 1, 3); % imagesc accuracies
                    toplot = Linear_Classification(:,2); scale = toplot'.*mask;scale(scale==0)=NaN;
                    imagesc(timevect,1,scale);xlabel('Time in ms');
                    color_images_(scale,LIMO);
                    ylabel(' '); set(gca,'YTickLabel',{''});
                end

                if strncmp(FileName,'Quadratic_Classification',24)
                    Quadratic_Classification = load(fullfile(LIMO.dir,'Quadratic_Classification'));
                    Quadratic_Classification = Quadratic_Classification.Quadratic_Classification;
                    timevect = linspace(LIMO.data.start,LIMO.data.end,size(Quadratic_Classification,1));
                    figure;set(gcf,'Color','w');
                    subplot(3,1,[1 2]); % lineplot
                    plot(timevect,Quadratic_Classification(:,2),'LineWidth',3);title('CV quadratic decoding accuracies +/- 2SD', 'Fontsize', 18);
                    ylabel('decoding accuracies', 'Fontsize', 14);grid on; axis tight; hold on
                    plot(timevect, Quadratic_Classification(:,2) + 2*Quadratic_Classification(:,3), 'k-','LineWidth',1); hold on;
                    plot(timevect, Quadratic_Classification(:,2) - 2*Quadratic_Classification(:,3), 'k-','LineWidth',1)
                    line([0,0],[0,1], 'color', 'black')
                    subplot(3,1, 3); % imagesc plot training accuracies
                    scale = Quadratic_Classification(:,2)'; scale(scale==0)=NaN;
                    imagesc(timevect,1,scale);xlabel('Time in ms');
                    color_images_(scale,LIMO);
                    ylabel(' '); set(gca,'YTickLabel',{''});
                end
            end


        case{2}

            %--------------------------
            % topoplot
            %--------------------------

            % univariate results from 1st level analysis
            % ------------------------------------------

            if strcmpi(LIMO.Analysis,'Time-Frequency')
                warndlg('topoplot not supported for 3D data')
            else
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                if strcmpi(LIMO.Analysis,'Time')
                    EEG.xmin  = LIMO.data.start / 1000;% in msec
                    EEG.xmax  = LIMO.data.end / 1000;  % in msec
                    EEG.times = LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000; % in sec;
                    if length(EEG.times) > 2
                        EEG.times = [EEG.times(1) EEG.times(end)];
                    end
                elseif strcmpi(LIMO.Analysis,'Frequency')
                    EEG.xmin = LIMO.data.freqlist(1);
                    EEG.xmax = LIMO.data.freqlist(end);
                    freqlist = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
                    if isempty(freqlist)
                        return
                    else
                        if contains(cell2mat(freqlist),':')
                            EEG.freq = eval(cell2mat(freqlist));
                        else
                            EEG.freq = str2double(cell2mat(freqlist));
                        end

                        if min(EEG.freq)<EEG.xmin || max(EEG.freq)>EEG.xmax
                            limo_errordlg('selected frequency out of bound'); return
                        end
                    end
                end

                if contains(FileName,'R2','IgnoreCase',true)
                    if size(LIMO.design.X,2)==1
                        EEG.data    = squeeze(toplot(:,:,1));
                        EEG.setname = 'R2 values for the mean';
                    else
                        EEG.data    = squeeze(toplot(:,:,2));
                        EEG.setname = 'R2 - F values';
                    end
                    call_topolot(EEG,FileName,LIMO.Analysis)
                elseif contains(FileName,'Condition_effect','IgnoreCase',true) || ...
                        contains(FileName,'Covariate_effect','IgnoreCase',true) || ...
                        contains(FileName,'Interaction_effect','IgnoreCase',true)  || ...
                        contains(FileName,'Condition_effect','IgnoreCase',true)
                    EEG.data    = squeeze(toplot(:,:,1));
                    call_topolot(EEG,FileName,LIMO.Analysis)
                elseif contains(FileName,'con','IgnoreCase',true) || contains(FileName,'ess','IgnoreCase',true)
                    EEG.data    = squeeze(toplot(:,:,end-1));
                    call_topolot(EEG,FileName,LIMO.Analysis)
                elseif contains(FileName,'semi partial_coef.mat','IgnoreCase',true)
                    regressor = str2double(cell2mat(inputdlg('which regressor(s) to plot (e.g. 1:3)','Plotting option')));
                    if max(regressor) > size(toplot,3); errordlg('error in regressor number'); return; end
                    for b = regressor
                        EEG.data     = squeeze(toplot(:,:,b,1));
                        call_topolot(EEG,FileName,LIMO.Analysis)
                    end
                else
                    disp('file not supported');
                    return
                end

                if contains(FileName,'con','IgnoreCase',true)
                    assignin('base','T_values',EEG.data);
                else
                    assignin('base','F_values',EEG.data);
                end
            end

        case{3}

            %--------------------------
            % Time course / Power
            %--------------------------

            % which variable(s) to plot
            % ----------------------
            if isempty(g.regressor)
                input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2));
                regressor   = inputdlg(input_title,'Plotting option');
            else
                regressor = g.regressor;
            end

            if isempty(regressor); disp('selection aborded'); return; end
            regressor = cell2mat(regressor);
            if isempty(regressor); disp('selection aborded'); return; end
            if ~contains(regressor,'['); regressor=['[' regressor ']']; end
            if ischar(regressor); regressor=str2num(regressor); end %#ok<ST2NM>
            regressor = sort(regressor);

            if max(regressor) > size(LIMO.design.X,2)
                limo_errordlg('invalid regressor number');
            end

            categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
            if max(regressor) == size(LIMO.design.X,2)
                tmp = regressor(1:end-1);
            else
                tmp = regressor;
            end

            cat = sum(tmp<=categorical); cont = sum(tmp>categorical);
            if cat >=1 && cont >=1
                errordlg2('you can''t plot categorical and continuous regressors together'); return
            end

            % which data type to make
            % ------------------------
            if isempty(g.plot3type) && ~any(strcmpi(g.plot3type,{'Original','Modelled','Adjusted'}))
                extra = questdlg('Which data type to plot?','Options','Original','Modelled','Adjusted','Adjusted');
            else
                extra = g.plot3type;
            end
            if isempty(extra)
                return
            elseif strcmpi(extra,'Original')
                if regressor == size(LIMO.design.X,2)
                    limo_errordlg('you can''t plot adjusted mean for original data'); return
                end
            end

            % timing /frequency info
            % -----------------------
            if strcmpi(LIMO.Analysis,'Time')
                timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
            elseif strcmpi(LIMO.Analysis,'Frequency')
                freqvect=LIMO.data.freqlist';
            elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                timevect = linspace(LIMO.data.start,LIMO.data.end,LIMO.data.size4D(3));
                freqvect = linspace(LIMO.data.lowf,LIMO.data.highf,LIMO.data.size4D(2));
            end

            % which channel/frequency to plot
            % --------------------------------
            if isempty(g.channels)
                channel = inputdlg('which channel to plot','Plotting option');
            else
                channel = g.channels;
            end
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                disp('loading the 4D data ...')
                frequency = inputdlg('which Frequency to plot','Plotting option');
            else
                frequency = [];
            end

            if strcmpi(channel,'') || strcmpi(frequency,'')
                disp('looking for max');
                R2 = load(fullfile(LIMO.dir,'R2.mat'));
                R2 = R2.R2;
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    tmp = squeeze(R2(:,:,:,1)); clear R2
                    [e,f,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                    if length(e) ~= 1; e = e(1); f = f(1); end

                    if strcmpi(channel,'')
                        channel = e;
                    else
                        channel = eval(cell2mat(channel));
                    end
                    if size(channel) > 1
                        errordlg('invalid channel choice'); return
                    elseif channel > size(LIMO.data.chanlocs,2) || channel < 1
                        errordlg('invalid channel number'); return
                    end

                    if strcmpi(frequency,'')
                        freq_index = f;
                        frequency = freqvect(freq_index);
                    else
                        frequency = eval(cell2mat(frequency));
                    end
                    if size(frequency) > 1
                        errordlg('invalid frequency choice'); return
                    elseif frequency > LIMO.data.tf_freqs(end) || frequency < LIMO.data.tf_freqs(1)
                        errordlg('invalid frequency number'); return
                    end
                else
                    tmp = squeeze(R2(:,:,1)); clear R2
                    [channel,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                end
                clear tmp
            else
                channel = eval(cell2mat(channel));
                if size(channel) > 1
                    limo_errordlg('invalid channel choice'); return
                elseif channel > size(LIMO.data.chanlocs,2) || channel < 1
                    limo_errordlg('invalid channel number'); return
                end

                if ~isempty(frequency)
                    frequency = eval(cell2mat(frequency));
                    if size(frequency) > 1
                        errordlg('invalid frequency choice');
                    elseif frequency > freqvect(end) || frequency < freqvect(1)
                        errordlg('invalid frequency number');
                    end
                    % pick the nearest frequency index
                    [~, freq_index] = min(abs(freqvect-frequency ));
                    frequency = freqvect(freq_index);
                end
            end

            % down to business
            % ----------------------
            data_cached = 0;
            if isfield(LIMO,'cache')
                if strcmpi(LIMO.Analysis,'Time-Frequency') && isfield(LIMO.cache,'ERPplot')

                    if mean([LIMO.cache.Courseplot.channel == channel ...
                            LIMO.cache.Courseplot.regressor == regressor ...
                            LIMO.cache.Courseplot.frequency == frequency]) == 1 ...
                            && strcmpi('LIMO.cache.Courseplot.extra',extra)

                        if sum(regressor <= categorical) == length(regressor)
                            average = LIMO.cache.Courseplot.average;
                            ci = LIMO.cache.Courseplot.ci;
                            mytitle = LIMO.cache.Courseplot.title;
                            disp('using cached data');
                            data_cached = 1;
                        else
                            continuous = LIMO.cache.Courseplot.continuous;
                            mytitle = LIMO.cache.Courseplot.title;
                            disp('using cached data');
                            data_cached = 1;
                        end
                    end

                elseif strcmpi(LIMO.Analysis,'Time') && isfield(LIMO.cache,'Courseplot') || ...
                        strcmpi(LIMO.Analysis,'Frequency') && isfield(LIMO.cache,'Courseplot')

                    if length(LIMO.cache.Courseplot.regressor) == length(channel)
                        if mean([LIMO.cache.Courseplot.channel == channel ...
                                LIMO.cache.Courseplot.regressor == regressor]) == 1  ...
                                && strcmpi('LIMO.cache.Courseplot.extra',extra)

                            if sum(regressor <= categorical) == length(regressor)
                                average = LIMO.cache.Courseplot.average;
                                ci = LIMO.cache.Courseplot.ci;
                                mytitle = LIMO.cache.Courseplot.title;
                                disp('using cached data');
                                data_cached = 1;
                            else
                                continuous = LIMO.cache.Courseplot.continuous;
                                mytitle = LIMO.cache.Courseplot.title;
                                disp('using cached data');
                                data_cached = 1;
                            end
                        end
                    end
                end
            end

            % no cache = compute
            if data_cached == 0

                probs = [p/2; 1-p/2];
                z = norminv(probs);

                if strcmpi(extra,'Original')
                    Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                data     = squeeze(Yr(channel,freq_index,:,index{i}));
                                mytitle  = sprintf('Original ERSP at \n channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                data     = squeeze(Yr(channel,:,index{i}));
                            end
                            average(i,:) = nanmean(data,2);
                            se           = nanstd(data,0,2) ./ sqrt(numel(index{i}));
                            ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Yr,2));
                        end

                        if strcmpi(LIMO.Analysis,'Time')
                            mytitle = sprintf('Original ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmpi(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Original Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                continuous(i,:,:) = Yr(channel,freq_index,:,sorting_values);
                                mytitle{i}        = sprintf('Original single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            else
                                continuous(i,:,:) = Yr(channel,:,sorting_values);
                                mytitle{i}        = sprintf('Original single trials \n sorted by regressor %g \n channel %s (%g) at %s Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            end
                        end
                    end
                    clear Yr
                elseif strcmpi(extra,'Modelled')
                    Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                    Betas = Betas.Betas;
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        Betas = squeeze(Betas(channel,freq_index,:,:));
                    else
                        Betas = squeeze(Betas(channel,:,:));
                    end
                    Yh = (LIMO.design.X*Betas')'; % modelled data

                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            Yr = squeeze(Yr(:,freq_index,:,:));
                        end
                        R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));

                        for i=length(regressor):-1:1
                            index{i}     = find(LIMO.design.X(:,regressor(i)));
                            data         = squeeze(Yh(:,index{i}));
                            average(i,:) = mean(data,2);
                            var          = diag(((R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')) / LIMO.model.model_df(2));
                            CI           = sqrt(var/size(index{i},1))*z';
                            ci(i,:,:)    = (repmat(mean(data,2),1,2)+CI)';
                        end

                        if strcmpi(LIMO.Analysis,'Time')
                            mytitle = sprintf('Modelled ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmpi(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Modelled Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Modelled ERSP \n channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i}                         = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values] = sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:)                = Yh(:,sorting_values);

                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g \n channel %s (%g) at %g Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            end
                        end
                    end
                else % Adjusted
                    allvar = 1:size(LIMO.design.X,2)-1;
                    allvar(regressor)=[];
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
                        Yr    = squeeze(Yr.Yr(channel,freq_index,:,:));
                        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                        Betas = squeeze(Betas.Betas(channel,freq_index,:,:));
                    else
                        Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
                        Yr    = squeeze(Yr.Yr(channel,:,:));
                        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                        Betas = squeeze(Betas.Betas(channel,:,:));
                    end
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                    Ya        = Yr - confounds; clear Yr Betas confounds;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=length(regressor):-1:1
                            index{i}     = find(LIMO.design.X(:,regressor(i)));
                            data         = squeeze(Ya(:,index{i}));
                            average(i,:) = nanmean(data,2);
                            se           = nanstd(data,0,2) ./ sqrt(numel(index{i}));
                            ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Ya,1));
                        end
                        if strcmpi(LIMO.Analysis,'Time')
                            mytitle = sprintf('Adjusted ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmpi(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Adjusted Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Adjusted ERSP channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i}                         = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values] = sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:)                = Ya(:,sorting_values);
                            if strcmpi(LIMO.Analysis,'Time-Frequency')
                                mytitle{i} = sprintf('Adjusted single trials \n sorted by regressor \n %g channel %s (%g) at %g Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                mytitle{i} = sprintf('Adjusted single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            end
                        end
                    end
                end
            end

            % make the figure(s)
            % ------------------
            figure;set(gcf,'Color','w')
            if sum(regressor <= categorical) == length(regressor)
                for i=1:size(average,1)
                    if i==1
                        colorOrder = get(gca, 'ColorOrder');
                        colorOrder = repmat(colorOrder,ceil(size(average,1)/size(colorOrder,1)),1);
                    end

                    if strcmpi(LIMO.Analysis,'Frequency')
                        try
                            plot(freqvect,average(i,:),'LineWidth',1.5,'Color',colorOrder(i,:)); hold on
                        catch
                            freqvect = linspace(LIMO.data.start,LIMO.data.end,size(average,2));
                            plot(freqvect,average(i,:),'LineWidth',1.5,'Color',colorOrder(i,:)); hold on
                        end
                    else
                        plot(timevect,average(i,:),'LineWidth',1.5,'Color',colorOrder(i,:)); hold on
                    end

                    x = squeeze(ci(i,1,:)); y = squeeze(ci(i,2,:));
                    if strcmpi(LIMO.Analysis,'Frequency')
                        fillhandle = patch([reshape(freqvect, 1, numel(freqvect)) fliplr(reshape(freqvect, 1, numel(freqvect)))], [x' fliplr(y')], colorOrder(i,:));
                    else
                        fillhandle = patch([reshape(timevect, 1, numel(timevect)) fliplr(reshape(timevect, 1, numel(timevect)))], [x',fliplr(y')], colorOrder(i,:));
                    end
                    set(fillhandle,'EdgeColor',colorOrder(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);
                end

                % if regressor spans columns of an effect, plot significant time frames
                index = 1; index2 = LIMO.design.nb_conditions(1);
                for i=1:length(LIMO.design.nb_conditions)
                    effect = index:index2;
                    if length(regressor) == length(effect)
                        if mean(regressor == effect) == 1
                            name = sprintf('Condition_effect_%g.mat',i);
                            % load(name);
                            if isfield(LIMO,'cache') && isfield(LIMO,'fig')
                                if strcmpi(LIMO.cache.fig.name,name) && ...
                                        LIMO.cache.fig.MCC == MCC && ...
                                        LIMO.cache.fig.threshold == p
                                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                                        sig = single(LIMO.cache.fig.mask(channel,freq_index,:)); sig(sig==0)=NaN;
                                    else
                                        sig = single(LIMO.cache.fig.mask(channel,:)); sig(sig==0)=NaN;
                                    end
                                end
                            else
                                if strcmpi(LIMO.Analysis,'Time-Frequency')
                                    [~, mask] = limo_stat_values(name,p,MCC,LIMO);
                                    sig = single(squeeze(mask(channel,freq_index,:))); sig(sig==0)=NaN;
                                else
                                    [~, mask] = limo_stat_values(name,p,MCC,LIMO);
                                    sig = single(mask(channel,:)); sig(sig==0)=NaN;
                                end
                            end
                            h = axis;
                            if strcmpi(LIMO.Analysis,'Frequency')
                                plot(freqvect,sig.*h(3),'r*','LineWidth',2)
                            else
                                plot(timevect,sig.*h(3),'r*','LineWidth',2)
                            end
                            break
                        end
                    else
                        index = index+LIMO.design.nb_conditions(i);
                        if i<length(LIMO.design.nb_conditions)
                            index2 = LIMO.design.nb_conditions(i)+LIMO.design.nb_conditions(i+1);
                        end
                    end
                end

                if LIMO.design.nb_interactions ~= 0
                    index = sum(LIMO.design.nb_conditions)+1; index2 = sum(LIMO.design.nb_conditions)+LIMO.design.nb_interactions(1);
                    for i=1:length(LIMO.design.nb_interactions)
                        effect = index:index2;
                        if length(regressor) == length(effect)
                            if mean(regressor == effect) == 1
                                name = sprintf('Interaction_effect_%g.mat',i);
                                % load(name);
                                if isfield(LIMO,'cache') && isfield(LIMO,'fig')
                                    if strcmpi(LIMO.cache.fig.name,name) && ...
                                            LIMO.cache.fig.MCC == MCC && ...
                                            LIMO.cache.fig.threshold == p
                                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                                            sig = single(squeeze(LIMO.cache.fig.mask(channel,freq_index,:))); sig(sig==0)=NaN;
                                        else
                                            sig = single(LIMO.cache.fig.mask(channel,:)); sig(sig==0)=NaN;
                                        end
                                    end
                                else
                                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                                        [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                        sig = single(squeeze(mask(channel,freq_index,:))); sig(sig==0)=NaN;
                                    else
                                        [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                        sig = single(mask(channel,:)); sig(sig==0)=NaN;
                                    end
                                end
                                h = axis;
                                if strcmpi(LIMO.Analysis,'Frequency')
                                    plot(freqvect,sig.*h(3),'r*','LineWidth',2)
                                else
                                    plot(timevect,sig.*h(3),'r*','LineWidth',2)
                                end
                                break
                            end
                        else
                            index = index+LIMO.design.nb_interactions(i);
                            if i<length(LIMO.design.nb_interactions)
                                index2 = LIMO.design.nb_interactions(i)+LIMO.design.nb_interactions(i+1);
                            end
                        end
                    end
                end

                % --
                axis tight; grid on; box on
                title(mytitle,'FontSize',19); drawnow;
                assignin('base','Plotted_data', average)
                set(gca,'FontSize',14,'Layer','Top')
                if strcmpi(LIMO.Analysis,'Frequency')
                    xlabel('Freq in Hz','FontSize',16)
                    ylabel('Power Spectrum in {\mu}V^2/Hz','FontSize',16);
                else
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude in {\mu}V','FontSize',16)
                end

                LIMO.cache.Courseplot.extra     = extra;
                LIMO.cache.Courseplot.average   = average;
                LIMO.cache.Courseplot.channel   = channel;
                LIMO.cache.Courseplot.regressor = regressor;
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.Courseplot.frequency = frequency;
                end
                LIMO.cache.Courseplot.ci        = ci;
                LIMO.cache.Courseplot.title     = mytitle;
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')

            else
                for i=1:size(continuous,1)
                    if i > 1; figure;set(gcf,'Color','w'); end
                    index = find(~isnan(squeeze(continuous(i,1,:))));

                    if strcmpi(LIMO.Analysis,'Frequency')
                        try
                            surf(index,freqvect,squeeze(continuous(i,:,index)));shading interp
                        catch
                            freqvect = linspace(LIMO.data.start,LIMO.data.end,size(continuous,2));
                            surf(index,freqvect,squeeze(continuous(i,:,index)));shading interp
                        end
                        ylabel('Frequency in Hz','FontSize',16)
                        zlabel('Power Spectrum in {\mu}V^2/Hz','FontSize',16)
                    else
                        surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Time in ms','FontSize',16)
                        zlabel('Amplitude in {\mu}V','FontSize',16)
                    end
                    % --
                    axis tight; title(mytitle{i},'FontSize',14); drawnow;
                    xlabel('Sorted trials','FontSize',16)
                    try %#ok<TRYNC>
                        set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                    end
                end

                LIMO.cache.Courseplot.continuous = continuous;
                LIMO.cache.Courseplot.channel    = channel;
                LIMO.cache.Courseplot.regressor  = regressor;
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.Courseplot.frequency = frequency;
                end
                LIMO.cache.Courseplot.title      = mytitle;
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
            end

    end % closes switch

    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % -------------------                               LEVEL 2
    % -------------------                            GROUP EFFECTS
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------

elseif LIMO.Level == 2

    if contains(FileName,'LIMO') || contains(FileName,'Y')
        limo_warndlg('select a statistical result file - plot aborded')
        return
    end

    % if previously plotted, recover from the cache
    data_cached = 0;
    if isfield(LIMO,'cache')
        try
            if strcmpi(LIMO.cache.fig.name, FileName) && ...
                    LIMO.cache.fig.MCC == MCC && ...
                    LIMO.cache.fig.threshold == p

                disp('using cached data');
                mask = LIMO.cache.fig.mask;
                if isempty(mask)
                    data_cached = 0;
                elseif sum(mask(:)) == 0
                    limo_errordlg('  no values under threshold  ','no significant effect','modal');
                    return
                else
                    toplot      = LIMO.cache.fig.stats;
                    M           = LIMO.cache.fig.pval;
                    mask        = LIMO.cache.fig.mask;
                    mytitle     = LIMO.cache.fig.title;
                    data_cached = 1;
                    assignin('base','stat_values',toplot)
                    assignin('base','p_values',M)
                    assignin('base','mask',mask)
                end
            end
        catch no_cache
            data_cached = 0;
            limo_errordlg(no_cache,'failed to chache data %s',no_cache.message)
        end
    end

    % if there is no cached data, compute and plot
    % -------------------------------------------
    skip_stat_values = false;
    if Type == 3 % quickly check if central tendency file
        tmp = load(FileName);
        if isfield(tmp,'Data')
            skip_stat_values = true;
        end
    end

    if data_cached == 0 && ~skip_stat_values

        [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO);

        if isempty(mask) || sum(mask(:)) == 0
            p = 1; mask = zeros(size(M));
            assignin('base','p_values',squeeze(M))
            limo_warndlg('  no values under threshold  ','no significant effect','modal');
        else
            assignin('base','p_values',squeeze(M))
            assignin('base','mask',squeeze(mask))
        end

        if strcmpi(LIMO.Analysis,'Time-Frequency') || strcmpi(LIMO.Analysis,'ITC')
            if contains(FileName,'R2') || ...
                    contains(FileName,'semi_partial')
                toplot = squeeze(toplot(:,:,:,1));
            elseif contains(FileName,'ttest','IgnoreCase',true) || ...
                    contains(FileName,'LI_Map','IgnoreCase',true)
                toplot = squeeze(toplot(:,:,:,4));
            elseif strncmp(FileName,'con_',4)
                toplot = squeeze(toplot(:,:,:,4));
            elseif strncmp(FileName,'ess_',4)
                if ~exist('ess','var')
                    effect_nb = eval(FileName(5:end-4)); %#ok<NASGU>
                end
                toplot = squeeze(toplot(:,:,:,end-1));
            elseif contains(FileName,'Condition') || ...
                    contains(FileName,'Covariate') || ...
                    contains(FileName,'Rep_ANOVA')
                toplot = squeeze(toplot(:,:,:,1));
            else
                disp('file no supported'); return
            end
        else
            if contains(FileName,'R2') || ...
                    contains(FileName,'semi_partial')
                toplot = squeeze(toplot(:,:,1));
            elseif contains(FileName,'ttest','IgnoreCase',true) || ...
                    contains(FileName,'LI_Map','IgnoreCase',true)
                toplot = squeeze(toplot(:,:,4));
            elseif strncmp(FileName,'con_',4)
                toplot = squeeze(toplot(:,:,4));
            elseif strncmp(FileName,'ess_',4)
                if ~exist('ess','var')
                    effect_nb = eval(FileName(5:end-4)); %#ok<NASGU>
                end
                toplot = squeeze(toplot(:,:,end-1));
            elseif contains(FileName,'Condition') || ...
                    contains(FileName,'Covariate') || ...
                    contains(FileName,'Rep_ANOVA')
                toplot = squeeze(toplot(:,:,1));
            else
                disp('file no supported'); return
            end
        end
        assignin('base','stat_values',toplot)
        data_cached = 0;
    end

    % ------------------------------
    %      Image and topoplot
    % ----------------------------
    if Type == 1 || Type == 2

        % cache the results for next time
        % ------------------------------
        if data_cached == 0
            LIMO.cache.fig.name       = FileName;
            LIMO.cache.fig.MCC        = MCC;
            LIMO.cache.fig.stats      = toplot;
            LIMO.cache.fig.threshold  = p;
            LIMO.cache.fig.pval       = squeeze(M);
            LIMO.cache.fig.mask       = squeeze(mask);
            LIMO.cache.fig.title      = mytitle;
            if exist(LIMO.dir,'dir')
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
            else
                disp('Cached data in LIMO.mat cannot be updated - LIMO dir doesn''t exist (likely moved files)')
            end
        end

        % replace toplot by user selection
        if ~isempty(g.regressor) && ~isequal(g.regressor, 0)
            if ~exist('freq_index', 'var'), freq_index = []; end
            toplot = limo_get_model_data(LIMO, g.regressor, extra, p, freq_index);
        end

        % image all results
        % ------------------
        if Type == 1 && ~strcmpi(LIMO.Analysis,'Time-Frequency') && ~strcmpi(LIMO.Analysis,'ITC')
            limo_display_image(LIMO,toplot,mask,mytitle,flag)

        elseif Type == 1 && strcmpi(LIMO.Analysis,'Time-Frequency') || ...
                Type == 1 && strcmpi(LIMO.Analysis,'ITC')
            if ndims(toplot)==3
                limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
            else
                limo_display_image(LIMO,squeeze(toplot),squeeze(mask),mytitle,flag)
            end

        elseif Type == 2
            %--------------------------
            % topoplot
            %--------------------------

            if strcmpi(LIMO.Analysis,'Time-Frequency')
                errordlg('topoplot not supported for time-frequency analyses')
            else
                if isfield(LIMO.design,'channel')  % not full scalp
                    if ~isempty(LIMO.design.electrode)
                        msgbox('Only one channel found','No topoplot')
                        return
                    end
                end
            end

            if sum(mask(:)) == 0
                limo_errordlg('no values under threshold','no significant effect');
            else
                EEG.data     = toplot;
                EEG.setname  = mytitle;
                EEG.chanlocs = LIMO.data.chanlocs;

                if size(toplot,2) == 1
                    opt = {'maplimits','maxmin','verbose','off'};
                    if isfield(LIMO,'Type')
                        if strcmpi(LIMO.Type,'Components')
                            opt = {'maplimits','absmax','electrodes','off','verbose','off'};
                        end
                    end
                    figure; set(gcf,'Color','w','InvertHardCopy','off');
                    topoplot(toplot(:,1),EEG.chanlocs,opt{:});
                    title('Topoplot','FontSize',12)
                else
                    if strcmpi(LIMO.Analysis,'Time')
                        EEG.xmin  = LIMO.data.start/1000; % in sec
                        EEG.xmax  = LIMO.data.end/1000;   % in sec
                        EEG.times = (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
                        call_topolot(EEG,FileName,LIMO.Analysis)
                        % pop_topoplot(EEG);
                    elseif strcmpi(LIMO.Analysis,'Frequency')
                        EEG.xmin  = LIMO.data.freqlist(1);
                        EEG.xmax  = LIMO.data.freqlist(end);
                        freqlist  = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
                        if isempty(freqlist)
                            return
                        else
                            EEG.freq = str2double(cell2mat(freqlist));
                            if min(EEG.freq)<EEG.xmin || max(EEG.freq)>EEG.xmax
                                errordlg('selected frequency out of bound'); return
                            end
                        end
                        call_topolot(EEG,FileName,LIMO.Analysis)
                    end
                    assignin('base','Plotted_data',EEG.data)
                end
            end
        end

    elseif Type == 3

        %--------------------------
        % Course plot
        %--------------------------

        if contains(FileName,'one_sample','IgnoreCase',true) || contains(FileName,'two_samples','IgnoreCase',true) || ...
                contains(FileName,'paired_samples','IgnoreCase',true) || contains(FileName,'con_','IgnoreCase',true) || ...
                contains(FileName,'ess_','IgnoreCase',true)
            % ------------------------------------------------------------------------------------------------------------
            % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
            % H0 file dim = (electrodes,frames,[t, p],nboot)

            data = load(fullfile(PathName,FileName));
            data = data.(cell2mat(fieldnames(data)));
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [~,channel,freq,time] = limo_display_reducedim(data(:,:,:,[4 5]),LIMO,g.channels,g.restrict,g.dimvalue);
                data                  = squeeze(data(channel,freq,time,:,:)); % 2D
                sig                   = squeeze(single(mask(channel,freq,time))); %1D
            else
                [~,channel,freq,time] = limo_display_reducedim(data(:,:,[4 5]),LIMO,g.channels);
                data                  = squeeze(data(channel,time,:));
                sig                   = single(mask(channel,:));
            end
            sig(sig==0)=NaN;

            % compute
            trimci      = NaN(size(data,1),3);
            trimci(:,2) = data(:,1); % mean values
            if contains(FileName,'ess','IgnoreCase',true)
                start_at = max(strfind(FileName,'_'))+1;
                C = LIMO.contrast{eval(FileName(start_at:end-4))}.C;
                df = rank(C); % rank of the relevant contrast
                trimci(:,1) = squeeze(trimci(:,2))-(finv(1-p./2*size(C,1),df,data(:,3)).*data(:,2));
                trimci(:,3) = squeeze(trimci(:,2))+(finv(1-p./2*size(C,1),df,data(:,3)).*data(:,2));
            else
                trimci(:,1) = squeeze(trimci(:,2))-(tinv(1-p./2,data(:,3)).*data(:,2));
                trimci(:,3) = squeeze(trimci(:,2))+(tinv(1-p./2,data(:,3)).*data(:,2));
            end

            % plot
            if strcmpi(LIMO.Analysis,'Time')
                if isfield(LIMO.data,'timevect')
                    xvect = LIMO.data.timevect;
                else
                    xvect = [];
                end

                if size(xvect,2) ~= size(toplot,2)
                    xvect              = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
                    LIMO.data.timevect = xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end

            elseif strcmpi(LIMO.Analysis,'Frequency')
                if isfield(LIMO.data,'freqlist')
                    xvect=LIMO.data.freqlist;
                else
                    xvect = [];
                end

                if size(xvect,2) ~= size(toplot,2)
                    xvect              = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
                    LIMO.data.freqlist = xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end

            elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                if length(time) > 1 && isfield(LIMO.data,'tf_times')
                    xvect = LIMO.data.tf_times;
                elseif length(freq) > 1 && isfield(LIMO.data,'tf_freqs')
                    xvect = LIMO.data.tf_freqs;
                else
                    xvect = [];
                end

                if length(time) > 1&& size(xvect,2) ~= size(data,1)
                    xvect              = linspace(LIMO.data.start,LIMO.data.end,size(data,1));
                    LIMO.data.tf_times =  xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                elseif length(freq) > 1 && size(xvect,2) ~= size(data,1)
                    xvect              = linspace(LIMO.data.lowf,LIMO.data.highf,size(data,1));
                    LIMO.data.tf_freqs =  xvect;
                    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
                end
            end

            figure;
            set(gcf,'Color','w')
            plot(xvect,squeeze(trimci(:,2)),'LineWidth',3);
            fillhandle = patch([xvect,fliplr(xvect)], [trimci(:,1)' fliplr(trimci(:,3)')], [1 0 0]);
            set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);% set edge color
            grid on; box on; axis tight
            h = axis;  hold on;
            plot(xvect,(sig./10+1).*h(3),'k.','MarkerSize',20)
            if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
                xlabel('Time in ms','FontSize',14)
                ylabel('Amplitude (A.U.)','FontSize',14)
            elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
                xlabel('Frequency in Hz','FontSize',14)
                ylabel('Spectral Power (A.U.)','FontSize',14)
            end
            if isempty(LIMO.design.electrode)
                title(sprintf('%s \n%s %s %s (%g)',mytitle,'Mean values',LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel),'FontSize',16); drawnow;
            else
                title(sprintf('%s \n%s virtual %s',mytitle,'Mean values',LIMO.Type(1:end-1)),'FontSize',16); drawnow;
            end
            assignin('base','Plotted_data',trimci);


        elseif contains(LIMO.design.name,'regression','IgnoreCase',true) && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
                contains(LIMO.design.name,'ANOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
                contains(LIMO.design.name,'ANCOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
            % --------------------------------------------------------------------------------

            % which variable(s) to plot
            % ----------------------
            if size(LIMO.design.X,2) >= 2
                if contains(FileName,'Condition_effect_') || ...
                        contains(FileName,'Covariate_effect_')
                    regressor = eval(FileName(18:end-4));
                    if contains(FileName,'Covariate_effect_')
                        regressor = regressor+sum(LIMO.design.nb_conditions);
                    end
                else
                    input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2)-1);
                    regressor = inputdlg(input_title,'Plotting option');
                end

                if isempty(regressor)
                    warning on
                    limo_errordlg('couldn''t figure out the regressor number/column, plot aborded')
                    return
                end

                try
                    if iscell(regressor)
                        regressor = sort(eval(cell2mat(regressor)));
                    end
                    if max(regressor) > size(LIMO.design.X,2)
                        errordlg('invalid regressor number');
                    end
                catch reginput_error
                    fprintf('error with regressor numbers/columns line 1373:\n %s',reginput_error.message)
                    return
                end
            else
                regressor = 1;
            end

            categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
            if max(regressor) == size(LIMO.design.X,2)
                tmp = regressor(1:end-1);
            else
                tmp = regressor;
            end

            if sum(tmp<=categorical) >=1 && sum(tmp>categorical) >=1
                limo_errordlg('you can''t plot categorical and continuous regressors together'); return
            end

            % load the effect
            % --------------
            data = load(FileName);
            data = data.(cell2mat(fieldnames(data)));
            if numel(size(data)) == 3 && size(data,2) == 1
                limo_errordlg('single time point detected, plot aborded'); return
            elseif numel(size(data)) == 4 && size(data,3) == 1
                limo_errordlg('single time point detected, plot aborded'); return
            end

            % which course plot to make
            % -------------------------
            if isempty(g.plot3type)
                extra = questdlg('Plotting data','Options','Original data','Modelled data','Adjusted data','Modelled data');
            else
                extra = g.plot3type;
                % allow typos
                if contains(extra,'orig','Ignorecase',true)
                    extra = 'Original';
                elseif contains(extra,'Model','Ignorecase',true)
                    extra = 'Modelled';
                elseif contains(extra,'Adj','Ignorecase',true)
                    extra = 'Adjusted';
                else
                    limo_errordlg(sprintf('input option ''%s'' invalid',extra)); return
                end
            end

            if isempty(extra)
                return
            elseif strcmpi(extra,'Original data')
                if regressor == size(LIMO.design.X,2)
                    limo_errordlg('you can''t plot adjusted mean for original data'); return
                end
            end

            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [~,channel,freq,time] = limo_display_reducedim(data,LIMO,g.channels,g.restrict,g.dimvalue);
                if length(freq) == 1; g.restrict = 'time';
                else; g.restrict = 'frequency'; end
                sig         = squeeze(single(mask(channel,freq,time))); %1D
            else
                [~,channel] = limo_display_reducedim(data,LIMO,g.channels);
                sig         = single(mask(channel,:));
            end
            sig(sig==0)=NaN;
            clear data

            % down to business
            % ----------------------
            probs = [p/2; 1-p/2];
            z     = norminv(probs);
            Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
            Yr    = Yr.Yr;

            if contains(extra,'Original','Ignorecase',true)
                if regressor <= length(LIMO.design.nb_conditions) && ...
                        LIMO.design.nb_conditions ~= 0 % for categorical variables
                    if length(LIMO.design.nb_conditions) == 1
                        start = 1;
                    else
                        start = sum(LIMO.design.nb_conditions(1:regressor-1));
                    end

                    for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
                        index{i} = find(LIMO.design.X(:,i));
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            data = squeeze(Yr(channel,freq,:,index{i}));
                        else
                            data = squeeze(Yr(channel,:,index{i}));
                        end
                        average(i,:) = nanmean(data,2);
                        se           = (nanstd(data,0,2) ./ sqrt(numel(index{i})));
                        ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Yr,2));
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Original subjects'' parameters at optimized channel');
                    end
                else % continuous variable
                    for i=max(regressor):-1:min(regressor)
                        index{i}           = find(LIMO.design.X(:,i));
                        [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
                        reg_values(i,:)    = LIMO.design.X(sorting_values,i);
                        if strcmpi(LIMO.Analysis,'Time-Frequency')
                            if strcmpi(g.restrict,'time')
                                continuous(i,:,:) = Yr(channel,freq,:,sorting_values);
                            else
                                continuous(i,:,:) = Yr(channel,:,time,sorting_values);
                            end
                        else
                            continuous(i,:,:) = Yr(channel,:,sorting_values);
                        end
                        clear mytitle
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Original subjects'' parameters \n sorted by regressor %g channel %s (%g)', i, LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle{i} = sprintf('Original subjects'' parameters, \n sorted by regressor %g at optimized channel', i);
                        end
                    end
                    remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
                    continuous(remove,:,:) = [];
                    reg_values(remove,:)   = [];
                end

            elseif contains(extra,{'Modelled','Modeled'},'Ignorecase',true)
                if exist('Betas.mat','file') % OLS & IRLS GLM
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        Betas = load('Betas.mat');
                        if strcmpi(g.restrict,'time')
                            Betas = squeeze(Betas.Betas(channel,freq,:,:));
                        else
                            Betas = squeeze(Betas.Betas(channel,:,time,:));
                        end
                        R     = eye(size(Yr,4)) - (LIMO.design.X*pinv(LIMO.design.X));
                    else
                        Betas = load('Betas.mat');
                        Betas = squeeze(Betas.Betas(channel,:,:));
                        R     = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                    end
                    Yh = (LIMO.design.X*Betas')'; % modelled data
                else % strcmpi(LIMO.design.method,'Generalized Welch's method')
                    [~,~,Yh,~,dfe] = limo_robust_1way_anova(squeeze(Yr(channel,:,:)),LIMO.design.X);
                    Res            = squeeze(Yr(channel,:,:))-Yh;
                end

                if regressor <= length(LIMO.design.nb_conditions) && ...
                        LIMO.design.nb_conditions ~= 0 % for categorical variables
                    if length(LIMO.design.nb_conditions) == 1
                        start = 1;
                    else
                        start = sum(LIMO.design.nb_conditions(1:regressor-1));
                    end

                    for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
                        index{i}     = find(LIMO.design.X(:,i));
                        data         = squeeze(Yh(:,index{i}));
                        average(i,:) = nanmean(data,2);
                        index{i}     = index{i}(find(~isnan(squeeze(Yr(channel,1,index{i}))))); %#ok<FNDSB>
                        if exist('R','var')
                            var      = diag(((R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')) / LIMO.model.model_df(2));
                        else
                            var      = diag(Res*Res')./dfe;
                        end
                        CI           = sqrt(var/size(index{i},1))*z';
                        ci(i,:,:)    = (repmat(nanmean(data,2),1,2)+CI)';
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Modelled subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Modelled subjects'' parameters at optimized channel');
                    end
                    clear Yr
                else % continuous variable
                    for i=max(regressor):-1:min(regressor)
                        index{i}           = find(LIMO.design.X(:,i));
                        [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
                        reg_values(i,:)    = LIMO.design.X(sorting_values,i);
                        continuous(i,:,:)  = Yh(:,sorting_values);
                        clear mytitle
                        if isempty(LIMO.design.electrode)
                            mytitle = sprintf('Modelled subjects'' parameters \n sorted by regressor %g channel %s (%g)', ...
                                i, LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Modelled subjects'' parameters \n sorted by regressor %g at optimized channel', i);
                        end
                    end
                    remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
                    continuous(remove,:,:) = [];
                    reg_values(remove,:)   = [];
                end

            elseif contains(extra,'Adjusted','Ignorecase',true)
                if length(LIMO.design.nb_conditions) == 1 && LIMO.design.nb_continuous == 0
                    warning on;
                    if exist('warndlg2','file')
                        warndlg2('Only one condition detected, no adjusted data possible');return
                    else
                        warndlg('Only one condition detected, no adjusted data possible');return
                    end
                end

                allvar = 1:size(LIMO.design.X,2)-1;
                allvar(regressor)=[]; % all but constant and columns of interest
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    Betas     = load('Betas.mat');
                    if strcmpi(g.restrict,'time')
                        Yr    = squeeze(Yr(channel,freq,:,:));
                        Betas = squeeze(Betas.Betas(channel,freq,:,:));
                    else
                        Yr    = squeeze(Yr(channel,:,time,:));
                        Betas = squeeze(Betas.Betas(channel,:,time,:));
                    end
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                else
                    Yr        = squeeze(Yr(channel,:,:));
                    Betas     = load('Betas.mat');
                    Betas     = squeeze(Betas.Betas(channel,:,:));
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                end
                Ya = Yr - confounds;
                clear Yr Betas confounds;

                if regressor <= length(LIMO.design.nb_conditions) && ...
                        LIMO.design.nb_conditions ~= 0 % for categorical variables
                    if length(LIMO.design.nb_conditions) == 1
                        start = 1;
                    else
                        start = sum(LIMO.design.nb_conditions(1:regressor-1));
                    end

                    for i = (start+LIMO.design.nb_conditions(regressor)-1):-1:start
                        index{i}     = find(LIMO.design.X(:,i));
                        data         = squeeze(Ya(:,index{i})); % use adjusted data
                        average(i,:) = nanmean(data,2);
                        se           = nanstd(data,0,2) ./ sqrt(numel(index{i}));
                        ci(i,:,:)    = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Ya,1));
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Adjusted subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Adjusted subjects'' parameters at  at optimized channel');
                    end
                    clear Yr

                else % continuous variable ; regressor value already + sum(LIMO.design.nb_conditions)
                    for i=max(regressor):-1:min(regressor)
                        index{i}           = find(LIMO.design.X(:,i));
                        [~,sorting_values] = sort(LIMO.design.X(index{i},i));  % continuous variable 3D plot
                        reg_values(i,:)    = LIMO.design.X(sorting_values,i);
                        continuous(i,:,:)  = Ya(:,sorting_values);
                        clear mytitle
                        if isempty(LIMO.design.electrode)
                            mytitle = sprintf('Adjusted subjects'' parameters \n sorted by regressor %g channel %s (%g)', i, LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Adjusted subjects'' parameters \n sorted by regressor %g at optimized channel', i);
                        end
                    end
                    remove                 = find(sum(reg_values == 0,2) == size(continuous,3));
                    continuous(remove,:,:) = [];
                    reg_values(remove,:)   = [];
                end
            else
                error('unspecified data type to plot ''Original'',''Modelled'' or ''Adjusted'' ')
            end

            % make the figure(s)
            % ------------------
            figure;set(gcf,'Color','w')
            if regressor <= length(LIMO.design.nb_conditions) && ...
                    LIMO.design.nb_conditions ~= 0 % for categorical variables
                brewcolours = limo_color_images(size(average,1));
                for i=1:size(average,1)
                    if strcmpi(LIMO.Analysis,'Time') || strcmpi(g.restrict,'time')
                        timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                        plot(timevect,squeeze(average(i,:)),'Color',brewcolours(i,:),'LineWidth',1.5); hold on
                        xlabel('Time in ms','FontSize',14)
                        ylabel('Amplitude (A.U.)','FontSize',14)
                    elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(g.restrict,'frequency')
                        freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                        plot(freqvect,squeeze(average(i,:)),brewcolours(i,:),'LineWidth',1.5); hold on
                        xlabel('Frequency in Hz','FontSize',14)
                        ylabel('Spectral Power (A.U.)','FontSize',14)
                    else
                        error('couldn''t figure out what dimension to plot')
                    end
                    x = squeeze(ci(i,1,:)); y = squeeze(ci(i,2,:));
                    fillhandle = patch([timevect fliplr(timevect)], [x',fliplr(y')], brewcolours(i,:));
                    set(fillhandle,'EdgeColor',brewcolours(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                end

                h = axis;
                if strcmpi(LIMO.Analysis,'Time')
                    plot(timevect,(sig./10+1).*h(3),'r*','LineWidth',2)
                else
                    plot(freqvect,(sig./10+1).*h(3),'r*','LineWidth',2)
                end

                axis tight; grid on; box on
                title(mytitle,'FontSize',16); drawnow;
                assignin('base','Plotted_data', average)
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                if strcmpi(LIMO.Analysis,'Time')
                    ylabel('Amplitude (A.U.)','FontSize',16)
                    xlabel('Time in ms','FontSize',16)
                else
                    ylabel('Spectral Power (A.U.)','FontSize',16)
                    xlabel('Frequency in Hz','FontSize',16)
                end

            else % 3D plots
                for i=1:size(continuous,1)
                    if i > 1; figure;set(gcf,'Color','w'); end
                    index = find(~isnan(squeeze(continuous(i,1,:))));
                    if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency')
                        if strcmpi(LIMO.Analysis,'Time')
                            timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in msec
                        else
                            timevect = LIMO.data.tf_times;
                        end
                        surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Time in ms','FontSize',16)
                        zlabel('Amplitude (A.U.)','FontSize',16)
                    else
                        surf(index,LIMO.data.freqlist,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Frequency in Hz','FontSize',16)
                        zlabel('Spectral Power (A.U.)','FontSize',16)
                    end
                    % --
                    axis tight; title(mytitle,'FontSize',14); drawnow;
                    xlabel('Sorted variable','FontSize',14)
                    try
                        set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                    catch label_err
                        warning on; warning('could not set X-labels:\n%s',label_err)
                    end
                end
            end


        elseif contains(LIMO.design.name,'Repeated','IgnoreCase',true)   % All stuffs for repeated measures ANOVA
            % -----------------------------------------------------------------------------

            if contains(FileName,'LIMO')
                error('Select summary stat file, nothing to infer from LIMO file')
            end

            % which summary stat
            % -------------------
            if ~isempty(g.sumstats) && any(strcmpi(g.sumstats,{'Mean','Trimmed'}))
                extra = g.sumstats;
            else
                if contains(LIMO.design.name,'robust','Ignorecase',true)
                    extra = 'Trimmed Mean';
                else
                    extra = 'Mean';
                end
                % let's not give GUI option and follow the design
                % extra = questdlg('Summarize data using:','Data plot option','Mean','Trimmed Mean','Mean');
                % if isempty(extra)
                %     return
                % end
            end

            if ~contains(FileName,'Rep_ANOVA_Interaction') && ...
                    ~contains(FileName,'Rep_ANOVA_Gp')
                % contains(FileName,{'Rep_ANOVA_Main','Rep_ANOVA'})
                % --------------------------------------------------

                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                if contains(FileName,'Main_effect','IgnoreCase',true)
                    index1     = strfind(FileName,'Main_effect')+length('Main_effect')+1;
                    index2     = max(strfind(FileName,'_'))-1;
                    effect_nb  = eval(FileName(index1:index2));
                elseif contains(FileName,'Interaction','IgnoreCase',true)
                    index1     = strfind(FileName,'Interaction')+length('Interaction')+1;
                    index2     = max(strfind(FileName,'_'))-1;
                    effect_nb  = eval(FileName(index1:index2));
                else
                    index1     = strfind(FileName,'Factor')+length('Factor')+1;
                    effect_nb  = eval(FileName(index1:end));
                end
                C              = LIMO.design.C{effect_nb};
                Data           = load(fullfile(LIMO.dir,'Yr.mat'));
                Data           = Data.(cell2mat(fieldnames(Data)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    [~,channel,freq,time] = limo_display_reducedim(Data,LIMO,g.channels,g.restrict,g.dimvalue);
                    Data              = squeeze(Data(channel,freq,time,:,:));
                    sig               = squeeze(single(mask(channel,freq,time)));
                else
                    [~,channel,freq,time] = limo_display_reducedim(Data,LIMO,g.channels);
                    Data              = squeeze(Data(channel,time,:,:)); % note freq/time variables have the same values
                    sig               = single(mask(channel,:));
                end
                sig(sig==0)=NaN;

                % compute differences between pairs using C and Cov
                n = size(Data,2);
                if strcmpi(extra,'Mean')
                    for time_or_freq = size(Data,1):-1:1
                        avg(time_or_freq,:) = nanmean(C*squeeze(Data(time_or_freq,:,:))',2);
                        S(time_or_freq,:,:) = nancov(squeeze(Data(time_or_freq,:,:)));
                    end
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original %s \n %s %s (%g)',mytitle,LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Original %s \n virtual %s',mytitle,LIMO.Type(1:end-1));
                    end
                else
                    g=floor((20/100)*n); %% compute for 20% trimmed mean
                    for time_or_freq = size(Data,1):-1:1
                        [v,indices]          = sort(squeeze(Data(time_or_freq,:,:))); % sorted data
                        TD(time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
                        avg(time_or_freq,:)  = nanmean(C*squeeze(TD(time_or_freq,:,:))',2);
                        v(1:g+1,:)           = repmat(v(g+1,:),g+1,1);
                        v(n-g:end,:)         = repmat(v(n-g,:),g+1,1); % winsorized data
                        [~,reorder]          = sort(indices);
                        for j = size(Data,3):-1:1
                            SD(:,j)          = v(reorder(:,j),j);
                        end % restore the order of original data
                        S(time_or_freq,:,:)  = cov(SD); % winsorized covariance
                    end
                    clear mytitle
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
                    end
                end

                % CI
                dfe = size(Data,2)-size(Data,3)+1;
                % c = avg + 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C'))); % uses Bonferoni inequality
                % b = avg - 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C')));
                bound = (abs(tinv(p./(2*size(C,1)),dfe)).*diag((sqrt(C*squeeze(S(time_or_freq,:,:))*C'))));
                c = avg + repmat(bound', [length(avg),1]);
                b = avg - repmat(bound', [length(avg),1]);

                % do the figure
                brewcolours = limo_color_images(size(avg,2));
                if strcmpi(LIMO.Analysis,'Time')
                    if isfield(LIMO.data,'timevect')
                        xvect = LIMO.data.timevect;
                    else
                        LIMO.data.timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                        save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.timevect;
                    end
                elseif strcmpi(LIMO.Analysis,'Frequency')
                    if isfield(LIMO.data,'timevect')
                        xvect = LIMO.data.freqlist;
                    else
                        LIMO.data.freqlist = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,2));
                        save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.freqlist;
                    end
                elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                    if length(time) > 1
                        if isfield(LIMO.data,'tf_times')
                            xvect = LIMO.data.tf_times;
                        else
                            LIMO.data.tf_times = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                            save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.tf_times;
                        end
                    elseif length(freq) > 1
                        if isfield(LIMO.data,'tf_freqs')
                            xvect = LIMO.data.tf_freqs;
                        else
                            LIMO.data.tf_freqs = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,2));
                            save(fullfile(LIMO.dir,'LIMO'),'LIMO'); xvect = LIMO.data.tf_freqs;
                        end
                    end
                end

                figure;set(gcf,'Color','w')
                for cond = 1:size(c,2)
                    plot(xvect,avg(:,cond)','LineWidth',3,'color',brewcolours(cond,:));
                    fillhandle = patch([xvect fliplr(xvect)], [c(:,cond)',fliplr(b(:,cond)')], brewcolours(cond,:));
                    set(fillhandle,'EdgeColor',brewcolours(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                    hold on
                end
                grid on; box on; axis tight; hold on;
                h = axis;  plot(xvect,(sig./10+1).*h(3),'k.','MarkerSize',20)
                if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
                    xlabel('Time in ms','FontSize',14)
                    ylabel('Amplitude (A.U.)','FontSize',14)
                elseif strcmpi(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
                    xlabel('Frequency in Hz','FontSize',14)
                    ylabel('Spectral Power (A.U.)','FontSize',14)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                trimci = [c ; avg ; b];
                assignin('base','Plotted_data',trimci);

                % ----------------------
            elseif contains(FileName,'Rep_ANOVA_Gp')  %% plot pairs of gp differences

                % -------------------
                % which ERP to make
                % ------------------
                extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Original');
                if isempty(extra)
                    return;
                end
                % -----------------------
                Rep_ANOVA_Gp_effect = load(FileName);
                Rep_ANOVA_Gp_effect = Rep_ANOVA_Gp_effect.(cell2mat(fieldnames(Rep_ANOVA_Gp_effect)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    [~,channel,freq,time] = limo_display_reducedim(Rep_ANOVA_Gp_effect,LIMO,g.channels,g.restrict,g.dimvalue);
                    Rep_ANOVA_Gp_effect   = squeeze(Rep_ANOVA_Gp_effect(channel,freq,time,:,:));
                    sig                   = squeeze(single(mask(channel,freq,time)));
                else
                    [~,channel,freq,time] = limo_display_reducedim(Rep_ANOVA_Gp_effect,LIMO,g.channels); %#ok<ASGLU>
                    Rep_ANOVA_Gp_effect   = squeeze(Rep_ANOVA_Gp_effect(channel,time,:,:)); % note freq/time variables have the same values
                    sig                   = single(mask(channel,:));
                end

                % check channel to plot
                if size(Rep_ANOVA_Gp_effect,1) > 1
                    if isempty(channel)
                        [v,e] = max(Rep_ANOVA_Gp_effect(:,:,1));
                        [~,c] = max(v); channel = e(c);
                    else
                        if ischar(channel)
                            channel = str2double(channel);
                        end

                        if length(channel) > 1
                            limo_errordlg('1 channel only can be plotted'); return
                        elseif channel > size(Rep_ANOVA_Gp_effect,1)
                            limo_errordlg('channel number invalid'); return
                        end
                    end
                else
                    if length(LIMO.design.electrode) == 1
                        channel = LIMO.design.electrode;
                    else
                        channel = 1;  % accomodates the fact that all matrices have the channel dim (even = 1)
                    end
                end
                clear Rep_ANOVA_Gp_effect

                % compute the pair-wise differences and plot
                Yr = load(fullfile(LIMO.dir,'Yr.mat'));
                Yr = Yr.Yr;
                if strcmpi(extra,'Original')
                    Data = mean(squeeze(Yr(channel,:,:,:)),3);
                    combinations = nchoosek(1:LIMO.design.nb_conditions,2);
                    for d = size(combinations,1):-1:1
                        Effect(:,d) = nanmean(Data(:,find(LIMO.data.Cat == combinations(d,1))),2) - nanmean(Data(:,find(LIMO.data.Cat == combinations(d,2))),2); %#ok<FNDSB>
                    end
                end

                % design matrix
                X = zeros(size(Yr,3),LIMO.design.nb_conditions+1);
                X(:,end) = 1;
                for i=1:LIMO.design.nb_conditions
                    X(find(LIMO.data.Cat == i),i) = 1; %#ok<FNDSB>
                end

                % data again
                Y = nanmean(squeeze(Yr(channel,:,:,:)),3);
                X = X(find(~isnan(Y(1,:))),:); %#ok<FNDSB>
                Y = Y(:,find(~isnan(Y(1,:))))'; %#ok<FNDSB>
                if strcmpi(extra,'Modelled')
                    beta = pinv(X)*Y; Yhat = X*beta;
                    combinations = nchoosek(1:LIMO.design.nb_conditions,2);
                    for d = 1:size(combinations,1)
                        Effect(:,d) = nanmean(Yhat(find(X(:,combinations(d,1))),:),1)' - nanmean(Yhat(find(X(:,combinations(d,2))),:),1)'; %#ok<FNDSB>
                    end
                end

                Res    = (Y'*(eye(size(Y,1)) - (X*pinv(X)))*Y);
                df     = size(Y,1)-rank(X); t = tcdf(1-p,df);
                sigma2 = sum((Res.^2./df),2);
                v      = t.*sqrt(sigma2 ./ norm(X(:,1:end-1)).^2);
                b      = Effect - v(channel,:);
                c      = Effect + v(channel,:);
                if channel == 1
                    if length(LIMO.design.electrode) == 1
                        if strcmpi(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        end
                    else
                        if strcmpi(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n optimized channel');
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n optimized channel');
                        end
                    end
                else
                    if strcmpi(extra,'Original')
                        mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
                    end
                end

                figure;set(gcf,'Color','w'); hold on
                brewcolours = limo_color_images(size(Effect,2));
                if strcmpi(LIMO.Analysis,'Time')
                    xvect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
                else
                    xvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                end

                for d = size(combinations,1):-1:1
                    plot(xvect,Effect(:,d),'LineWidth',3,'Color',brewcolours(d,:));
                    fillhandle = patch([xvect fliplr(xvect)], [c(:,d)',fliplr(b(:,d)')], brewcolours(d,:));
                    set(fillhandle,'EdgeColor',brewcolours(d,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                    Gp_difference{d} = [c(:,d)' ; Effect(:,d)' ; b(:,d)'];
                end
                grid on; box on; axis tight; hold on;
                h = axis; sig(sig==0)=NaN;
                plot(xvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                assignin('base','Plotted_data',Gp_difference);

                % -------------------------
            elseif contains(FileName,'Rep_ANOVA_Interaction') % Gp * Repeated measures - plot differences btween condition per gp
                % ------------------------

                Rep_ANOVA_Interaction_with_gp = load(FileName);
                Rep_ANOVA_Interaction_with_gp = Rep_ANOVA_Interaction_with_gp.(cell2mat(fieldnames(Rep_ANOVA_Interaction_with_gp)));
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    [~,channel,freq,time]         = limo_display_reducedim(Rep_ANOVA_Interaction_with_gp,LIMO,g.channels,g.restrict,g.dimvalue);
                    Rep_ANOVA_Interaction_with_gp = squeeze(Rep_ANOVA_Interaction_with_gp(channel,freq,time,:,:));
                    sig                           = squeeze(single(mask(channel,freq,time)));
                else
                    [~,channel,freq,time]         = limo_display_reducedim(Rep_ANOVA_Interaction_with_gp,LIMO,g.channels); %#ok<ASGLU>
                    Rep_ANOVA_Interaction_with_gp = squeeze(Rep_ANOVA_Interaction_with_gp(channel,time,:,:)); % note freq/time variables have the same values
                    sig                           = single(mask(channel,:));
                end
                sig(sig==0)=NaN;

                if size(Rep_ANOVA_Interaction_with_gp,1) > 1
                    if isempty(channel)
                        [v,e] = max(Rep_ANOVA_Interaction_with_gp(:,:,1));
                        [~,c] = max(v); channel = e(c);
                    else
                        if ischar(channel)
                            channel = str2double(channel);
                        end

                        if length(channel) > 1
                            limo_error('1 channel only can be plotted'); return
                        elseif channel > size(Rep_ANOVA_Interaction_with_gp,1)
                            limo_error('channel number invalid'); return
                        end
                    end
                else
                    channel = 1;
                end

                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                effect_nb = str2double(FileName(strfind(FileName,'Factor_')+7:strfind(FileName,'Factor_')+6+strfind(FileName(strfind(FileName,'Factor_')+6:end),'_')));
                C         = LIMO.design.C{effect_nb};
                Yr        = load(fullfile(LIMO.dir,'Yr.mat'));
                Yr        = Yr.Yr;
                Data      = squeeze(Yr(channel,:,:,:));

                % compute differences between pairs using C and Cov
                if strcmpi(extra,'Mean')
                    for gp = 1:LIMO.design.nb_conditions
                        index = find(LIMO.data.Cat==gp);
                        for time=size(Data,1):-1:1
                            avg(gp,time,:) = nanmean(C*squeeze(Data(time,index,:))',2);
                            S(gp,time,:,:) = cov(squeeze(Data(time,index,:)));
                        end
                    end
                    if ~isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Original %s \n optimized channel',mytitle);
                    end
                else
                    for gp = 1:LIMO.design.nb_conditions
                        index = find(LIMO.data.Cat==gp);
                        n     = length(index);
                        g     = floor((20/100)*n); %% compute for 20% trimmed mean
                        for time_or_freq=size(Data,1):-1:1
                            [v,indices]             = sort(squeeze(Data(time_or_freq,index,:))); % sorted data
                            TD(gp,time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
                            avg(gp,time_or_freq)    = nanmean(C*squeeze(TD(gp,time_or_freq,:,:))',2);
                            v(1:g+1,:)              = repmat(v(g+1,:),g+1,1);
                            v(n-g:end,:)            = repmat(v(n-g,:),g+1,1); % winsorized data
                            [~,reorder]             = sort(indices);
                            for j = 1:size(Data,3)
                                SD(:,j)             = v(reorder(:,j),j); % restore the order of original data
                            end
                            S(gp,time_or_freq,:,:)  = cov(SD); % winsorized covariance
                        end
                        clear SD
                    end
                    if ~isempty(LIMO.design.electrode)
                        mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
                    end
                end

                figure; set(gcf,'Color','w'); hold on
                if strcmpi(LIMO.Analysis,'Time')
                    xvect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
                else
                    xvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                end

                colorindex = 1;
                trimci      = cell(size(avg,1)*size(avg,3),1);
                brewcolours = limo_color_images(size(avg,1)*size(avg,3));
                for gp = 1:LIMO.design.nb_conditions
                    % get the variance per comparison
                    for frame = size(avg,2):-1:1
                        varc(:,frame) = diag(sqrt(C*squeeze(S(gp,frame,:,:))*C'));
                    end
                    % plot each comparison
                    for c = 1:size(avg,3)
                        plot(xvect,squeeze(avg(gp,:,c)),'Color',brewcolours(colorindex,:),'LineWidth',3);
                        index      = find(LIMO.data.Cat==gp);
                        dfe        = size(Data,2)-length(index)+1;
                        up         = avg(gp,:,c) + tinv(p./(2*size(C,1)),dfe).* varc(c,:);
                        down       = avg(gp,:,c) - tinv(p./(2*size(C,1)),dfe).* varc(c,:);
                        fillhandle = patch([xvect fliplr(xvect)], [up,fliplr(down)], brewcolours(colorindex,:));
                        set(fillhandle,'EdgeColor',brewcolours(colorindex,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                        trimci{colorindex} = [up ; squeeze(avg(gp,:,c)); down];
                        colorindex = colorindex + 1;
                    end
                end

                grid on; box on; axis tight
                h = axis;  hold on;
                plot(xvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                ylabel('Amplitude (A.U.)','FontSize',14)
                if strcmpi(LIMO.Analysis,'Time')
                    xlabel('Time in ms','FontSize',14)
                else
                    ylabel('Frequency in Hz','FontSize',14)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                assignin('base','Plotted_data',trimci);
            end

            % -------------------------
        else % might be a summary data file
            % -------------------------
            try
                limo_add_plots({FileName,''},LIMO);
            catch no_plot
                limo_errordlg('course plot failed because %s',no_plot.message)
                return
            end
        end
    end % closes type

    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % -------------------                            LATERALIZATION STUFF
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------

elseif strcmpi(LIMO.Level,'LI')

    [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,[],[]);

    if Type == 1
        %--------------------------
        % imagesc of the results
        %--------------------------
        if sum(mask(:)) == 0
            limo_errordlg('no values under threshold parameter','no significant effect');
        else
            scale = M.*mask;
            if min(scale(:))<0
                scale(scale==0)=min(scale(:))+(min(scale(:))/10);
            else
                scale(scale==0)=NaN;
            end

            figure; set(gcf,'Color','w');
            timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(M,2));
            imagesc(timevect,1:size(M,1),scale);
            title(mytitle,'FontSize',18);
            color_images_(scale,LIMO);
            assignin('base','Plotted_data',scale)
            assignin('base','Mask_of_sig',mask)
        end

    elseif Type == 2
        %--------------------------
        % topoplot
        %--------------------------
        if sum(mask(:)) == 0
            warndlg('no values under threshold','no significant effect');
        else
            EEG.data = M.*mask;
            EEG.setname = 'Lateralization Map';
            EEG.pnts = size(EEG.data,2);
            EEG.xmin = LIMO.data.start/1000; % in sec
            EEG.xmax = LIMO.data.end/1000;   % in sec
            EEG.times =  (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
            EEG.trials = 1;
            EEG.chanlocs = LIMO.data.chanlocs;
            EEG.nbchan = size(EEG.data,1);
            pop_topoplot(EEG);
            assignin('base','Plotted_data',EEG.data)
        end

    elseif Type == 3
        disp('no ERP Plots to Lateralization'); % we could but nobody asked for it
    end
end % closes if LIMO.level
end % closes the function

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% color map
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function color_images_(scale,LIMO)

scale(scale==0) = NaN;
cc              = limo_color_images(scale); % get a color map commensurate to that
colormap(cc);

set(gca,'XMinorTick','on','LineWidth',2)
if isfield(LIMO.data,'expected_chanlocs')
    set(gca,'YTick',1:length(LIMO.data.expected_chanlocs));
else
    set(gca,'YTick',1:length(LIMO.data.chanlocs));
end

ylabel('Electrodes','FontSize',14);
if strcmpi(LIMO.Analysis,'Time')
    xlabel('Time in ms','FontSize',16)
elseif strcmpi(LIMO.Analysis,'Frequency')
    xlabel('Frequency in Hz','FontSize',16)
end

if LIMO.Level == 1
    if strcmpi(LIMO.data.chanlocs,'Components')
        label_electrodes = [];
    else
        for i = length(LIMO.data.chanlocs):-1:1
            if isfield(LIMO.data,'expected_chanlocs')
                label_electrodes{i} = LIMO.data.expected_chanlocs(i).labels;
            else
                label_electrodes{i} = LIMO.data.chanlocs(i).labels;
            end
        end
    end
else
    if isempty(LIMO.design.electrode)
        for i = length(LIMO.data.chanlocs):-1:1
            label_electrodes{i} = LIMO.data.chanlocs(i).labels;
        end
    else
        if length(LIMO.design.electrode) == 1
            label_electrodes = LIMO.design.electrode;
        else
            label_electrodes = ' ';
            ylabel('optimized channel','FontSize',14);
        end
    end
end
set(gca,'YTickLabel', label_electrodes);
end

function call_topolot(EEG,FileName,Domain)

EEG.pnts     = size(EEG.data,2);
EEG.nbchan   = size(EEG.data,1);
EEG.trials   = 1;

if strcmpi(FileName,'R2.mat') || strcmpi(FileName,'R2')
    newname = 'R^2';
else
    newname = [FileName(1:min(strfind(FileName,'_'))-1),FileName(max(strfind(FileName,'_'))+1:end-4)];
end

if strcmpi(Domain,'Time')
    if ~isfield(EEG,'setname')
        if contains(FileName,'con','IgnoreCase',true)
            EEG.setname = sprintf('%s - T values',newname);
        else
            EEG.setname = sprintf('%s - F values',newname);
        end
    end
    pop_topoplot(EEG);
    % set(gca,'Colormap',limo_color_images(EEG.data),'CLim',[min(EEG.data(:)),max(EEG.data(:))])
else % freq
    N = size(EEG.freq,2);
    figure;
    for f=1:N
        if N<=6
            subplot(1,N,f)
        else
            subplot(ceil(N/6),6,f);
        end
        [~,ind] = min(abs(EEG.freq-EEG.freq(f)));
        opt = {'electrodes','on','maplimits','maxmin','verbose','off','colormap', limo_color_images(EEG.data(:,ind))};
        topoplot(EEG.data(:,ind),EEG.chanlocs,opt{:});
        if isfield(EEG,'setname')
            title(sprintf('Frequency %g Hz from \n%s',round(EEG.freq(ind)),EEG.setname));
        else
            title(sprintf('Frequency %g Hz from %s - F values',round(EEG.freq(ind)),newname));
        end
    end
end
end


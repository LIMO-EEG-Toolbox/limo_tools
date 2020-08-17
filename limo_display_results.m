function limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag,varargin)

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
%   flag      = indicates to allow surfing the figure (1) or not (0)
%
% OPTIONAL INPUTS  (Usage: {''key'', value, ... })
% 'channels' : Provide the index of the channel to be used.
% 'regressor': Provide the index of the regressor to be used.
% 'plot3type': Type of plots to show when 'Type' is 3. Select between {'Original', 'Modeled', 'Adjusted'}
%
% Although the function is mainly intented to be used via the GUI, some figures
% can be generated automatically, for instance limo_display_results(1,'R2.mat',pwd,0.05,5,LIMO,0);
% would load the R2.mat file from the current directory, and plot all
% electrodes/time frames F values thresholded using tfce at alpha 0.05
% topoplot and ERP like figures can't be automated since they require user
% input
%
% Cyril Pernet, Guillaume Rousselet, Carl Gaspar,  
% Nicolas Chauveau, Andrew Stewart, Ramon Martinez-Cancino
%
% see also limo_stat_values limo_display_image topoplot limo_course_plot
% ----------------------------------------------------------------------
%  Copyright (C) LIMO Team 2019

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
    disp('limo_display_results() error: calling convention {''key'', value, ... } error'); return;
end

try g.channels;   catch, g.channels  = [];  end % No default values
try g.regressor;  catch, g.regressor = [];  end % No default values
try g.plot3type;  catch, g.plot3type = [];  end % No default values

toplot = load(fullfile(PathName,FileName));
toplot = toplot.(cell2mat(fieldnames(toplot)));
if nargin <= 6
    flag = 1;
end

choice = 'use theoretical p values'; % threshold based on what is computed since H0 is used for clustering
                     % see limo_stat_values - discontinuated empirical threshold (misleading)

if LIMO.design.bootstrap == 0
    if MCC == 2
        errordlg2('Clustering thresholding necessitates boostrap - invalid choice');
    elseif MCC == 3
        errordlg2('TFCE thresholding necessitates boostrap - invalid choice');
    elseif MCC == 4
        errordlg2('Maximum stat thresholding necessitates bootstrap - invalid choice');
    end
    MCC = 1;
end

if LIMO.design.bootstrap == 1 && LIMO.design.tfce == 0 && MCC == 3
    errordlg2('TFCE thresholding hasn''t been computed - invalid choice');
    MCC =1;
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
            
            if  strcmp(LIMO.design.type_of_analysis,'Mass-univariate')
                
                % univariate results from 1st level analysis
                % ------------------------------------------
                
                % if previously plotted recover data from the cache
                data_cached = 0;
                if isfield(LIMO,'cache')
                    try
                        if strcmp(LIMO.cache.fig.name, FileName) && ...
                                LIMO.cache.fig.MCC == MCC && ...
                                LIMO.cache.fig.threshold == p
                            
                            disp('using cached data');
                            mask = LIMO.cache.fig.mask;
                            if isempty(mask)
                                data_cached = 0;
                            elseif sum(mask(:)) == 0
                                warndlg('  no values under threshold  ','no significant effect','modal');
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
                    
                    [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,choice);
                    
                    if isempty(mask)
                        disp('no values computed'); return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        LIMO.cache.fig.name       = FileName;
                        LIMO.cache.fig.MCC        = MCC;
                        LIMO.cache.fig.stats      = [];
                        LIMO.cache.fig.threshold  = p;
                        LIMO.cache.fig.pval       = M;
                        LIMO.cache.fig.mask       = mask;
                        LIMO.cache.fig.title      = mytitle;
                        % do an exception for designs with just the constant 
                        if strcmp(FileName,'R2.mat') && size(LIMO.design.X,2)==1
                            mask = ones(size(mask)); LIMO.cache.fig.mask = mask;
                            mytitle = 'R^2 Coef unthresholded'; LIMO.cache.fig.title = mytitle;
                            save LIMO LIMO
                        else
                            save LIMO LIMO
                            return
                        end
                    else
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                    end
                    
                    if contains(FileName,'R2','IgnoreCase',true)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,1)); % plot R2 values instead of F
                        else
                            toplot = squeeze(toplot(:,:,1));
                        end
                        assignin('base','R2_values',squeeze(toplot(:,:,2)))
                        
                    elseif contains(FileName,'Condition_effect','IgnoreCase',true) || ...
                            contains(FileName,'Covariate_effect','IgnoreCase',true) || ...
                            contains(FileName,'Interaction_effect','IgnoreCase',true) || ...
                            contains(FileName,'semi_partial_coef','IgnoreCase',true)
                        
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,1)); % plot F values
                        else
                            toplot = squeeze(toplot(:,:,1));
                        end
                        
                        if contains(FileName,'semi_partial_coef','IgnoreCase',true)
                            assignin('base','semi_partial_coef',toplot)
                        else
                            assignin('base','F_values',toplot)
                        end
                        
                    elseif strcmp(FileName(1:4),'con_')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,4)); % plot T values
                        else
                            toplot = squeeze(toplot(:,:,4));
                        end
                        assignin('base','T_values',toplot)
                        
                    elseif strcmp(FileName(1:4),'ess_')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(toplot(:,:,:,end-1)); % plot F values
                        else
                            toplot = squeeze(toplot(:,:,end-1));
                        end
                        assignin('base','F_values',toplot)
                        
                    else
                        errordlg2('file not supported');
                        return
                    end
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
                        save LIMO LIMO
                    end
                    
                    if ndims(toplot)==3
                        limo_display_image_tf(LIMO,toplot,mask,mytitle,flag);
                    else
                        limo_display_image(LIMO,toplot,mask,mytitle,flag)
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
                            warndlg('  no values under threshold  ','no significant effect','modal');
                            return
                        else
                            toplot = squeeze(toplot(:,1)); % plot R2 values instead of F
                            assignin('base','F_values',F_values)
                            assignin('base','p_values',M)
                            assignin('base','mask',mask)
                            clear R2
                        end

                    else

                        if strcmp(FileName(end-6:end),'_EV.mat')
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
                            warndlg('  no values under threshold  ','no significant effect','modal');
                            return
                        else
                            if strcmp(choice,'Roy')
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
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                warndlg('topoplot not supported for 3D data')
            else
                EEG.trials = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                if strcmp(LIMO.Analysis,'Time')
                    EEG.xmin = LIMO.data.start / 1000;% in msec
                    EEG.xmax = LIMO.data.end / 1000;  % in msec
                    EEG.times = LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000; % in sec;
                    if length(EEG.times) > 2
                       EEG.times = [EEG.times(1) EEG.times(end)];
                    end
                elseif strcmp(LIMO.Analysis,'Frequency')
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
                            errordlg('slected frequency out of bound'); return
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
                    EEG.pnts    = size(EEG.data,2);
                    EEG.nbchan  = size(EEG.data,1);
                    if strcmp(LIMO.Analysis,'Time')
                        pop_topoplot(EEG);
                    else
                        N = size(EEG.freq,2); figure;
                        for f=1:N
                            if N<=6
                                subplot(1,N,f)
                            else
                                subplot(ceil(N/6),6,f);
                            end
                            [~,ind] = min(abs(LIMO.data.freqlist-EEG.freq(f)));
                            topoplot(EEG.data(:,ind),EEG.chanlocs);
                            title([num2str(LIMO.data.freqlist(ind)) ' Hz'],'FontSize',12)
                        end
                    end
                    assignin('base',EEG.setname(1:2),EEG.data);
                elseif contains(FileName,'Condition_effect','IgnoreCase',true) || ...
                        contains(FileName,'Covariate_effect','IgnoreCase',true) || ...
                        contains(FileName,'Interaction_effect','IgnoreCase',true)  || ...
                        contains(FileName,'Condition_effect','IgnoreCase',true) 
                    EEG.data    = squeeze(toplot(:,:,1));
                    EEG.pnts    = size(EEG.data,2);
                    EEG.nbchan  = size(EEG.data,1);
                    EEG.setname = sprintf('%s %s - F values',FileName(1:min(strfind(FileName,'_'))-1),FileName(max(strfind(FileName,'_'))+1:end-4));
                    if strcmp(LIMO.Analysis,'Time')
                        pop_topoplot(EEG);
                    else
                        N = size(EEG.freq,2); figure;
                        for f=1:N
                            if N<=6
                                subplot(1,N,f)
                            else
                                subplot(ceil(N/6),6,f);
                            end
                            [~,ind] = min(abs(LIMO.data.freqlist-EEG.freq(f)));
                            topoplot(EEG.data(:,ind),EEG.chanlocs);
                            title([num2str(LIMO.data.freqlist(ind)) ' Hz'],'FontSize',12)
                        end
                    end
                    assignin('base','F values',EEG.data);
                elseif strcmp(FileName,'semi partial_coef.mat')
                    regressor = str2double(cell2mat(inputdlg('which regressor(s) to plot (e.g. 1:3)','Plotting option')));
                    if max(regressor) > size(toplot,3); errordlg('error in regressor number'); return; end
                    for b = regressor
                        EEG.data     = squeeze(toplot(:,:,b,1));
                        EEG.pnts     = size(EEG.data,2);
                        EEG.nbchan   = size(EEG.data,1);
                        EEG.setname  = sprintf('semi partial coef R2 values variable %g',b);
                        EEG.pnts     = size(EEG.data,2);
                        EEG.times    = LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000;
                        EEG.trials   = 1;
                        EEG.chanlocs = LIMO.data.chanlocs;
                        EEG.nbchan   = size(EEG.data,1);
                        if strcmp(LIMO.Analysis,'Time')
                            EEG.xmin = LIMO.data.start/1000;
                            EEG.xmax = LIMO.data.end/1000;
                            pop_topoplot(EEG);
                        else
                            EEG.xmin = LIMO.data.freqlist(1);
                            EEG.xmax = LIMO.data.freqlist(end);
                            N = size(EEG.freq,2); figure;
                            for f=1:N
                                if N<=6
                                    subplot(1,N,f)
                                else
                                    subplot(ceil(N/6),6,f);
                                end
                                [~,ind] = min(abs(LIMO.data.freqlist-EEG.freq(f)));
                                topoplot(EEG.data(:,ind),EEG.chanlocs);
                                title([num2str(LIMO.data.freqlist(ind)) ' Hz'],'FontSize',12)
                            end
                        end
                        tmp_name = sprintf('semi_partial_coef_%g',b);
                        assignin('base',tmp_name,EEG.data);
                    end
                elseif strcmp(FileName(1:4),'con_')
                    EEG.data    = squeeze(toplot(:,:,4));
                    EEG.pnts    = size(EEG.data,2);
                    EEG.nbchan  = size(EEG.data,1);
                    EEG.setname = ['Contrast ',FileName(5:end-4),' -- T values'];
                    if strcmp(LIMO.Analysis,'Time')
                        pop_topoplot(EEG);
                    else
                        N = size(EEG.freq,2); figure;
                        for f=1:N
                            if N<=6
                                subplot(1,N,f)
                            else
                                subplot(ceil(N/6),6,f);
                            end
                            [~,ind] = min(abs(LIMO.data.freqlist-EEG.freq(f)));
                            topoplot(EEG.data(:,ind),EEG.chanlocs);
                            title([num2str(LIMO.data.freqlist(ind)) ' Hz'],'FontSize',12)
                        end
                    end
                    assignin('base',EEG.setname(1:8),EEG.data);
                elseif strcmp(FileName(1:4),'ess_')
                    EEG.data    = squeeze(toplot(:,:,end-1));
                    EEG.pnts    = size(EEG.data,2);
                    EEG.nbchan  = size(EEG.data,1);
                    EEG.setname = ['Contrast ',FileName(5:end-4),' -- F values'];
                    if strcmp(LIMO.Analysis,'Time')
                        pop_topoplot(EEG);
                    else
                        N = size(EEG.freq,2); figure;
                        for f=1:N
                            if N<=6
                                subplot(1,N,f)
                            else
                                subplot(ceil(N/6),6,f);
                            end
                            [~,ind] = min(abs(LIMO.data.freqlist-EEG.freq(f)));
                            topoplot(EEG.data(:,ind),EEG.chanlocs);
                            title([num2str(LIMO.data.freqlist(ind)) ' Hz'],'FontSize',12)
                        end
                    end
                    assignin('base',EEG.setname(1:8),EEG.data);
                else
                    disp('file not supported');
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
                regressor = inputdlg(input_title,'Plotting option');
            else
                regressor = g.regressor;
            end
            
            if isempty(regressor); disp('selection aborded'); return; end
            regressor = cell2mat(regressor);
            if isempty(regressor); disp('selection aborded'); return; end
            
            try regressor = sort(str2double(regressor));
                if max(regressor) > size(LIMO.design.X,2); errordlg('invalid regressor number'); end
            catch ME
                error('can''t select this regressor: %s',ME.message); 
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
            
            % which ERP to make
            % ------------------
            if isempty(g.plot3type) && ~any(strcmp(g.plot3type,{'Original','Modelled','Adjusted'}))
                extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Adjusted','Adjusted');
            else
                extra = g.plot3type;
            end
            if isempty(extra)
                return;
            elseif strcmp(extra,'Original')
                if regressor == size(LIMO.design.X,2)
                    errordlg('you can''t plot adjusted mean for original data'); return
                end
            end
            
            % timing /frequency info
            % -----------------------
            if strcmp(LIMO.Analysis,'Time')
                timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
            elseif strcmp(LIMO.Analysis,'Frequency')
                freqvect=LIMO.data.freqlist';
            elseif strcmp(LIMO.Analysis,'Time-Frequency')
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
            if strcmp(LIMO.Analysis,'Time-Frequency')
                disp('loading the 4D data ...')
                frequency = inputdlg('which Frequency to plot','Plotting option');
            else
                frequency = [];
            end
            
            if strcmp(channel,'') || strcmp(frequency,'')
                disp('looking for max'); 
                R2 = load(fullfile(LIMO.dir,'R2.mat'));
                R2 = R2.R2;
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    tmp = squeeze(R2(:,:,:,1)); clear R2
                    [e,f,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                    if length(e) ~= 1; e = e(1); f = f(1); end
                    
                    if strcmp(channel,'') 
                        channel = e;
                    else
                        channel = eval(cell2mat(channel)); 
                    end
                    if size(channel) > 1
                        errordlg('invalid channel choice'); return
                    elseif channel > size(LIMO.data.chanlocs,2) || channel < 1
                        errordlg('invalid channel number'); return
                    end
                    
                    if strcmp(frequency,'')
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
                    errordlg('invalid channel choice'); return
                elseif channel > size(LIMO.data.chanlocs,1) || channel < 1
                    errordlg('invalid channel number'); return
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
                if strcmp(LIMO.Analysis,'Time-Frequency') && isfield(LIMO.cache,'ERPplot')
                    
                    if mean([LIMO.cache.ERPplot.channel == channel ...
                            LIMO.cache.ERPplot.regressor == regressor ...
                            LIMO.cache.ERPplot.frequency == frequency]) == 1 ...
                            && strcmp('LIMO.cache.ERPplot.extra',extra)
                        
                        if sum(regressor <= categorical) == length(regressor)
                            average = LIMO.cache.ERPplot.average;
                            ci = LIMO.cache.ERPplot.ci;
                            mytitle = LIMO.cache.ERPplot.title;
                            disp('using cached data');
                            data_cached = 1;
                        else
                            continuous = LIMO.cache.ERPplot.continuous;
                            mytitle = LIMO.cache.ERPplot.title;
                            disp('using cached data');
                            data_cached = 1;
                        end
                    end
                    
                elseif strcmp(LIMO.Analysis,'Time') && isfield(LIMO.cache,'ERPplot') || ...
                        strcmp(LIMO.Analysis,'Frequency') && isfield(LIMO.cache,'ERPplot')
                    
                    if length(LIMO.cache.ERPplot.regressor) == length(channel)
                        if mean([LIMO.cache.ERPplot.channel == channel ...
                                LIMO.cache.ERPplot.regressor == regressor]) == 1  ...
                                && strcmp('LIMO.cache.ERPplot.extra',extra)
                            
                            if sum(regressor <= categorical) == length(regressor)
                                average = LIMO.cache.ERPplot.average;
                                ci = LIMO.cache.ERPplot.ci;
                                mytitle = LIMO.cache.ERPplot.title;
                                disp('using cached data');
                                data_cached = 1;
                            else
                                continuous = LIMO.cache.ERPplot.continuous;
                                mytitle = LIMO.cache.ERPplot.title;
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
                
                if strcmp(extra,'Original')
                    Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            if strcmp(LIMO.Analysis,'Time-Frequency')
                                data = squeeze(Yr(channel,freq_index,:,index{i}));
                                mytitle = sprintf('Original ERSP at \n channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                data = squeeze(Yr(channel,:,index{i}));
                            end
                            average(i,:) = mean(data,2);
                            se = (std(data') ./ sqrt(numel(index{i})));
                            ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Yr,2));
                        end
                        
                        if strcmp(LIMO.Analysis,'Time')
                            mytitle = sprintf('Original ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmp(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Original Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            if strcmp(LIMO.Analysis,'Time-Frequency')
                                continuous(i,:,:) = Yr(channel,freq_index,:,sorting_values);
                                mytitle{i} = sprintf('Original single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            else
                                continuous(i,:,:) = Yr(channel,:,sorting_values);
                                mytitle{i} = sprintf('Original single trials \n sorted by regressor %g \n channel %s (%g) at %s Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            end
                        end
                    end
                    clear Yr
                elseif strcmp(extra,'Modelled')
                    Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                    Betas = Betas.Betas;
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        Betas = squeeze(Betas(channel,freq_index,:,:));
                    else
                        Betas = squeeze(Betas(channel,:,:));
                    end
                    Yh = (LIMO.design.X*Betas')'; % modelled data
                    
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            Yr = squeeze(Yr(:,freq_index,:,:));
                        end
                        R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                        
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            data = squeeze(Yh(:,index{i}));
                            average(i,:) = mean(data,2);
                            var   = diag(((R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')) / LIMO.model.model_df(2));
                            CI = sqrt(var/size(index{i},1))*z';
                            ci(i,:,:) = (repmat(mean(data,2),1,2)+CI)';
                        end
                        
                        if strcmp(LIMO.Analysis,'Time')
                            mytitle = sprintf('Modelled ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmp(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Modelled Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Modelled ERSP \n channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:) = Yh(:,sorting_values);
                            
                            if strcmp(LIMO.Analysis,'Time-Frequency')
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g \n channel %s (%g) at %g Hz', regressor(i), LIMO.data.chanlocs(channel).labels, channel, frequency);
                            else
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                            end
                        end
                    end
                else % Adjusted
                    allvar = 1:size(LIMO.design.X,2)-1;
                    allvar(regressor)=[];
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        Yr = load(fullfile(LIMO.dir,'Yr.mat')); 
                        Yr = squeeze(Yr.Yr(channel,freq_index,:,:));
                        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                        Betas = squeeze(Betas.Betas(channel,freq_index,:,:));
                    else
                        Yr = load(fullfile(LIMO.dir,'Yr.mat'));
                        Yr = squeeze(Yr.Yr(channel,:,:));
                        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
                        Betas = squeeze(Betas.Betas(channel,:,:));
                    end
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                    Ya = Yr - confounds; clear Yr Betas confounds;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            data = squeeze(Ya(:,index{i}));
                            average(i,:) = mean(data,2);
                            se = std(data') ./ sqrt(numel(index{i}));
                            ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Ya,1));
                        end
                        if strcmp(LIMO.Analysis,'Time')
                            mytitle = sprintf('Adjusted ERP at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        elseif strcmp(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Adjusted Power Spectrum at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle = sprintf('Adjusted ERSP channel %s (%g) at %g Hz', LIMO.data.chanlocs(channel).labels, channel, frequency);
                        end
                    else % continuous variable
                        for i=length(regressor):-1:1
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:) = Ya(:,sorting_values);
                            if strcmp(LIMO.Analysis,'Time-Frequency')
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
                    
                    if strcmp(LIMO.Analysis,'Frequency')
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
                    if strcmp(LIMO.Analysis,'Frequency')
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
                                if strcmp(LIMO.cache.fig.name,name) && ...
                                        LIMO.cache.fig.MCC == MCC && ...
                                        LIMO.cache.fig.threshold == p
                                    if strcmp(LIMO.Analysis,'Time-Frequency')
                                        sig = single(LIMO.cache.fig.mask(channel,freq_index,:)); sig(find(sig==0)) = NaN;
                                    else
                                        sig = single(LIMO.cache.fig.mask(channel,:)); sig(find(sig==0)) = NaN;
                                    end
                                end
                            else
                                if strcmp(LIMO.Analysis,'Time-Frequency')
                                    [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                    sig = single(squeeze(mask(channel,freq_index,:))); sig(find(sig==0)) = NaN;
                                else
                                    [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                    sig = single(mask(channel,:)); sig(find(sig==0)) = NaN;
                                end
                            end
                            h = axis;
                            if strcmp(LIMO.Analysis,'Frequency')
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
                                    if strcmp(LIMO.cache.fig.name,name) && ...
                                            LIMO.cache.fig.MCC == MCC && ...
                                            LIMO.cache.fig.threshold == p
                                        if strcmp(LIMO.Analysis,'Time-Frequency')
                                            sig = single(squeeze(LIMO.cache.fig.mask(channel,freq_index,:))); sig(find(sig==0)) = NaN;
                                        else
                                            sig = single(LIMO.cache.fig.mask(channel,:)); sig(find(sig==0)) = NaN;
                                        end
                                    end
                                else
                                    if strcmp(LIMO.Analysis,'Time-Frequency')
                                        [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                        sig = single(squeeze(mask(channel,freq_index,:))); sig(find(sig==0)) = NaN;
                                    else
                                        [~, mask] = limo_stat_values(name,p,MCC,LIMO,choice);
                                        sig = single(mask(channel,:)); sig(find(sig==0)) = NaN;
                                    end
                                end
                                h = axis;
                                if strcmp(LIMO.Analysis,'Frequency')
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
                if strcmp(LIMO.Analysis,'Frequency')
                    xlabel('Freq in Hz','FontSize',16)
                    ylabel('Power Spectrum in {\mu}V^2/Hz','FontSize',16);
                else
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude in {\mu}V','FontSize',16)
                end
                
                LIMO.cache.ERPplot.extra     = extra;
                LIMO.cache.ERPplot.average   = average;
                LIMO.cache.ERPplot.channel   = channel;
                LIMO.cache.ERPplot.regressor = regressor;
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.ERPplot.frequency = frequency;
                end
                LIMO.cache.ERPplot.ci        = ci;
                LIMO.cache.ERPplot.title     = mytitle;
                save LIMO LIMO
                
            else
                for i=1:size(continuous,1)
                    if i > 1; figure;set(gcf,'Color','w'); end
                    index = find(~isnan(squeeze(continuous(i,1,:))));
                    
                    if strcmp(LIMO.Analysis,'Frequency')
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
                    try
                        set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                    end
                end
                
                LIMO.cache.ERPplot.continuous = continuous;
                LIMO.cache.ERPplot.channel    = channel;
                LIMO.cache.ERPplot.regressor  = regressor;
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.ERPplot.frequency = frequency;
                end
                LIMO.cache.ERPplot.title      = mytitle;
                save LIMO LIMO
            end
            
    end % closes switch
    
    
    
    
    
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % -------------------                               LEVEL 2
    % -------------------                            GROUP EFFECTS
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
elseif LIMO.Level == 2
    
    if ~contains(FileName,'LIMO') % in all cases but course plot
        
        % if previously plotted, recover from the cache
        data_cached = 0;
        if isfield(LIMO,'cache') 
            try
                if strcmp(LIMO.cache.fig.name, FileName) && ...
                        LIMO.cache.fig.MCC == MCC && ...
                        LIMO.cache.fig.threshold == p
                    
                    disp('using cached data');
                    mask = LIMO.cache.fig.mask;
                    if isempty(mask)
                        data_cached = 0;
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        return
                    else
                        toplot      = LIMO.cache.fig.stats;
                        M           = LIMO.cache.fig.pval;
                        mask        = LIMO.cache.fig.mask;
                        mytitle     = LIMO.cache.fig.title;
                        data_cached = 1;
                        assignin('base','stat_values',M)
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                    end
                end
            catch no_cache
                data_cached = 0;
            end
        end
        
        % if there is no cached data, compute and plot
        % -------------------------------------------
        if data_cached == 0 
            
            [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,choice,[]);
            
            if isempty(mask)
                return
            elseif sum(mask(:)) == 0
                warndlg('  no values under threshold  ','no significant effect','modal');
                return
            else
                assignin('base','p_values',squeeze(M))
                assignin('base','mask',squeeze(mask))
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                if contains(FileName,'R2')
                    toplot = squeeze(toplot(:,:,:,1));
                elseif contains(FileName,'ttest','IgnoreCase',true)
                    toplot = squeeze(toplot(:,:,:,4));
                elseif strncmp(FileName,'con_',4)
                    toplot = squeeze(toplot(:,:,:,4));
                elseif strncmp(FileName,'ess_',4)
                    if ~exist('ess','var')
                        effect_nb = eval(FileName(22:end-4));
                        try
                            ess = eval(['ess' num2str(effect_nb)]);
                            clear(['ess' num2str(effect_nb)])
                          catch
                            ess = ess1; clear ess1
                        end
                    end
                    toplot = squeeze(toplot(:,:,:,4));
                elseif contains(FileName,'Condition') || ...
                        contains(FileName,'Covariate') || ...
                        contains(FileName,'Rep_ANOVA')
                    toplot = squeeze(toplot(:,:,:,1));
                else
                    disp('file no supported'); return
                end
            else
                if contains(FileName,'R2')
                    toplot = squeeze(toplot(:,:,1));
                elseif contains(FileName,'ttest','IgnoreCase',true)
                    toplot = squeeze(toplot(:,:,4));
                elseif strncmp(FileName,'con_',4)
                    toplot = squeeze(toplot(:,:,4));
                elseif strncmp(FileName,'ess_',4)
                    toplot = squeeze(toplot(:,:,4));
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
            save LIMO LIMO
        end
        
        % image all results
        % ------------------
        if Type == 1 && ~strcmp(LIMO.Analysis,'Time-Frequency') && ~strcmp(LIMO.Analysis,'ITC')
            limo_display_image(LIMO,toplot,mask,mytitle,flag)
                        
        elseif Type == 1 && strcmp(LIMO.Analysis,'Time-Frequency') || ...
                Type == 1 && strcmp(LIMO.Analysis,'ITC')
            if ndims(toplot)==3
                limo_display_image_tf(LIMO,toplot,mask,mytitle);
            else
                limo_display_image(LIMO,squeeze(toplot),squeeze(mask),mytitle)
            end
            
        elseif Type == 2
            %--------------------------
            % topoplot
            %--------------------------
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
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
                warndlg('no values under threshold','no significant effect');
            else
                EEG.data     = toplot;
                EEG.setname  = mytitle;
                EEG.pnts     = size(EEG.data,2);
                EEG.trials   = 1;
                EEG.chanlocs = LIMO.data.chanlocs;
                EEG.nbchan   = size(EEG.data,1);
                
                if size(toplot,2) == 1
                    opt = {'maplimits','maxmin','verbose','off'};
                    if isfield(LIMO,'Type')
                        if strcmp(LIMO.Type,'Components')
                            opt = {'maplimits','absmax','electrodes','off','verbose','off'};
                        end
                    end
                    figure; set(gcf,'Color','w','InvertHardCopy','off');
                    topoplot(toplot(:,1),EEG.chanlocs,opt{:});
                    title('Topoplot','FontSize',12)
                else
                    if strcmp(LIMO.Analysis,'Time')
                        EEG.xmin  = LIMO.data.start/1000; % in sec
                        EEG.xmax  = LIMO.data.end/1000;   % in sec
                        EEG.times = (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
                        pop_topoplot(EEG);
                    elseif strcmp(LIMO.Analysis,'Frequency')
                        EEG.xmin  = LIMO.data.freqlist(1);
                        EEG.xmax  = LIMO.data.freqlist(end);
                        freqlist  = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
                        if isempty(freqlist)
                            return
                        else
                            try
                                EEG.freq = eval(cell2mat(freqlist));
                            catch NO_NUM
                                EEG.freq = str2num(cell2mat(freqlist));
                            end
                            
                            if min(EEG.freq)<EEG.xmin || max(EEG.freq)>EEG.xmax
                                errordlg('selected frequency out of bound'); return
                            end
                        end
                        
                        N=length(EEG.freq);
                        figure;
                        for f=1:N
                            if N<=6
                                subplot(1,N,f)
                            else
                                subplot(ceil(N/6),6,f);
                            end
                            [~,ind] = min(abs(LIMO.data.freqlist-EEG.freq(f)));
                            topoplot(EEG.data(:,ind),EEG.chanlocs);
                            title([num2str(LIMO.data.freqlist(ind)) ' Hz'],'FontSize',12)
                        end
                        EEG.times = LIMO.data.freqlist;
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
                contains(FileName,'paired_samples','IgnoreCase',true) || contains(FileName,'con','IgnoreCase',true) || ...
                contains(FileName,'ess','IgnoreCase',true)
            % ------------------------------------------------------------------------------------------------------------
            % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
            % H0 file dim = (electrodes,frames,[t, p],nboot)
            
            data = load(fullfile(PathName,FileNme));
            data = data.(cell2mat(fieldnames(data)));            
            if strcmpi(LIMO.Analysis,'Time-Frequency')
                [~,channel,freq,time] = limo_display_reducedim(squeeze(data(:,:,:,[4 5])),LIMO);
                data                  = squeeze(data(channel,freq,time,:,:)); % 2D
                sig                   = squeeze(single(mask(channel,freq,time))); %1D
            else
                [~,channel,freq,time] = limo_display_reducedim(squeeze(data(:,:,[4 5])),LIMO);
                data                  = squeeze(data(channel,time,:));
                sig                   = single(mask(channel,:));
            end
            sig(find(sig==0)) = NaN;
            
            % compute
            trimci      = NaN(size(data,1),3);
            trimci(:,2) = data(:,1);
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
                xvect = LIMO.data.timevect;
            elseif strcmpi(LIMO.Analysis,'Frequency')
                xvect=LIMO.data.freqlist;
            elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                if length(time) > 1
                    xvect = LIMO.data.tf_times;
                elseif length(freq) > 1
                    xvect = LIMO.data.tf_freqs;
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
            if strcmp(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
                xlabel('Time in ms','FontSize',14)
                ylabel('Amplitude (A.U.)','FontSize',14)
            elseif strcmp(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
                xlabel('Frequency in Hz','FontSize',14)
                ylabel('Spectral Power (A.U.)','FontSize',14)
            end
            if size(data,1)>1
                title(sprintf('%s \n %s %s (%g)',mytitle,LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel),'FontSize',16); drawnow;
            else
                title(sprintf('%s \n virtual %s',mytitle,LIMO.Type(1:end-1)),'FontSize',16); drawnow;
            end
            assignin('base','Plotted_data',trimci);
            
            
        elseif strncmp(LIMO.design.name,'regression analysis',19) || strncmp(LIMO.design.name,'ANOVA',5) || strncmp(LIMO.design.name,'ANCOVA',6)
            % --------------------------------------------------------------------------------------------------------------------------------
            
            % which variable(s) to plot
            % ----------------------
            if size(LIMO.design.X,2) > 2
                input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2));
                regressor = inputdlg(input_title,'Plotting option');
                if isempty(regressor); return; end
                try regressor = sort(eval(cell2mat(regressor)));
                    if max(regressor) > size(LIMO.design.X,2);
                        errordlg('invalid regressor number');
                    end
                catch ME
                    return
                end
            else
                regressor = 1;
            end
            
            categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
            if max(regressor) == size(LIMO.design.X,2); tmp = regressor(1:end-1); else tmp = regressor; end
            cat = sum(tmp<=categorical); cont = sum(tmp>categorical);
            if cat >=1 && cont >=1
                errordlg('you can''t plot categorical and continuous regressors together'); return
            end
            
            % load the effect
            % --------------
            name = sprintf('Covariate_effect_%g.mat',regressor);
            data = load(name); data = data.(cell2mat(fieldnames(data)));
            
            % which ERP to make
            % ------------------
            extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Adjusted','Adjusted');
            if isempty(extra)
                return;
            elseif strcmp(extra,'Original')
                if regressor == size(LIMO.design.X,2)
                    errordlg('you can''t plot adjusted mean for original data'); return
                end
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                if size(data,1) > 1
                    if isempty(g.channels)
                        channel = inputdlg('which channel to plot','Plotting option');
                    else
                        channel = g.channels;
                    end
                    if isempty(channel) || strcmp(cell2mat(channel),'')
                        tmp = squeeze(data(:,:,:,1));
                        [channel,freq_index,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                        clear tmp ; frequency = LIMO.data.tf_freqs(freq_index);
                    else
                        channel = eval(cell2mat(channel));
                        if length(channel) > 1
                            error('1 channel only can be plotted')
                        elseif channel > size(data,1)
                            error('channel number invalid')
                        end
                    end
                else
                    channel = 1;
                end
                
                frequency = inputdlg('which frequency to plot','Plotting option');
                if isempty(frequency) || strcmp(cell2mat(frequency),'')
                    if ~exist('freq_index','var')
                        [v,f] = max(squeeze(data(channel,:,:,1)));
                        [v,c]=max(v); freq_index = f(c);
                        frequency = LIMO.data.tf_freqs(freq_index);
                    end
                else
                    frequency = eval(cell2mat(frequency));
                    [~,freq_index] = min(abs(LIMO.data.tf_freqs-frequency));
                    if length(frequency) > 1
                        error('1 frequency only can be plotted')
                    elseif freq_index > size(data,2)
                        error('frequency number invalid')
                    end
                end
                
            else % Time or Freq
                if size(data,1) > 1
                    if isempty(g.channels)
                        channel = inputdlg('which channel to plot','Plotting option');
                    else
                        channel = g.channels;
                    end
                    if isempty(channel) || strcmp(cell2mat(channel),'')
                        [v,e] = max(data(:,:,1)); [v,c]=max(v); channel = e(c);
                    else
                        channel = eval(cell2mat(channel));
                        if length(channel) > 1
                            error('1 channel only can be plotted')
                        elseif channel > size(data,1)
                            error('channel number invalid')
                        end
                    end
                else
                    channel = 1;
                end
            end
            clear data
            
            % down to business
            % ----------------------
            probs = [p/2; 1-p/2];
            z = norminv(probs);
            Yr = load(fullfile(LIMO.dir,'Yr.mat'));
            Yr = Yr.Yr;
            
            if strcmp(extra,'Original')
                if sum(regressor <= categorical) == length(regressor) % for categorical variables
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            data = squeeze(Yr(channel,freq_index,:,index{i}));
                        else
                            data = squeeze(Yr(channel,:,index{i}));
                        end
                        average(i,:) = nanmean(data,2);
                        se = (nanstd(data') ./ sqrt(numel(index{i})));
                        ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Yr,2));
                    end
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Original subjects'' parameters at optimized channel');
                    end
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [~,sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        reg_values(i,:) = LIMO.data.Cont(sorting_values,i);
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            continuous(i,:,:) = Yr(channel,freq_index,:,sorting_values);
                        else
                            continuous(i,:,:) = Yr(channel,:,sorting_values);
                        end
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Original subjects'' parameters \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle{i} = sprintf('Original subjects'' parameters, \n sorted by regressor %g at optimized channel', regressor(i));
                        end
                    end
                end
            elseif strcmp(extra,'Modelled')
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    load Betas; Betas = squeeze(Betas(channel,freq_index,:,:));
                    R = eye(size(Yr,4)) - (LIMO.design.X*pinv(LIMO.design.X));
                else
                    load Betas; Betas = squeeze(Betas(channel,:,:));
                    R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                end
                Yh = (LIMO.design.X*Betas')'; % modelled data
                
                if sum(regressor <= categorical) == length(regressor) % for categorical variables
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        data = squeeze(Yh(:,index{i}));
                        average(i,:) = nanmean(data,2);
                        index{i} = index{i}(find(~isnan(squeeze(Yr(channel,1,index{i})))));
                        var   = diag(((R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(channel,:,index{i}))')) / LIMO.model.model_df(2));
                        CI = sqrt(var/size(index{i},1))*z';
                        ci(i,:,:) = (repmat(nanmean(data,2),1,2)+CI)';
                    end
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Modelled subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Modelled subjects'' parameters at optimized channel');
                    end
                    clear Yr
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [~,sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        reg_values(i,:) = LIMO.data.Cont(sorting_values);
                        continuous(i,:,:) = Yh(:,sorting_values);
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Modelled single subjects'' parameters \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle{i} = sprintf('Modelled single subjects'' parameters \n sorted by regressor %g at optimized channel', regressor(i));
                        end
                    end
                end
            else % Adjusted
                allvar = [1:size(LIMO.design.X,2)-1]; allvar(regressor)=[];
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    Yr = squeeze(Yr(channel,freq_index,:,:));
                    load Betas; Betas = squeeze(Betas(channel,freq_index,:,:));
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                else
                    Yr = squeeze(Yr(channel,:,:));
                    load Betas; Betas = squeeze(Betas(channel,:,:));
                    confounds = (LIMO.design.X(:,allvar)*Betas(:,allvar)')';
                end
                Ya = Yr - confounds; clear Yr Betas confounds;
                
                if sum(regressor <= categorical) == length(regressor) % for categorical variables
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        data = squeeze(Ya(:,index{i}));
                        average(i,:) = nanmean(data,2);
                        se = nanstd(data') ./ sqrt(numel(index{i}));
                        ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Ya,1));
                    end
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Adjusted subjects'' parameters at channel %s (%g)', LIMO.data.chanlocs(channel).labels, channel);
                    else
                        mytitle = sprintf('Adjusted subjects'' parameters at  at optimized channel');
                    end
                    clear Yr
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [~,sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        reg_values(i,:) = LIMO.data.Cont(sorting_values);
                        continuous(i,:,:) = Ya(:,sorting_values);
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Adjusted single subjects'' parameters \n sorted by regressor %g channel %s (%g)', regressor(i), LIMO.data.chanlocs(channel).labels, channel);
                        else
                            mytitle{i} = sprintf('Adjusted single subjects'' parameters \n sorted by regressor %g at optimized channel', regressor(i));
                        end
                    end
                end
            end
            
            % make the figure(s)
            % ------------------
            figure;set(gcf,'Color','w')
            if sum(regressor <= categorical) == length(regressor)
                for i=1:size(average,1)
                    if strcmp(LIMO.Analysis,'Time')
                        timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                        plot(timevect,squeeze(average(i,:)),'LineWidth',1.5); hold on
                        xlabel('Time in ms','FontSize',14)
                        ylabel('Amplitude (A.U.)','FontSize',14)
                    else
                        freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                        plot(freqvect,squeeze(average(i,:)),'LineWidth',1.5); hold on
                        xlabel('Frequency in Hz','FontSize',14)
                        ylabel('Spectral Power (A.U.)','FontSize',14)
                    end
                    
                    if i==1
                        colorOrder = get(gca, 'ColorOrder');
                        colorOrder = repmat(colorOrder,ceil(size(average,1)/size(colorOrder,1)),1);
                    end
                    x = squeeze(ci(i,1,:)); y = squeeze(ci(i,2,:));
                    fillhandle = patch([timevect fliplr(timevect)], [x',fliplr(y')], colorOrder(i,:));
                    set(fillhandle,'EdgeColor',colorOrder(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                end
                
                % if regressor spans columns of an effect, plot significant time frames
                index = 1; index2 = LIMO.design.nb_conditions(1);
                for i=1:length(LIMO.design.nb_conditions)
                    effect = index:index2;
                    if length(regressor) == length(effect)
                        if mean(regressor == effect) == 1
                            name = sprintf('Condition_effect_%g.mat',i);
                            load(name); [M, mask, mytitle2] = limo_stat_values(name,p,MCC,LIMO,choice);
                            sig = single(mask(channel,:)); sig(find(sig==0)) = NaN;
                            h = axis;
                            if strcmp(LIMO.Analysis,'Time')
                                plot(timevect,(sig./10+1).*h(3),'r*','LineWidth',2)
                            else
                                plot(freqvect,(sig./10+1).*h(3),'r*','LineWidth',2)
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
                                load(name); [M, mask, mytitle2] = limo_stat_values(name,p,MCC,LIMO,choice);
                                sig = single(mask(channel,:)); sig(find(sig==0)) = NaN;
                                h = axis;
                                if LIMO.analysis_flag == 1
                                    plot(timevect,(sig./10+1).*h(3),'r*','LineWidth',2)
                                else
                                    plot(freqvect,(sig./10+1).*h(3),'r*','LineWidth',2)
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
                title(mytitle,'FontSize',16); drawnow;
                assignin('base','Plotted_data', average)
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                if LIMO.analysis_flag == 1
                    ylabel('Amplitude (A.U.)','FontSize',16)
                    xlabel('Time in ms','FontSize',16)
                else
                    ylabel('Spectral Power (A.U.)','FontSize',16)
                    xlabel('Frequency in Hz','FontSize',16)
                end
                
            else
                for i=1:size(continuous,1)
                    if i > 1; figure;set(gcf,'Color','w'); end
                    index = find(~isnan(squeeze(continuous(i,1,:))));
                    if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'Time-Frequency')
                        if strcmp(LIMO.Analysis,'Time')
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
                    axis tight; title(mytitle{i},'FontSize',14); drawnow;
                    xlabel('Sorted subjects','FontSize',16)
                    try
                        set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                    end
                end
            end
            
            
            
        elseif contains(FileName,'Rep_ANOVA')   % All stuffs for repeated measures ANOVA
            % -----------------------------------------------------------------------------
            
            
            if contains(FileName,'Rep_ANOVA_Main')
                % -----------------------------------
                
                % which summary stat
                % -------------------
                extra = questdlg('Summarize data using:','Data plot option','Mean','Trimmed','Trimmean');
                if isempty(extra)
                    return
                end
                % -----------------------
                
                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                if contains(FileName,'Main_effect','IgnoreCase',true)
                    index1 = strfind(FileName,'Main_effect')+12;
                elseif contains(FileName,'Interaction','IgnoreCase',true)
                    index1 = strfind(FileName,'Interaction')+12;
                end
                index2                = max(strfind(FileName,'_'))-1;
                effect_nb             = eval(FileName(index1:index2));
                C                     = LIMO.design.C{effect_nb};
                Data                  = load(fullfile(LIMO.dir,'Yr.mat'));
                Data                  = Data.(cell2mat(fieldnames(Data)));
                [~,channel,freq,time] = limo_display_reducedim(Rep_ANOVA,LIMO);
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    Data              = squeeze(Data(channel,freq,time,:,:));
                    sig               = squeeze(single(mask(channel,freq,time))); 
                else
                    Data              = squeeze(Data(channel,time,:,:)); % note freq/time variables have the same values
                    sig               = single(mask(channel,:)); 
                end
                sig(find(sig==0)) = NaN;
                
                % compute differences between pairs using C and Cov
                n = size(Data,2);
                if strcmp(extra,'Mean')
                    for time_or_freq=1:size(Data,1)
                        avg(time_or_freq,:) = nanmean(C*squeeze(Data(time_or_freq,:,:))',2);
                        S(time_or_freq,:,:) = nancov(squeeze(Data(time_or_freq,:,:)));
                    end
                    if size(Rep_ANOVA,1)>1
                        mytitle = sprintf('Original %s \n %s %s (%g)',mytitle,LIMO.Type(1:end-1),LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Original %s \n virtual %s',mytitle,LIMO.Type(1:end-1));
                    end
                else
                    g=floor((20/100)*n); %% compute for 20% trimmed mean
                    for time_or_freq=1:size(Data,1)
                        [v,indices] = sort(squeeze(Data(time_or_freq,:,:))); % sorted data
                        TD(time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
                        avg(time_or_freq,:) = nanmean(C*squeeze(TD(time_or_freq,:,:))',2);
                        v(1:g+1,:)=repmat(v(g+1,:),g+1,1);
                        v(n-g:end,:)=repmat(v(n-g,:),g+1,1); % winsorized data
                        [~,reorder] = sort(indices);
                        for j = 1:size(Data,3), SD(:,j) = v(reorder(:,j),j); end % restore the order of original data
                        S(time_or_freq,:,:) = cov(SD); % winsorized covariance
                    end
                    if size(Data,1)>1
                        mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
                    end
                end
                
                % CI
                df = rank(C);
                dfe = size(Data,2)-size(Data,3)+1;
                % c = avg + 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C'))); % uses Bonferoni inequality
                % b = avg - 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time_or_freq,:,:))*C')));
                bound = (abs(tinv(p./(2*size(C,1)),dfe)).*diag((sqrt(C*squeeze(S(time_or_freq,:,:))*C'))));
                c = avg + repmat(bound', [length(avg),1]);
                b = avg - repmat(bound', [length(avg),1]);

                % do the figure
                colours = limo_color_images(size(avg,2));
                if strcmpi(LIMO.Analysis,'Time')
                    xvect = LIMO.data.timevect;
                elseif strcmpi(LIMO.Analysis,'Frequency')    
                    xvect=LIMO.data.freqlist;
                elseif strcmpi(LIMO.Analysis,'Time-Frequency') 
                    if length(time) > 1
                        xvect = LIMO.data.tf_times;
                    elseif length(freq) > 1
                        xvect = LIMO.data.tf_freqs;
                    end
                end
                
                figure;set(gcf,'Color','w')
                for cond = 1:size(c,2)
                    plot(xvect,avg(:,cond)','LineWidth',3,'color',colours(cond,:));
                    fillhandle = patch([xvect fliplr(xvect)], [c(:,cond)',fliplr(b(:,cond)')], colours(cond,:));
                    set(fillhandle,'EdgeColor',colours(cond,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                    hold on
                end
                grid on; box on; axis tight
                h = axis;  hold on;
                plot(xvect,(sig./10+1).*h(3),'k.','MarkerSize',20)
                if strcmp(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(time) > 1
                    xlabel('Time in ms','FontSize',14)
                    ylabel('Amplitude (A.U.)','FontSize',14)
                elseif strcmp(LIMO.Analysis,'Frequency') || strcmpi(LIMO.Analysis,'Time-Frequency') && length(freq) > 1
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
                [e,f,d]=size(Rep_ANOVA_Gp_effect);
                
                % check channel to plot
                if e > 1
                    if isempty(g.channels)
                        channel = inputdlg('which channel to plot','Plotting option');
                    else
                        channel = g.channels;
                    end
                    if isempty(channel) || strcmp(cell2mat(channel),'')
                        clear channel; [v,e] = max(Rep_ANOVA_Gp_effect(:,:,1)); [v,c]=max(v); channel = e(c);
                    else
                        channel = eval(cell2mat(channel));
                        if length(channel) > 1
                            error('1 channel only can be plotted')
                        elseif channel > e
                            error('channel number invalid')
                        end
                    end
                else
                    if length(LIMO.design.electrode) == 1
                        channel = LIMO.design.electrode;
                    else
                        channel = 1;  % accomodates the fact that all matrices have the channel dim (even = 1)
                    end
                end
                
                % compute the pair-wise differences and plot
                load Yr;
                if strcmp(extra,'Original')
                    Data = mean(squeeze(Yr(channel,:,:,:)),3);
                    combinations = nchoosek(1:LIMO.design.nb_conditions,2);
                    for d = 1:size(combinations,1)
                        Effect(:,d) = nanmean(Data(:,find(LIMO.data.Cat == combinations(d,1))),2) - nanmean(Data(:,find(LIMO.data.Cat == combinations(d,2))),2);
                    end
                end
                
                % design matrix
                X = zeros(size(Yr,3),LIMO.design.nb_conditions+1);
                X(:,end) = 1;
                for i=1:LIMO.design.nb_conditions
                    X(find(LIMO.data.Cat == i),i) = 1;
                end
                
                % data again
                Y = nanmean(squeeze(Yr(channel,:,:,:)),3);
                X = X(find(~isnan(Y(1,:))),:);
                Y = Y(:,find(~isnan(Y(1,:))))';
                if strcmp(extra,'Modelled')
                    beta = pinv(X)*Y; Yhat = X*beta;
                    combinations = nchoosek(1:LIMO.design.nb_conditions,2);
                    for d = 1:size(combinations,1)
                        Effect(:,d) = mean(Yhat(:,find(X(:,combinations(d,1)))),2) - nanmean(Yhat(:,find(X(:,combinations(d,2)))),2);
                    end
                end
                
                Res = (Y'*[eye(size(Y,1)) - (X*pinv(X))]*Y);
                df = size(Y,1)-rank(X); t = tcdf(1-p,df);
                sigma2 = sum((Res.^2./df),2);
                v = t.*sqrt(sigma2 ./ norm(X(:,1:end-1)).^2);
                b = Effect - v(channel,:);
                c = Effect + v(channel,:);
                if channel == 1
                    if length(LIMO.design.electrode) == 1
                        if strcmp(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        end
                    else
                        if strcmp(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n optimized channel');
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n optimized channel');
                        end
                    end
                else
                    if strcmp(extra,'Original')
                        mytitle = sprintf('Mean parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Modelled parameter difference between groups \n at channel %s (%g)', LIMO.data.chanlocs(channel).labels,channel);
                    end
                end
                
                figure;set(gcf,'Color','w')
                if strcmp(LIMO.Analysis,'Time')
                    timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
                    plot(timevect,Effect,'LineWidth',3);
                    fillhandle = patch([timevect fliplr(timevect)], [c',fliplr(b')], [1 0 0]);
                else
                    freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                    plot(freqvect,Effect,'LineWidth',3);
                    fillhandle = patch([freqvectt fliplr(freqvect)], [c',fliplr(b')], [1 0 0]);
                end
                set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                grid on; box on; axis tight
                sig = single(mask(channel,:)); sig(find(sig==0)) = NaN;
                h = axis;  hold on;
                if strcmp(LIMO.Analysis, 'Time');
                    plot(timevect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                    xlabel('Time in ms','FontSize',14)
                    ylabel('Amplitude (A.U.)','FontSize',14)
                else
                    plot(freqvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                    xlabel('Frequency in Hz','FontSize',14)
                    ylabel('Spectral Power (A.U.)','FontSize',14)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                Gp_difference = [c' ; Effect' ; b'];
                assignin('base','Plotted_data',Gp_difference);
                
                
                
                % -------------------------
            elseif contains(FileName,'Rep_ANOVA_Interaction') % Gp * Repeated measures - plot differences btween condition per gp
                % ------------------------
                
                
                % which ERP to make
                % ------------------
                extra = questdlg('Plotting ERP','ERP Options','Original','Trimmed','Original');
                if isempty(extra)
                    return;
                end
                % -----------------------
                
                [e,f,d]=size(Rep_ANOVA_Interaction_with_gp);
                
                if e > 1
                    if isempty(g.channels)
                        channel = inputdlg('which channel to plot','Plotting option');
                    else
                        channel = g.channels;
                    end
                    if isempty(channel) || strcmp(cell2mat(channel),'')
                        [v,e] = max(Rep_ANOVA_Interaction_with_gp(:,:,1)); [v,c]=max(v); channel = e(c);
                    else
                        channel = eval(cell2mat(channel));
                        if length(channel) > 1
                            error('1 channel only can be plotted')
                        elseif channel > e
                            error('channel number invalid')
                        end
                    end
                else
                    channel = 1;
                end
                
                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                effect_nb = eval(FileName(end-4));
                C = LIMO.design.C{effect_nb};
                load Yr; Data = squeeze(Yr(channel,:,:,:));
                
                % compute differences between pairs using C and Cov
                n = size(Data,2);
                if strcmp(extra,'Original')
                    for gp = 1:LIMO.design.nb_conditions
                        index = find(LIMO.data.Cat==gp);
                        for time=1:size(Data,1)
                            avg(gp,time,:) = nanmean(C*squeeze(Data(time,index,:))',2);
                            S(gp,time,:,:) = cov(squeeze(Data(time,index,:)));
                        end
                    end
                    if e>1
                        mytitle = sprintf('Original %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Original %s \n optimized channel',mytitle);
                    end
                else
                    for gp = 1:LIMO.design.nb_conditions
                        index = find(LIMO.data.Cat==gp);
                        n = length(index);
                        g=floor((20/100)*n); %% compute for 20% trimmed mean
                        for time_or_freq=1:size(Data,1)
                            [v,indices] = sort(squeeze(Data(time_or_freq,index,:))); % sorted data
                            TD(gp,time_or_freq,:,:) = v((g+1):(n-g),:); % trimmed data - doesn't matter if relationship was kept since we only compute means
                            avg(gp,time_or_freq) = nanmean(C*squeeze(TD(gp,time_or_freq,:,:))',2);
                            v(1:g+1,:)=repmat(v(g+1,:),g+1,1);
                            v(n-g:end,:)=repmat(v(n-g,:),g+1,1); % winsorized data
                            [~,reorder] = sort(indices);
                            for j = 1:size(Data,3), SD(:,j) = v(reorder(:,j),j); end % restore the order of original data
                            S(gp,time_or_freq,:,:) = cov(SD); % winsorized covariance
                        end
                        clear SD
                    end
                    if e>1
                        mytitle = sprintf('Trimmed %s \n channel %s (%g)',mytitle,LIMO.data.chanlocs(channel).labels,channel);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized channel',mytitle);
                    end
                end
                
                
                figure;set(gcf,'Color','w')
                df = rank(C);
                timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
                for gp = 1:LIMO.design.nb_conditions
                    index = find(LIMO.data.Cat==gp);
                    dfe = size(Data,2)-length(index)+1;
                    if gp==1
                        if strcmp(LIMO.Analysis,'Time')
                            plot(timevect,avg(gp,:),'LineWidth',3);
                        elseif LIMO.analysis_flag == 2
                            plot(freqvect,avg(gp,:),'LineWidth',3);
                        end
                        colorOrder = get(gca, 'ColorOrder');
                    else
                        if strcmp(LIMO.Analysis,'Time')
                            plot(timevect,avg(gp,:),'Color',colorOrder(gp,:),'LineWidth',3);
                        else
                            plot(freqvect,avg(gp,:),'Color',colorOrder(gp,:),'LineWidth',3);
                        end
                    end

                    % there is an error in the following 2 formulas I cannot fix -- @disbeat
                    c = avg(gp,:,:) + tinv(p./(2*size(C,1)),dfe).*(sqrt(C*squeeze(S(gp,:,:,:))*C'));
                    b = avg(gp,:,:) - tinv(p./(2*size(C,1)),dfe).*(sqrt(C*squeeze(S(gp,:,:,:))*C'));
                    if strcmp(LIMO.Analysis,'Time')
                        fillhandle = patch([timevect fliplr(timevect)], [c,fliplr(b)], colorOrder(gp,:));
                        xlabel('Time in ms','FontSize',14)
                        ylabel('Amplitude (A.U.)','FontSize',14)
                    else
                        fillhandle = patch([freqvect fliplr(freqvect)], [c,fliplr(b)], colorOrder(gp,:));
                        xlabel('Frequency in Hz','FontSize',14)
                        ylabel('Spectral Power (A.U.)','FontSize',14)
                    end
                    set(fillhandle,'EdgeColor',colorOrder(gp,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                    hold on
                end
                
                grid on; box on; axis tight
                sig = single(mask(channel,:)); sig(find(sig==0)) = NaN;
                h = axis;  hold on;
                if strcmp(LIMO.Analysis,'Time')
                    plot(timevect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                else
                    plot(freqvect,(sig./10+1).*h(3),'r.','MarkerSize',20)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                trimci = [c ; avg ; b];
                assignin('base','Plotted_data',trimci);
                
            end
        else
            errordlg('this file is not supported for this kind of plot','Nothing plotted')
        end
    end % closes type
    
elseif strcmp(LIMO.Level,'LI')
    
    [M, mask, mytitle] = limo_stat_values(FileName,p,MCC,LIMO,[],[]);
    
    if Type == 1
        %--------------------------
        % imagesc of the results
        %--------------------------
        if sum(mask(:)) == 0
            warndlg('no values under threshold parameter','no significant effect');
        else
            scale = M.*mask;
            v = max(scale(:)); [e,f]=find(scale==v);
            if min(scale(:))<0
                scale(scale==0)=min(scale(:))+(min(scale(:))/10);
            else
                scale(scale==0)=NaN;
            end
            
            figure; set(gcf,'Color','w');
            timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(M,2));
            imagesc(timevect,1:size(M,1),scale);
            title([mytitle],'FontSize',18);
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
try
    set(gca,'YTick',1:length(LIMO.data.expected_chanlocs));
catch ME
    set(gca,'YTick',1:length(LIMO.data.chanlocs));
end

ylabel('Electrodes','FontSize',14);
if strcmp(LIMO.Analysis,'Time')
    xlabel('Time in ms','FontSize',16)
elseif strcmp(LIMO.Analysis,'Frequency')
    xlabel('Frequency in Hz','FontSize',16)
end

if LIMO.Level == 1
    if strcmp(LIMO.data.chanlocs,'Components')
        label_electrodes = [];
    else
        for i = 1 : length(LIMO.data.chanlocs)
            try
                label_electrodes{i} = LIMO.data.expected_chanlocs(i).labels;
            catch ME
                label_electrodes{i} = LIMO.data.chanlocs(i).labels;
            end
        end
    end
else
    if isempty(LIMO.design.electrode)
        for i = 1 : length(LIMO.data.chanlocs)
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



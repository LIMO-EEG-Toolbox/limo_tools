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
% Cyril Pernet, Guillaume Rousselet v3 06-05-2009
% Carl Gaspar 03-09-2009 - fixed some axis issues for 3D plots (see subfunction time_vect_)
% Cyril P. v4 09-09-2009 allows random effect results to be displayed (+ some clean up)
% Cyril P. v5. 10-03-2010 split the whole file into 2 parts based on LIMO.level (1 or 2)
% Guillaume Rousselet v4 06-07-2010 added the max(T)/max(F) and cluster stats for random effect
% Cyril Pernet v4 16-05-2010 fixed the random effect to automatically load bootstrap and get the neighbouring matrix for clusters
% Nicolas Chauveau 08-12-2011 fixed the ERP plot of gp*repeated measures (for levels>2)
% Cyril Pernet v5 10-10-2012 added tfce and redesigned CI with filling
% Andrew Stewart 10-11-2013 added options for spectral power and time-freq
% Cyril Pernet 21-03-2014 made time-freq to work with the new display +
% changed limo_stat values to take timne-freq
% Cyril Pernet & Ramon Martinez-Cancino 23-10-2014 updates for components (ICA)
%
% see also limo_stat_values topoplot
% ----------------------------------
%  Copyright (C) LIMO Team 2010

try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g = []; end;
catch
    disp('limo_display_results() error: calling convention {''key'', value, ... } error'); return;
end;

try g.channels;   catch, g.channels  = [];  end; % No default values
try g.regressor;  catch, g.regressor = [];  end; % No default values
try g.plot3type;  catch, g.plot3type = [];  end; % No default values

cd(PathName)
load (FileName);
if nargin <= 6
    flag = 1;
end

choice = 'use theoretical p values'; % threshold based on what is computed since H0 is used for clustering
% see limo_stat_values

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
                                toplot = []; return
                            else
                                M = LIMO.cache.fig.pval;
                                mytitle = LIMO.cache.fig.title;
                                toplot = LIMO.cache.fig.stats;
                                data_cached = 1;
                                assignin('base','p_values',M)
                                assignin('base','mask',mask)
                            end
                        end
                    catch no_cache
                        data_cached = 0;
                    end
                end
                
                % ------------------
                % compute the plot
                % ------------------
                if data_cached == 0;
                    
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        [M, mask, mytitle] = limo_stat_values_tf(Type,FileName,p,MCC,LIMO,choice);
                    else
                        [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,choice);
                    end
                    
                    if isempty(mask)
                        return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        LIMO.cache.fig.name       = FileName;
                        LIMO.cache.fig.MCC        = MCC;
                        LIMO.cache.fig.stats      = [];
                        LIMO.cache.fig.threshold  = p;
                        LIMO.cache.fig.pval       = M;
                        LIMO.cache.fig.mask       = mask;
                        LIMO.cache.fig.title      = mytitle;
                        save LIMO LIMO
                        toplot = []; return
                    else
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                    end
                    
                    if strcmp(FileName,'R2.mat')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(R2(:,:,:,1)); % plot R2 values instead of F
                        else
                            toplot = squeeze(R2(:,:,1));
                        end
                        assignin('base','R2_values',squeeze(R2(:,:,2)))
                        clear R2
                        
                    elseif strncmp(FileName,'Condition_effect',16)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(Condition_effect(:,:,:,1)); % plot F values
                        else
                            toplot = squeeze(Condition_effect(:,:,1));
                        end
                        assignin('base','F_values',toplot)
                        clear Condition_effect
                        
                    elseif strncmp(FileName,'Covariate_effect',16)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(Covariate_effect(:,:,:,1)); % plot F values
                        else
                            toplot = squeeze(Covariate_effect(:,:,1));
                        end
                        assignin('base','F_values',toplot)
                        clear Covariate_effect
                        
                    elseif strncmp(FileName,'Interaction_effect',18)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(Interaction_effect(:,:,:,1)); % plot F values
                        else
                            toplot = squeeze(Interaction_effect(:,:,1));
                        end
                        assignin('base','F_values',toplot)
                        clear Interaction_effect
                        
                    elseif strncmp(FileName,'semi_partial_coef',17)
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(semi_partial_coef(:,:,:,1)); % plot the coef rather than F
                        else
                            toplot = squeeze(semi_partial_coef(:,:,1)); % plot the coef rather than F
                        end
                        assignin('base','semi_partial_coef',toplot)
                        clear semi_partial_coef
                        
                    elseif strcmp(FileName(1:4),'con_')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(con(:,:,:,4)); % plot T values
                        else
                            toplot = squeeze(con(:,:,4));
                        end
                        assignin('base','T_values',toplot)
                        clear con
                        
                    elseif strcmp(FileName(1:4),'ess_')
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            toplot = squeeze(ess(:,:,:,end-1)); % plot F values
                        else
                            toplot = squeeze(ess(:,:,end-1));
                        end
                        assignin('base','F_values',toplot)
                        clear ess
                        
                    else
                        errordlg2('file not supported');
                        return
                    end
                    
                    update_cache = 1;
                end
                
% -------------------------------------------------------------------------
%              Actual plot takes place here
% -------------------------------------------------------------------------
                if ~isempty(toplot)
                    
                    % cache the results for next time
                    if data_cached == 0
                        LIMO.cache.fig.name       = FileName;
                        LIMO.cache.fig.MCC        = MCC;
                        LIMO.cache.fig.stats      = toplot;
                        LIMO.cache.fig.threshold  = p;
                        LIMO.cache.fig.pval       = M;
                        LIMO.cache.fig.mask       = mask;
                        LIMO.cache.fig.title      = mytitle;
                        save LIMO LIMO
                    end
                    
                    
                    if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'Frequency')
                        limo_display_image(LIMO,toplot,mask,mytitle)
                                               
                    else % strcmp(LIMO.Analysis,'Time-Frequency')  - 3D maps
                        limo_display_results_tf(LIMO,toplot,mask,mytitle);
                    end
                end
                
            else
                
                % mutivariate results from 1st level analysis
                % ------------------------------------------
                
                if strncmp(FileName,'R2',2)
                    
                    cd(LIMO.dir); load R2_EV.mat; EV = R2_EV(1:5,:); % no point plotting 0, just pick 5 1st Eigen values
                    test = sum(sum(R2_EV(1:5,:)>1,2)>1); % check if we have more than 1 EV>1
                    if test ==1; choice = 'Roy'; else choice = 'Pillai'; end
                    clear R2_EV;
                    
                    load R2.mat;
                    F_values(:,1) = squeeze(R2(:,2));
                    F_values(:,2) = squeeze(R2(:,4));
                    [M, mask, mytitle] = limo_statm_values(Type,FileName,p,MCC,LIMO,choice);
                    if isempty(mask)
                        return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        toplot = []; return
                    else
                        toplot = squeeze(R2(:,1)); % plot R2 values instead of F
                        assignin('base','F_values',F_values)
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                        clear R2
                    end
                    
                elseif strncmp(FileName,'Condition_effect',16)
                    
                    cd(LIMO.dir);
                    if strcmp(FileName(end-6:end),'_EV.mat'); FileName = [FileName(1:end-7) '.mat']; end
                    name = sprintf('Condition_effect_%g_EV',eval(FileName(18:end-4))); load(name);
                    EV = Condition_effect_EV(1:5,:); % no point plotting 0, just pick 5 1st Eigen values
                    test =  sum(sum(EV>1,2)>1); % check if we have more than 1 EV>1
                    if test <=1; choice = 'Roy'; else choice = 'Pillai'; end
                    clear Condition_effect_EV;
                    
                    load(FileName);
                    F_values(:,1) = squeeze(Condition_effect(:,1));
                    F_values(:,2) = squeeze(Condition_effect(:,3));
                    [M, mask, mytitle] = limo_statm_values(Type,FileName,p,MCC,LIMO,choice);
                    if isempty(mask)
                        return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        toplot = []; return
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
                    
                elseif strncmp(FileName,'Covariate_effect',16)
                    
                    cd(LIMO.dir);
                    if strcmp(FileName(end-6:end),'_EV.mat'); FileName = [FileName(1:end-7) '.mat']; end
                    name = sprintf('Covariate_effect_%g_EV',eval( FileName(18:end-4))); load(name);
                    EV = Covariate_effect_EV(1:5,:); % no point plotting 0, just pick 5 1st Eigen values
                    test = sum(sum(Covariate_effect_EV(1:5,:)>1,2)>1); % check if we have more than 1 EV>1
                    if test ==1; choice = 'Roy'; else choice = 'Pillai'; end
                    clear Covariate_effect_EV;
                    
                    load(FileName);
                    F_values(:,1) = squeeze(Covariate_effect(:,1));
                    F_values(:,2) = squeeze(Covariate_effect(:,3));
                    [M, mask, mytitle] = limo_statm_values(Type,FileName,p,MCC,LIMO,choice);
                    if isempty(mask)
                        return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        toplot = []; return
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
                h = subplot(3,3,[4 5 7 8]);
                timevect = linspace(LIMO.data.start,LIMO.data.end,size(EV,2));
                ratio = (LIMO.data.end*1000 - LIMO.data.start*1000) / size(EV,2);
                if LIMO.data.start < 0
                    frame_zeros = round(abs(LIMO.data.start*1000) / ratio);
                end
                scale = EV; scale(scale==0)=NaN;
                imagesc(timevect,1:size(EV,1),scale);
                color_images_(scale,LIMO);  colorbar
                ylabel('Eigen Values','Fontsize',14)
                set(gca,'YTickLabel',{'1','2','3','4','5'});
                title('5 first Eigen values','Fontsize',14)
                
                % imagesc effect values
                subplot(3,3,[1 2]);
                scale = toplot'.*mask; scale(scale==0)=NaN;
                imagesc(timevect,1,scale);
                v = max(toplot(:)); [~,f]=find(toplot==v);
                try
                    caxis([min(min(scale)), max(max(scale))]);
                end
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
                else strcmp(LIMO.Analysis,'Frequency')
                    EEG.xmin = LIMO.data.freqlist(1);
                    EEG.xmax = LIMO.data.freqlist(end);
                    freqlist = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
                    if isempty(freqlist)
                        return
                    else
                        try
                            EEG.freq = eval(cell2mat(freqlist));
                        catch NO_NUM
                            EEG.freq = str2num(cell2mat(freqlist));
                        end
                        
                        if min(EEG.freq)<EEG.xmin || max(EEG.freq)>EEG.xmax
                            errordlg('slected frequency out of bound'); return
                        end
                    end
                end
                
                if strcmp(FileName,'R2.mat')
                    EEG.data = squeeze(R2(:,:,2));
                    EEG.pnts = size(EEG.data,2);
                    EEG.nbchan = size(EEG.data,1);
                    EEG.setname = 'R2 - F values';
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
                elseif strncmp(FileName,'Condition_effect',16)
                    EEG.data = squeeze(Condition_effect(:,:,1));
                    EEG.pnts = size(EEG.data,2);
                    EEG.nbchan = size(EEG.data,1);
                    EEG.setname = sprintf('Condition %s - F values',FileName(18:end-4));
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
                    assignin('base',EEG.setname(1:9),EEG.data);
                elseif strncmp(FileName,'Covariate_effect',16)
                    EEG.data = squeeze(Covariate_effect(:,:,1));
                    EEG.pnts = size(EEG.data,2);
                    EEG.nbchan = size(EEG.data,1);
                    EEG.setname = sprintf('Covariate %s - F values',FileName(18:end-4));
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
                    assignin('base',EEG.setname(1:9),EEG.data);
                elseif strncmp(FileName,'Interaction_effect',18)
                    EEG.data = squeeze(Interaction_effect(:,:,1));
                    EEG.pnts = size(EEG.data,2);
                    EEG.nbchan = size(EEG.data,1);
                    EEG.setname = sprintf('Interaction %s - F values',FileName(20:end-4));
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
                    assignin('base',EEG.setname(1:11),EEG.data);
                elseif strcmp(FileName,'semi partial_coef.mat')
                    regressor = str2num(cell2mat(inputdlg('which regressor(s) to plot (e.g. 1:3)','Plotting option')));
                    if max(regressor) > size(Partial_coef,3); errordlg('error in regressor number'); return; end
                    for b = regressor
                        EEG.data = squeeze(Partial_coef(:,:,b,1));
                        EEG.pnts = size(EEG.data,2);
                        EEG.nbchan = size(EEG.data,1);
                        EEG.setname = sprintf('semi partial coef R2 values variable %g',b);
                        EEG.pnts  = size(EEG.data,2);
                        EEG.times = LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000;
                        EEG.trials = 1;
                        EEG.chanlocs = LIMO.data.chanlocs;
                        EEG.nbchan = size(EEG.data,1);
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
                    EEG.data = squeeze(con(:,:,4));
                    EEG.pnts = size(EEG.data,2);
                    EEG.nbchan = size(EEG.data,1);
                    EEG.setname = ['Contrast ',[FileName(5:end-4)],' -- T values'];
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
                    EEG.data = squeeze(ess(:,:,end-1));
                    EEG.pnts = size(EEG.data,2);
                    EEG.nbchan = size(EEG.data,1);
                    EEG.setname = ['Contrast ',[FileName(5:end-4)],' -- F values'];
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
            try regressor = sort(eval(cell2mat(regressor)));
                if max(regressor) > size(LIMO.design.X,2); errordlg('invalid regressor number'); end
            catch ME
                return
            end
            
            categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
            if max(regressor) == size(LIMO.design.X,2); tmp = regressor(1:end-1); else tmp = regressor; end
            cat = sum(tmp<=categorical); cont = sum(tmp>categorical);
            if cat >=1 && cont >=1
                errordlg('you can''t plot categorical and continuous regressors together'); return
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
            
            % which electrode/frequency to plot
            % --------------------------------
            if isempty(g.channels)
                electrode = inputdlg('which electrode to plot','Plotting option');
            else
                electrode = g.channels;
            end
            if strcmp(LIMO.Analysis,'Time-Frequency')
                disp('loading the 4D data ...')
                frequency = inputdlg('which Frequency to plot','Plotting option');
            else
                frequency = [];
            end
            
            if strcmp(electrode,'') || strcmp(frequency,'')
                disp('looking for max'); load R2;
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    tmp = squeeze(R2(:,:,:,1)); clear R2
                    [e,f,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                    if length(e) ~= 1; e = e(1); f = f(1); end
                    
                    if strcmp(electrode,''); electrode = e;
                    else electrode = eval(cell2mat(electrode)); end
                    if size(electrode) > 1
                        errordlg('invalid electrode choice'); return
                    elseif electrode > size(LIMO.data.chanlocs,2) || electrode < 1
                        errordlg('invalid electrode number'); return
                    end
                    
                    if strcmp(frequency,''); freq_index = f;
                        frequency = freqvect(freq_index);
                    else frequency = eval(cell2mat(frequency)); end
                    if size(frequency) > 1
                        errordlg('invalid frequency choice'); return
                    elseif frequency > LIMO.data.tf_freqs(end) || frequency < LIMO.data.tf_freqs(1)
                        errordlg('invalid frequency number'); return
                    end
                else
                    tmp = squeeze(R2(:,:,1)); clear R2
                    [electrode,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                end
                clear tmp
            else
                electrode = eval(cell2mat(electrode));
                if size(electrode) > 1
                    errordlg('invalid electrode choice'); return
                elseif electrode > size(LIMO.data.chanlocs,2) || electrode < 1
                    errordlg('invalid electrode number'); return
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
                    
                    if mean([LIMO.cache.ERPplot.electrode == electrode ...
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
                    
                    if mean([LIMO.cache.ERPplot.electrode == electrode ...
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
            
            % no cache = compute
            if data_cached == 0
                
                probs = [p/2; 1-p/2];
                z = norminv(probs);
                
                if strcmp(extra,'Original')
                    load Yr;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            if strcmp(LIMO.Analysis,'Time-Frequency')
                                data = squeeze(Yr(electrode,freq_index,:,index{i}));
                                mytitle = sprintf('Original ERSP at \n electrode %s (%g) at %g Hz', LIMO.data.chanlocs(electrode).labels, electrode, frequency);
                            else
                                data = squeeze(Yr(electrode,:,index{i}));
                            end
                            average(i,:) = mean(data,2);
                            se = (std(data') ./ sqrt(numel(index{i})));
                            ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Yr,2));
                        end
                        
                        if strcmp(LIMO.Analysis,'Time')
                            mytitle = sprintf('Original ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                        elseif strcmp(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Original Power Spectrum at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                        end
                    else % continuous variable
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            if strcmp(LIMO.Analysis,'Time-Frequency')
                                continuous(i,:,:) = Yr(electrode,freq_index,:,sorting_values);
                                mytitle{i} = sprintf('Original single trials \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                            else
                                continuous(i,:,:) = Yr(electrode,:,sorting_values);
                                mytitle{i} = sprintf('Original single trials \n sorted by regressor %g \n electrode %s (%g) at %s Hz', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode, frequency);
                            end
                        end
                    end
                    clear Yr
                elseif strcmp(extra,'Modelled')
                    load Betas;
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        Betas = squeeze(Betas(electrode,freq_index,:,:));
                    else
                        Betas = squeeze(Betas(electrode,:,:));
                    end
                    Yh = (LIMO.design.X*Betas')'; % modelled data
                    
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        load Yr;
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            Yr = squeeze(Yr(:,freq_index,:,:));
                        end
                        R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                        
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            data = squeeze(Yh(:,index{i}));
                            average(i,:) = mean(data,2);
                            var   = diag(((R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')) / LIMO.model.model_df(2));
                            CI = sqrt(var/size(index{i},1))*z';
                            ci(i,:,:) = (repmat(mean(data,2),1,2)+CI)';
                        end
                        
                        if strcmp(LIMO.Analysis,'Time')
                            mytitle = sprintf('Modelled ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                        elseif strcmp(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Modelled Power Spectrum at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                        else
                            mytitle = sprintf('Modelled ERSP \n electrode %s (%g) at %g Hz', LIMO.data.chanlocs(electrode).labels, electrode, frequency);
                        end
                    else % continuous variable
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:) = Yh(:,sorting_values);
                            
                            if strcmp(LIMO.Analysis,'Time-Frequency')
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g \n electrode %s (%g) at %g Hz', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode, frequency);
                            else
                                mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                            end
                        end
                    end
                else % Adjusted
                    all = [1:size(LIMO.design.X,2)-1]; all(regressor)=[];
                    if strcmp(LIMO.Analysis,'Time-Frequency')
                        load Yr; Yr = squeeze(Yr(electrode,freq_index,:,:));
                        load Betas; Betas = squeeze(Betas(electrode,freq_index,:,:));
                    else
                        load Yr; Yr = squeeze(Yr(electrode,:,:));
                        load Betas; Betas = squeeze(Betas(electrode,:,:));
                    end
                    confounds = (LIMO.design.X(:,all)*Betas(:,all)')';
                    Ya = Yr - confounds; clear Yr Betas confounds;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            data = squeeze(Ya(:,index{i}));
                            average(i,:) = mean(data,2);
                            se = std(data') ./ sqrt(numel(index{i}));
                            ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Ya,1));
                        end
                        if strcmp(LIMO.Analysis,'Time')
                            mytitle = sprintf('Adjusted ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                        elseif strcmp(LIMO.Analysis,'Frequency')
                            mytitle = sprintf('Adjusted Power Spectrum at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                        else
                            mytitle = sprintf('Adjusted ERSP electrode %s (%g) at %g Hz', LIMO.data.chanlocs(electrode).labels, electrode, frequency);
                        end
                    else % continuous variable
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:) = Ya(:,sorting_values);
                            if strcmp(LIMO.Analysis,'Time-Frequency')
                                mytitle{i} = sprintf('Adjusted single trials \n sorted by regressor \n %g electrode %s (%g) at %g Hz', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode, frequency);
                            else
                                mytitle{i} = sprintf('Adjusted single trials \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
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
                            name = sprintf('Condition_effect_%g.mat',i); load(name);
                            if isfield(LIMO,'cache') && isfield(LIMO,'fig')
                                if strcmp(LIMO.cache.fig.name,name) && ...
                                        LIMO.cache.fig.MCC == MCC && ...
                                        LIMO.cache.fig.threshold == p
                                    if strcmp(LIMO.Analysis,'Time-Frequency')
                                        sig = single(LIMO.cache.fig.mask(electrode,freq_index,:)); sig(find(sig==0)) = NaN;
                                    else
                                        sig = single(LIMO.cache.fig.mask(electrode,:)); sig(find(sig==0)) = NaN;
                                    end
                                end
                            else
                                if strcmp(LIMO.Analysis,'Time-Frequency')
                                    [M, mask, mytitle2] = limo_stat_values_tf(1,name,p,MCC,LIMO,choice);
                                    sig = single(squeeze(mask(electrode,freq_index,:))); sig(find(sig==0)) = NaN;
                                else
                                    [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                                    sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
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
                                name = sprintf('Interaction_effect_%g.mat',i); load(name);
                                if isfield(LIMO,'cache') && isfield(LIMO,'fig')
                                    if strcmp(LIMO.cache.fig.name,name) && ...
                                            LIMO.cache.fig.MCC == MCC && ...
                                            LIMO.cache.fig.threshold == p
                                        if strcmp(LIMO.Analysis,'Time-Frequency')
                                            sig = single(squeeze(LIMO.cache.fig.mask(electrode,freq_index,:))); sig(find(sig==0)) = NaN;
                                        else
                                            sig = single(LIMO.cache.fig.mask(electrode,:)); sig(find(sig==0)) = NaN;
                                        end
                                    end
                                else
                                    if strcmp(LIMO.Analysis,'Time-Frequency')
                                        [M, mask, mytitle2] = limo_stat_values_tf(1,name,p,MCC,LIMO,choice);
                                        sig = single(squeeze(mask(electrode,freq_index,:))); sig(find(sig==0)) = NaN;
                                    else
                                        [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                                        sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
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
                
                LIMO.cache.ERPplot.extra = extra;
                LIMO.cache.ERPplot.average = average;
                LIMO.cache.ERPplot.electrode = electrode;
                LIMO.cache.ERPplot.regressor = regressor;
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.ERPplot.frequency = frequency;
                end
                LIMO.cache.ERPplot.ci = ci;
                LIMO.cache.ERPplot.title = mytitle;
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
                LIMO.cache.ERPplot.electrode = electrode;
                LIMO.cache.ERPplot.regressor = regressor;
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    LIMO.cache.ERPplot.frequency = frequency;
                end
                LIMO.cache.ERPplot.title = mytitle;
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
    
    if ~strncmp(FileName,'LIMO',4) % in all cases but course plot
        
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
                        toplot = []; return
                    else
                        M = LIMO.cache.fig.pval;
                        mytitle = LIMO.cache.fig.title;
                        toplot = LIMO.cache.fig.stats;
                        data_cached = 1;
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                    end
                end
            catch no_cache
                data_cached = 0;
            end
        end
        
        if data_cached == 0 % compute and plot
            
            if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                [M, mask, mytitle] = limo_stat_values_tf(Type,FileName,p,MCC,LIMO,choice,[]);
            else
                [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,choice,[]);
            end
            
            if isempty(mask)
                return
            elseif sum(mask(:)) == 0
                warndlg('  no values under threshold  ','no significant effect','modal');
                toplot = []; return
            else
                assignin('base','p_values',M)
                assignin('base','mask',mask)
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                if strncmp(FileName,'R2',2)
                    toplot = squeeze(R2(:,:,:,1));
                    assignin('base','R_values',toplot)
                elseif strncmp(FileName,'one_sample',10)
                    toplot = squeeze(one_sample(:,:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'two_samples',11)
                    toplot = squeeze(two_samples(:,:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'paired',6)
                    toplot = squeeze(paired_samples(:,:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'Covariate',9);
                    toplot = squeeze(Covariate_effect(:,:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Condition',9);
                    toplot = squeeze(Condition_effect(:,:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'con_',4);
                    toplot = squeeze(con(:,:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'ess_',4);
                    if ~exist('ess','var')
                        effect_nb = eval(FileName(22:end-4));
                        try
                            ess = eval(['ess' num2str(effect_nb)]);
                            clear(['ess' num2str(effect_nb)])
                          catch
                            ess = ess1; clear ess1
                        end
                    end
                    toplot = squeeze(ess(:,:,:,4));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Rep_ANOVA_Interaction',21);
                    toplot = squeeze(Rep_ANOVA_Interaction_with_gp(:,:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Rep_ANOVA_Gp',12);
                    toplot = squeeze(Rep_ANOVA_Gp_effect(:,:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Rep_ANOVA',9);
                    toplot = squeeze(Rep_ANOVA(:,:,:,1));
                    assignin('base','F_values',toplot)
                else
                    disp('file no supported'); return
                end
            else
                if strncmp(FileName,'R2',2)
                    toplot = squeeze(R2(:,:,1));
                    assignin('base','R_values',toplot)
                elseif strncmp(FileName,'one_sample',10)
                    toplot = squeeze(one_sample(:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'two_samples',11)
                    toplot = squeeze(two_samples(:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'paired',6)
                    toplot = squeeze(paired_samples(:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'Covariate',9);
                    toplot = squeeze(Covariate_effect(:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Condition',9);
                    toplot = squeeze(Condition_effect(:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'con_',4);
                    toplot = squeeze(con(:,:,4));
                    assignin('base','T_values',toplot)
                elseif strncmp(FileName,'ess_',4);
                    if ~exist('ess','var')
                        effect_nb = eval(FileName(22:end-4));
                        try
                            ess = eval(['ess' num2str(effect_nb)]);
                            clear(['ess' num2str(effect_nb)])
                        catch
                            ess = ess1; clear ess1
                        end
                    end
                    toplot = squeeze(ess(:,:,4));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Rep_ANOVA_Interaction',21);
                    toplot = squeeze(Rep_ANOVA_Interaction_with_gp(:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Rep_ANOVA_Gp',12);
                    toplot = squeeze(Rep_ANOVA_Gp_effect(:,:,1));
                    assignin('base','F_values',toplot)
                elseif strncmp(FileName,'Rep_ANOVA',9);
                    toplot = squeeze(Rep_ANOVA(:,:,1));
                    assignin('base','F_values',toplot)
                else
                    disp('file no supported'); return
                end
            end
            data_cached = 0;
        end
    end
    % ------------------------------
    %      Image and topoplot
    % ----------------------------
    if Type == 1 || Type == 2
        
        % cache the results for next time
        if data_cached == 0
            LIMO.cache.fig.name       = FileName;
            LIMO.cache.fig.MCC        = MCC;
            LIMO.cache.fig.stats      = toplot;
            LIMO.cache.fig.threshold  = p;
            LIMO.cache.fig.pval       = M;
            LIMO.cache.fig.mask       = mask;
            LIMO.cache.fig.title      = mytitle;
            save LIMO LIMO
        end
        
        
        if Type == 1 && ~strcmp(LIMO.Analysis,'Time-Frequency') && ~strcmp(LIMO.Analysis,'ITC')
            limo_display_image(LIMO,toplot,mask,mytitle)
                        
        elseif Type == 1 && strcmp(LIMO.Analysis,'Time-Frequency') || ...
                Type == 1 && strcmp(LIMO.Analysis,'ITC')
            if ndims(toplot)==3
                limo_display_results_tf(LIMO,toplot,mask,mytitle);
            else
                % plot time*freq map
                mask = squeeze(mask);
                scale = toplot.*mask;
                v = max(scale(:)); [e,f]=find(scale==v);
                if min(scale(:))<0
                    scale(scale==0)=min(scale(:))+(min(scale(:))/10);
                else
                    scale(scale==0)=NaN;
                end
                
                % do the figure
                figure; set(gcf,'Color','w');
                ax(1) = subplot(3,3,[1 2 4 5 7 8]);
                imagesc(scale);
                xlabel('Time','Fontsize',10); ylabel('Frequency','Fontsize',10);
                set(gca,'XTick',LIMO.data.tf_times,'YTick',LIMO.data.tf_freqs)
                title('Time-Frequency Map','FontSize',16);
                color_images_(scale,LIMO);
                set(gca,'layer','top');
                ratio = (LIMO.data.tf_times(end) - LIMO.data.tf_times(1)) / size(toplot,2);
                if LIMO.data.tf_times(1) < 0
                    frame_zeros = round(abs(LIMO.data.tf_times(1) / ratio));
                else
                    frame_zeros = 1;
                end
                
                ax(2) = subplot(3,3,6);
                if length(f)>0
                    f=f(1);
                end
                plot(LIMO.data.tf_freqs,toplot(:,f),'LineWidth',3); grid on; axis tight
                mytitle = sprintf('freq spectrum @ %g ms)', LIMO.data.tf_times(f));
                title(mytitle,'FontSize',12)
                
                ax(3) = subplot(3,3,9);
                if length(e)>0
                    e=e(1);
                end
                plot(LIMO.data.tf_times,toplot(e,:),'LineWidth',3); grid on; axis tight
                mytitle = sprintf('time course @ %g Hz)', LIMO.data.tf_freqs(e));
                title(mytitle,'FontSize',12)
                
                update = 0;
                while update ==0
                    try % use try so that if figure deleted no complain
                        [x,y,button]=ginput(1);
                    catch
                        update = 1; break
                    end
                    if button > 1
                        update = 1;
                    end
                    clickedAx = gca;
                    if clickedAx ~=ax(1)
                        disp('right click to exit')
                    else
                        subplot(3,3,6,'replace');
                        frame = frame_zeros + round(x / ratio);
                        plot(LIMO.data.tf_freqs,toplot(:,frame),'LineWidth',3); grid on; axis tight
                        mytitle = sprintf('freq spectrum @ %g ms)', LIMO.data.tf_times(frame));
                        title(mytitle,'FontSize',12)
                        
                        subplot(3,3,9,'replace'); y = round(y);
                        plot(LIMO.data.tf_times,toplot(y,:),'LineWidth',3); grid on; axis tight
                        mytitle = sprintf('time course @ %g Hz)', LIMO.data.tf_freqs(y));
                        title(mytitle,'FontSize',12)
                        clear x y button
                    end
                end
            end
            
        elseif Type == 2
            %--------------------------
            % topoplot
            %--------------------------
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                errordlg('topoplot not supported for time-frequency analyses')
                
            else
                
                if ~isempty(LIMO.design.electrode)  % not full scalp
                    msgbox('Only one electrode found','No topoplot')
                elseif sum(mask(:)) == 0
                    warndlg('no values under threshold','no significant effect');
                else
                    EEG.data = toplot;
                    EEG.setname = mytitle;
                    EEG.pnts = size(EEG.data,2);
                    EEG.trials = 1;
                    EEG.chanlocs = LIMO.data.chanlocs;
                    EEG.nbchan = size(EEG.data,1);
                    
                    if size(toplot,2) == 1
                        opt = {'maplimits','maxmin','verbose','off'};
                        if isfield(LIMO,'Type')
                            if strcmp(LIMO.Type,'Components')
                                opt = {'maplimits','absmax','electrodes','off','verbose','off'};
                            end
                        end
                        figure;set(gcf,'Color','w','InvertHardCopy','off');
                        topoplot(toplot(:,1),EEG.chanlocs,opt{:});
                        title('Topoplot','FontSize',12)
                    else
                        
                        if strcmp(LIMO.Analysis,'Time')
                            EEG.xmin = LIMO.data.start/1000; % in sec
                            EEG.xmax = LIMO.data.end/1000;   % in sec
                            EEG.times =  (LIMO.data.start/1000:(LIMO.data.sampling_rate/1000):LIMO.data.end/1000); % in sec;
                            pop_topoplot(EEG);
                        elseif strcmp(LIMO.Analysis,'Frequency')
                            EEG.xmin = LIMO.data.freqlist(1);
                            EEG.xmax = LIMO.data.freqlist(end);
                            freqlist = inputdlg('specify frequency range e.g. [5:2:40]','Choose Frequencies to plot');
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
                            
                            N=length(EEG.freq); figure;
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
        end
        
    elseif Type == 3
        
        %--------------------------
        % Course plot
        %--------------------------
        
        if strncmp(FileName,'one_sample',10) || strncmp(FileName,'two_samples',11) || strncmp(FileName,'paired_samples',14) || ...
                strncmp(FileName,'con_',4) || strncmp(FileName,'ess_',4)
            % ------------------------------------------------------------------------------------------------------------
            % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
            % H0 file dim = (electrodes,frames,[t, p],nboot)
            
            if strncmp(FileName,'one_sample',10); data = one_sample; clear one_sample;
            elseif strncmp(FileName,'two_samples',11); data = two_samples; clear two_samples
            elseif  strncmp(FileName,'paired_samples',14);  data = paired_samples; clear paired_samples
            elseif  strncmp(FileName,'con',4); data = con; clear con
            elseif  strncmp(FileName,'ess',4); data = ess; clear ess
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                if size(data,1) > 1
                    if isempty(g.channels)
                        electrode = inputdlg('which electrode to plot','Plotting option');
                    else
                        electrode = g.channels;
                    end
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        tmp = squeeze(data(:,:,:,4));
                        [electrode,freq_index,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                        clear tmp ; frequency = LIMO.data.tf_freqs(freq_index);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > size(data,1)
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
                
                frequency = inputdlg('which frequency to plot','Plotting option');
                if isempty(frequency) || strcmp(cell2mat(frequency),'')
                    if ~exist('freq_index','var')
                        [v,f] = max(squeeze(data(electrode,:,:,4)));
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
                
            else % Time or Freq or ITC
                if size(data,1) > 1
                    if isempty(g.channels)
                        electrode = inputdlg('which electrode to plot','Plotting option');
                    else
                        electrode = g.channels;
                    end
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        [v,e] = max(data(:,:,4)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > size(data,1)
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency')
                trimci(:,2) = squeeze(data(electrode,freq_index,:,1));
                if strncmp(FileName,'ess',4)
                    start_at = max(strfind(FileName,'_'))+1;
                    C = LIMO.contrast{eval(FileName(start_at:end-4))}.C;
                    df = rank(C); % rank of the relevant contrast
                    trimci(:,1) = squeeze(trimci(:,2))-finv(p./2*size(C,1),df,squeeze(data(electrode,frequency,:,3))).*squeeze(data(electrode,frequency,:,2));
                    trimci(:,3) = squeeze(trimci(:,2))+finv(p./2*size(C,1),df,squeeze(data(electrode,frequency,:,3))).*squeeze(data(electrode,frequency,:,2));
                else
                    trimci(:,1) = squeeze(trimci(:,2))-tinv(p./2,squeeze(data(electrode,frequency,:,3))).*squeeze(data(electrode,frequency,:,2));
                    trimci(:,3) = squeeze(trimci(:,2))+tinv(p./2,squeeze(data(electrode,frequency,:,3))).*squeeze(data(electrode,frequency,:,2));
                end
                
            else
                trimci(:,2) = squeeze(data(electrode,:,1));
                if strncmp(FileName,'ess',4)
                    start_at = max(strfind(FileName,'_'))+1;
                    C = LIMO.contrast{eval(FileName(start_at:end-4))}.C;
                    df = rank(C); % rank of the relevant contrast
                    trimci(:,1) = squeeze(trimci(:,2))-finv(p./2*size(C,1),df,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
                    trimci(:,3) = squeeze(trimci(:,2))+finv(p./2*size(C,1),df,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
                else
                    trimci(:,1) = squeeze(trimci(:,2))-(tinv(p./2,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2)))';
                    trimci(:,3) = squeeze(trimci(:,2))+(tinv(p./2,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2)))';
                end
            end
            
            figure;
            set(gcf,'Color','w')
            if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'Time-Frequency')
                if isfield(LIMO.data,'tf_times')
                    timevect = LIMO.data.tf_times;
                    sig = squeeze(single(mask(electrode,freq_index,:))); sig(find(sig==0)) = NaN;
                else
                    timevect = LIMO.data.start:1000/LIMO.data.sampling_rate:LIMO.data.end;
                    sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                end
                plot(timevect,squeeze(trimci(:,2)),'LineWidth',3);
                try
                    fillhandle = patch([timevect fliplr(timevect)], [squeeze(trimci(:,1)),fliplr(squeeze(trimci(:,3)))], [1 0 0]);
                catch no_fill
                    fillhandle = patch([timevect fliplr(timevect)], [squeeze(trimci(:,1));flipud(squeeze(trimci(:,3)))]', [1 0 0]);
                end
                xlabel('Time in ms','FontSize',14)
                ylabel('Amplitude (A.U.)','FontSize',14)
            elseif strcmp(LIMO.Analysis,'Frequency')
                plot(LIMO.data.freqlist,squeeze(trimci(:,2)),'LineWidth',3);
                sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                if size(LIMO.data.freqlist,2) == 1
                    freq_axis = [LIMO.data.freqlist' flipud(LIMO.data.freqlist)'];
                else
                    freq_axis = [LIMO.data.freqlist fliplr(LIMO.data.freqlist)];
                end
                
                try
                    fillhandle = patch(freq_axis, [squeeze(trimci(:,1))' flipud(squeeze(trimci(:,3)))'], [1 0 0]);
                catch no_fill
                    fillhandle = patch(freq_axis, [squeeze(trimci(:,1)) fliplr(squeeze(trimci(:,3)))], [1 0 0]);
                end
                xlabel('Frequency in Hz','FontSize',14)
                ylabel('Spectral Power (A.U.)','FontSize',14)
            end
            
            set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
            grid on; box on; axis tight
            h = axis;  hold on;
            try
                plot(timevect,sig.*h(3),'r.','MarkerSize',20)
            catch no_timevect
                plot(LIMO.data.freqlist,sig.*h(3),'r.','MarkerSize',20)
            end
            
            set(gca,'FontSize',14,'layer','top');
            if size(data,1)>1
                title(sprintf('%s \n electrode %s (%g)',mytitle,LIMO.data.chanlocs(electrode).labels,electrode),'FontSize',16); drawnow;
            else
                title(sprintf('%s \n optimized electrode',mytitle),'FontSize',16); drawnow;
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
            data = load(name); data = getfield(data,cell2mat(fieldnames(data)));
            
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
                        electrode = inputdlg('which electrode to plot','Plotting option');
                    else
                        electrode = g.channels;
                    end
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        tmp = squeeze(data(:,:,:,1));
                        [electrode,freq_index,~] = ind2sub(size(tmp),find(tmp==max(tmp(:))));
                        clear tmp ; frequency = LIMO.data.tf_freqs(freq_index);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > size(data,1)
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
                
                frequency = inputdlg('which frequency to plot','Plotting option');
                if isempty(frequency) || strcmp(cell2mat(frequency),'')
                    if ~exist('freq_index','var')
                        [v,f] = max(squeeze(data(electrode,:,:,1)));
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
                        electrode = inputdlg('which electrode to plot','Plotting option');
                    else
                        electrode = g.channels;
                    end
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        [v,e] = max(data(:,:,1)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > size(data,1)
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
            end
            clear data
            
            % down to business
            % ----------------------
            probs = [p/2; 1-p/2];
            z = norminv(probs);
            load Yr;
            
            if strcmp(extra,'Original')
                if sum(regressor <= categorical) == length(regressor) % for categorical variables
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            data = squeeze(Yr(electrode,freq_index,:,index{i}));
                        else
                            data = squeeze(Yr(electrode,:,index{i}));
                        end
                        average(i,:) = nanmean(data,2);
                        se = (nanstd(data') ./ sqrt(numel(index{i})));
                        ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Yr,2));
                    end
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Original subjects'' parameters at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                    else
                        mytitle = sprintf('Original subjects'' parameters at optimized electrode');
                    end
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [~,sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        reg_values(i,:) = LIMO.data.Cont(sorting_values,i);
                        if strcmp(LIMO.Analysis,'Time-Frequency')
                            continuous(i,:,:) = Yr(electrode,freq_index,:,sorting_values);
                        else
                            continuous(i,:,:) = Yr(electrode,:,sorting_values);
                        end
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Original subjects'' parameters \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                        else
                            mytitle{i} = sprintf('Original subjects'' parameters, \n sorted by regressor %g at optimized electrode', regressor(i));
                        end
                    end
                end
            elseif strcmp(extra,'Modelled')
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    load Betas; Betas = squeeze(Betas(electrode,freq_index,:,:));
                    R = eye(size(Yr,4)) - (LIMO.design.X*pinv(LIMO.design.X));
                else
                    load Betas; Betas = squeeze(Betas(electrode,:,:));
                    R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                end
                Yh = (LIMO.design.X*Betas')'; % modelled data
                
                if sum(regressor <= categorical) == length(regressor) % for categorical variables
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        data = squeeze(Yh(:,index{i}));
                        average(i,:) = nanmean(data,2);
                        index{i} = index{i}(find(~isnan(squeeze(Yr(electrode,1,index{i})))));
                        var   = diag(((R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')) / LIMO.model.model_df(2));
                        CI = sqrt(var/size(index{i},1))*z';
                        ci(i,:,:) = (repmat(nanmean(data,2),1,2)+CI)';
                    end
                    if isempty(LIMO.design.electrode)
                        mytitle = sprintf('Modelled subjects'' parameters at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                    else
                        mytitle = sprintf('Modelled subjects'' parameters at optimized electrode');
                    end
                    clear Yr
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [~,sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        reg_values(i,:) = LIMO.data.Cont(sorting_values);
                        continuous(i,:,:) = Yh(:,sorting_values);
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Modelled single subjects'' parameters \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                        else
                            mytitle{i} = sprintf('Modelled single subjects'' parameters \n sorted by regressor %g at optimized electrode', regressor(i));
                        end
                    end
                end
            else % Adjusted
                all = [1:size(LIMO.design.X,2)-1]; all(regressor)=[];
                if strcmp(LIMO.Analysis,'Time-Frequency')
                    Yr = squeeze(Yr(electrode,freq_index,:,:));
                    load Betas; Betas = squeeze(Betas(electrode,freq_index,:,:));
                    confounds = (LIMO.design.X(:,all)*Betas(:,all)')';
                else
                    Yr = squeeze(Yr(electrode,:,:));
                    load Betas; Betas = squeeze(Betas(electrode,:,:));
                    confounds = (LIMO.design.X(:,all)*Betas(:,all)')';
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
                        mytitle = sprintf('Adjusted subjects'' parameters at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                    else
                        mytitle = sprintf('Adjusted subjects'' parameters at  at optimized electrode');
                    end
                    clear Yr
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [~,sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        reg_values(i,:) = LIMO.data.Cont(sorting_values);
                        continuous(i,:,:) = Ya(:,sorting_values);
                        if isempty(LIMO.design.electrode)
                            mytitle{i} = sprintf('Adjusted single subjects'' parameters \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                        else
                            mytitle{i} = sprintf('Adjusted single subjects'' parameters \n sorted by regressor %g at optimized electrode', regressor(i));
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
                            load(name); [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                            sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                            h = axis;
                            if strcmp(LIMO.Analysis,'Time')
                                plot(timevect,sig.*h(3),'r*','LineWidth',2)
                            else
                                plot(freqvect,sig.*h(3),'r*','LineWidth',2)
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
                                load(name); [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                                sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                                h = axis;
                                if LIMO.analysis_flag == 1
                                    plot(timevect,sig.*h(3),'r*','LineWidth',2)
                                else
                                    plot(freqvect,sig.*h(3),'r*','LineWidth',2)
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
            
            
            
        elseif strncmp(FileName,'Rep_ANOVA',9);   % All stuffs for repeated measures ANOVA
            % -----------------------------------------------------------------------------
            
            
            if strncmp(FileName,'Rep_ANOVA_Factor',16)
                % -----------------------------------
                
                % which ERP to make
                % ------------------
                extra = questdlg('Plotting ERP','ERP Options','Original','Trimmed','Original');
                if isempty(extra)
                    return;
                end
                % -----------------------
                
                [e,f,d]=size(Rep_ANOVA);
                
                if e > 1
                    if isempty(g.channels)
                        electrode = inputdlg('which electrode to plot','Plotting option');
                    else
                        electrode = g.channels;
                    end
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        [v,e] = max(Rep_ANOVA(:,:,1)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > e
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
                
                figure;set(gcf,'Color','w')
                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                effect_nb = eval(FileName(18:end-4));
                C = LIMO.design.C{effect_nb};
                load Yr; Data = squeeze(Yr(electrode,:,:,:));
                
                % compute differences between pairs using C and Cov
                n = size(Data,2);
                if strcmp(extra,'Original')
                    for time_or_freq=1:size(Data,1)
                        avg(time_or_freq,:) = nanmean(C*squeeze(Data(time_or_freq,:,:))',2);
                        S(time_or_freq,:,:) = cov(squeeze(Data(time_or_freq,:,:)));
                    end
                    if e>1
                        mytitle = sprintf('Original %s \n electrode %s (%g)',mytitle,LIMO.data.chanlocs(electrode).labels,electrode);
                    else
                        mytitle = sprintf('Original %s \n optimized electrode',mytitle);
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
                    if e>1
                        mytitle = sprintf('Trimmed %s \n electrode %s (%g)',mytitle,LIMO.data.chanlocs(electrode).labels,electrode);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized electrode',mytitle);
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
                
                if strcmp(LIMO.Analysis,'Time')
                    timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end;
                    plot(timevect,avg,'LineWidth',3);
                    fillhandle = patch([timevect fliplr(timevect)], [c',fliplr(b')], [1 0 0]);
                else
                    freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                    plot(freqvect,avg,'LineWidth',3);
                    fillhandle = patch([freqvect fliplr(freqvect)], [c',fliplr(b')], [1 0 0]);
                end
                
                set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                grid on; box on; axis tight
                sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                h = axis;  hold on;
                if strcmp(LIMO.Analysis,'Time')
                    plot(timevect,sig.*h(3),'r.','MarkerSize',20)
                    xlabel('Time in ms','FontSize',14)
                    ylabel('Amplitude (A.U.)','FontSize',14)
                elseif strcmp(LIMO.Analysis,'Frequency')
                    plot(freqvect,sig.*h(3),'r.','MarkerSize',20)
                    xlabel('Frequency in Hz','FontSize',14)
                    ylabel('Spectral Power (A.U.)','FontSize',14)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                trimci = [c ; avg ; b];
                assignin('base','Plotted_data',trimci);
                
                
                % ----------------------
            elseif strncmp(FileName,'Rep_ANOVA_Gp',12)  %% plot pairs of gp differences
                
                
                % -------------------
                % which ERP to make
                % ------------------
                extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Original');
                if isempty(extra)
                    return;
                end
                % -----------------------
                [e,f,d]=size(Rep_ANOVA_Gp_effect);
                
                % check electrode to plot
                if e > 1
                    if isempty(g.channels)
                        electrode = inputdlg('which electrode to plot','Plotting option');
                    else
                        electrode = g.channels;
                    end
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        clear electrode; [v,e] = max(Rep_ANOVA_Gp_effect(:,:,1)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > e
                            error('electrode number invalid')
                        end
                    end
                else
                    if length(LIMO.design.electrode) == 1
                        electrode = LIMO.design.electrode;
                    else
                        electrode = 1;  % accomodates the fact that all matrices have the electrode dim (even = 1)
                    end
                end
                
                % compute the pair-wise differences and plot
                load Yr;
                if strcmp(extra,'Original')
                    Data = mean(squeeze(Yr(electrode,:,:,:)),3);
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
                Y = nanmean(squeeze(Yr(electrode,:,:,:)),3);
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
                b = Effect - v(electrode,:);
                c = Effect + v(electrode,:);
                if electrode == 1
                    if length(LIMO.design.electrode) == 1
                        if strcmp(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n at electrode %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n at electrode %s (%g)', LIMO.data.chanlocs(LIMO.design.electrode).labels,LIMO.design.electrode);
                        end
                    else
                        if strcmp(extra,'Original')
                            mytitle = sprintf('Mean parameter difference between groups \n optimized electrode');
                        else
                            mytitle = sprintf('Modelled parameter difference between groups \n optimized electrode');
                        end
                    end
                else
                    if strcmp(extra,'Original')
                        mytitle = sprintf('Mean parameter difference between groups \n at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels,electrode);
                    else
                        mytitle = sprintf('Modelled parameter difference between groups \n at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels,electrode);
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
                sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                h = axis;  hold on;
                if LIMO.analysis_flag == 1
                    plot(timevect,sig.*h(3),'r.','MarkerSize',20)
                    xlabel('Time in ms','FontSize',14)
                    ylabel('Amplitude (A.U.)','FontSize',14)
                elseif LIMO.analysis_flag == 2
                    plot(freqvect,sig.*h(3),'r.','MarkerSize',20)
                    xlabel('Frequency in Hz','FontSize',14)
                    ylabel('Spectral Power (A.U.)','FontSize',14)
                end
                set(gca,'FontSize',14,'layer','top');
                title(mytitle,'FontSize',16); drawnow;
                Gp_difference = [c' ; Effect' ; b'];
                assignin('base','Plotted_data',Gp_difference);
                
                
                
                % -------------------------
            elseif strncmp(FileName,'Rep_ANOVA_Interaction',21); % Gp * Repeated measures - plot differences btween condition per gp
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
                        electrode = inputdlg('which electrode to plot','Plotting option');
                    else
                        electrode = g.channels;
                    end
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        [v,e] = max(Rep_ANOVA_Interaction_with_gp(:,:,1)); [v,c]=max(v); electrode = e(c);
                    else
                        electrode = eval(cell2mat(electrode));
                        if length(electrode) > 1
                            error('1 electrode only can be plotted')
                        elseif electrode > e
                            error('electrode number invalid')
                        end
                    end
                else
                    electrode = 1;
                end
                
                % the data to plot are the difference in Yr given LIMO.design.C (see limo_rep_anova)
                effect_nb = eval(FileName(end-4));
                C = LIMO.design.C{effect_nb};
                load Yr; Data = squeeze(Yr(electrode,:,:,:));
                
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
                        mytitle = sprintf('Original %s \n electrode %s (%g)',mytitle,LIMO.data.chanlocs(electrode).labels,electrode);
                    else
                        mytitle = sprintf('Original %s \n optimized electrode',mytitle);
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
                        mytitle = sprintf('Trimmed %s \n electrode %s (%g)',mytitle,LIMO.data.chanlocs(electrode).labels,electrode);
                    else
                        mytitle = sprintf('Trimmed %s \n optimized electrode',mytitle);
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
                sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                h = axis;  hold on;
                if strcmp(LIMO.Analysis,'Time')
                    plot(timevect,sig.*h(3),'r.','MarkerSize',20)
                else
                    plot(freqvect,sig.*h(3),'r.','MarkerSize',20)
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
    
    [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,[],[]);
    
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

if min(scale(:)) >= 0
    cc=cubehelixmap('increase',64);
elseif min(scale(:)) <= 0
    cc=cubehelixmap('decrease',64);   
else
    cc = zeros(64,3);
    tmp = scale.*(scale>0);
    cc(33:64,:)=cubehelixmap('increase',32);
    tmp = scale.*(scale<0);  
    cc(1:32,:)=cubehelixmap('decrease',32);
end

%cc=colormap(jet);
cc(1,:)=[.9 .9 .9]; % set NaNs to gray
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
            ylabel('optimized electrode','FontSize',14);
        end
    end
end
set(gca,'YTickLabel', label_electrodes);
end


%% time vector and label
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [timevect, label] = labels_time_(LIMO,ah)
interval = 50/(1000/LIMO.data.sampling_rate); % in frame
timevect = LIMO.data.start:(1000/LIMO.data.sampling_rate):LIMO.data.end; % in sec
zero_column = find(timevect == 0);
if isempty(zero_column) == 1 % in case it does not encompasses 0
    zero_column = 1;
end

if  LIMO.data.start < 0
    positive_label = timevect(zero_column:interval:end);
    negative_label = timevect(zero_column:-interval:1);
    if negative_label == 0
        negative_label = LIMO.data.start*1000;
    end
    label = [fliplr(negative_label(2:end)) positive_label];
else
    label = timevect(zero_column:interval:end);
end
set(ah,'XTick',find(timevect==label(1)):interval:find(timevect==label(end)),'XTickLabel',round(label));
end


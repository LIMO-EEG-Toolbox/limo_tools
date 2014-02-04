function limo_display_results(Type,FileName,PathName,p,MCC,LIMO,flag)

% This function displays various results
% The arguments specify cases for the
% different kind of figures, thresholds etc ..
%
% FORMAT:
% limo_display_results(Type,FileName,PathName,p,LIMO,flag)
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
%               1=None, 2=2D Cluster, 3=1D Cluster, 4=T max, 5=TFCE
%   LIMO      = LIMO structure
%   flag      = indicates to allow surfing the figure (1) or not (0)
%
% Although the function is mainly intented to be used via the GUI, some figures
% can be generated automatically, for instance limo_display_results(1,'R2.mat',pwd,0.05,5,LIMO,0);
% would load the R2.mar file from the current directory, and plot all
% electrodes/time frames F values threshiolded using tfce at alpha 0.05
% topoplot and ERP like figiures can't be automated since they require user
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
% -----------------------------
%  Copyright (C) LIMO Team 2010


cd(PathName)
load (FileName);
if nargin <= 6
    flag = 1;
end

choice = 'use theoretical p values'; % threshold based on what is computed since H0 is used for clustering
% see limo_stat_values

if LIMO.design.bootstrap == 0
    if MCC == 2 || MCC == 3
        errordlg2('Clustering thresholding necessitates boostrap - invalid choice');
    elseif MCC == 4
        errordlg2('Maximum stat thresholding necessitates bootstrap - invalid choice');
    elseif MCC == 5
        errordlg2('TFCE thresholding necessitates boostrap - invalid choice');
    end
    MCC = 1;
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
                if isfield(LIMO,'cache') && ...
                        strcmp(LIMO.cache.fig.name, FileName) && ...
                        LIMO.cache.fig.MCC == MCC && ...
                        LIMO.cache.fig.threshold == p
                    
                    mask = LIMO.cache.fig.mask;
                    if isempty(mask)
                        return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
                        toplot = []; return
                    else
                        M = LIMO.cache.fig.pval;
                        mytitle = LIMO.cache.fig.title;
                        toplot = LIMO.cache.fig.stats;
                        update_cache = 0;
                        assignin('base','p_values',M)
                        assignin('base','mask',mask)
                    end
                    
                else % compute the plot
                    
                    [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,choice);
                    if isempty(mask)
                        return
                    elseif sum(mask(:)) == 0
                        warndlg('  no values under threshold  ','no significant effect','modal');
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
                
                % do the figure
                % ---------------
                if ~isempty(toplot)
                    figure; set(gcf,'Color','w');
                    
                    % cache the results for next time
                    if update_cache == 1
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
                        % imagesc
                        ax(1) = subplot(3,3,[1 2 4 5 7 8]);
                        if strcmp(LIMO.Analysis,'Time') 
                            timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(toplot,2));
                            ratio =  timevect(2)-timevect(1); % (LIMO.data.end*1000 - LIMO.data.start*1000) / size(toplot,2);
                            if LIMO.data.start < 0
                                frame_zeros = find(timevect == 0);
                                if isempty(frame_zeros)
                                    frame_zeros = round(abs(LIMO.data.start*1000) / ratio)+1;
                                end
                            else
                                frame_zeros = 1;
                            end
                            scale = toplot.*mask; scale(scale==0)=NaN;
                            imagesc(timevect,1:size(toplot,1),scale);
                            
                        elseif strcmp(LIMO.Analysis,'Frequency') 
                            freqvect=linspace((1),LIMO.data.freqlist(end),size(toplot,2));
                            frame_zeros = 1; ratio =  freqvect(2)-freqvect(1);
                            scale = toplot.*mask; scale(scale==0)=NaN;
                            imagesc(freqvect,1:size(toplot,1),scale);
                        end
                        
                        v = max(toplot(:)); [e,f]=find(toplot==v);
                        try
                            caxis([min(min(scale)), max(max(scale))]);
                        catch caxiserror
                        end
                        color_images_(scale,LIMO);
                        title(mytitle,'Fontsize',16)
                        
                        if length(e)>1 % happen if we have multiple times the exact same max values
                            e = e(1); f = f(1); % then we take the 1st (usually an artefact but allows to see it)
                        end
                        
                        % ERP plot at best electrode
                        ax(3) = subplot(3,3,9);
                        
                        if strcmp(LIMO.Analysis,'Time') 
                            mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
                            plot(timevect,toplot(e,:),'LineWidth',3); grid on; axis tight
                        elseif strcmp(LIMO.Analysis,'Frequency') 
                            mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
                            plot(freqvect,toplot(e,:),'LineWidth',3); grid on; axis tight
                        end
                        title(mytitle2,'FontSize',12)
                        
                        % topoplot at max time
                        ax(2) = subplot(3,3,6);
                        chans = LIMO.data.chanlocs;
                        topoplot(toplot(:,f),chans,'maplimits','maxmin');
                        if strcmp(LIMO.Analysis,'Time') 
                            title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
                            set(gca,'XTickLabel', timevect);
                        elseif strcmp(LIMO.Analysis,'Frequency') 
                            title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
                            set(gca,'XTickLabel', LIMO.data.freqlist);
                        end;
                        
                        % update with mouse clicks
                        if flag == 1
                            update = 0;
                            while update ==0
                                try
                                    [x,y,button]=ginput(1);
                                catch
                                    update =1; break
                                end
                                if button > 1
                                    update = 1;
                                end
                                clickedAx = gca;
                                if clickedAx ~=ax(1)
                                    disp('right click to exit')
                                else
                                    frame = frame_zeros + round(x / ratio);
                                    % ERP plot at best electrode and topoplot
                                    % at max time or freq
                                    y = round(y);
                                    if strcmp(LIMO.Analysis,'Time') ;  
                                        subplot(3,3,6,'replace');
                                        topoplot(toplot(:,frame),LIMO.data.chanlocs);
                                        title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',12)

                                        subplot(3,3,9,'replace');
                                        plot(timevect,toplot(y,:),'LineWidth',3); grid on; axis tight
                                        mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
                                        title(mytitle2,'FontSize',12); 
                                        
                                    elseif strcmp(LIMO.Analysis,'Frequency') 
                                        subplot(3,3,6,'replace');
                                        topoplot(toplot(:,frame),LIMO.data.chanlocs);
                                        title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)

                                        subplot(3,3,9,'replace');
                                        plot(freqvect,toplot(y,:),'LineWidth',3); grid on; axis tight
                                        mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,LIMO.data.chanlocs(y).urchan);
                                        title(mytitle2,'FontSize',12); 
                                        
                                    end
                                    try
                                        mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.expected_chanlocs(y).labels,LIMO.data.expected_chanlocs(y).urchan);
                                    catch ME
                                        mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,LIMO.data.chanlocs(y).urchan);
                                    end
                                    title(mytitle2,'FontSize',12)
                                end
                            end
                        end
                        
                    else % strcmp(LIMO.Analysis,'Time-Frequency')  - 3D maps
                        limo_time_freq_display(toplot,timevect,freqvect,mytitle);
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
                
                % plots
                figure; set(gcf,'Color','w');
                % imagesc eigen values
                h = subplot(3,3,[4 5 7 8]);
                timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(EV,2));
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
                v = max(toplot(:)); [e,f]=find(toplot==v);
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
            
            EEG.pnts = size(EEG.data,2);
            EEG.nbchan = size(EEG.data,1);
            EEG.trials = 1;
            EEG.chanlocs = LIMO.data.chanlocs;
            if LIMO.analysis_flag == 1
                EEG.xmin = LIMO.data.start;
                EEG.xmax = LIMO.data.end;
                EEG.times = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
            else
                EEG.xmin = LIMO.data.freqlist(1);
                EEG.xmax = LIMO.data.freqlist(end);
                EEG.times = LIMO.data.freqlist;
            end
            
            
            if strcmp(FileName,'R2.mat')
                EEG.data = squeeze(R2(:,:,2));
                EEG.setname = 'R2 - F values';
                pop_topoplot(EEG); assignin('base',EEG.setname,EEG.data);
            elseif strncmp(FileName,'Condition_effect',16)
                EEG.data = squeeze(Condition_effect(:,:,1));
                EEG.setname = sprintf('Condition %s - F values',FileName(18:end-4));
                pop_topoplot(EEG); assignin('base',EEG.setname,EEG.data);
            elseif strncmp(FileName,'Covariate_effect',16)
                EEG.data = squeeze(Covariate_effect(:,:,1));
                EEG.setname = sprintf('Covariate %s - F values',FileName(18:end-4));
                pop_topoplot(EEG); assignin('base',EEG.setname,EEG.data);
            elseif strncmp(FileName,'Interaction_effect',18)
                EEG.data = squeeze(Interaction_effect(:,:,1));
                EEG.setname = sprintf('Interaction %s - F values',FileName(20:end-4));
                pop_topoplot(EEG); assignin('base',EEG.setname,EEG.data);
            elseif strcmp(FileName,'Partial_coef.mat')
                regressor = str2num(cell2mat(inputdlg('which regressor to plot (e.g. 1:3)','Plotting option')));
                if max(regressor) > size(Partial_coef,3); errordlg('error in regressor number'); return; end
                for b = regressor
                    EEG.data = squeeze(Partial_coef(:,:,b,1));
                    EEG.setname = sprintf('Partial coef R2 values variable %g',b);
                    EEG.pnts = size(EEG.data,2);
                    EEG.xmin = LIMO.data.start;
                    EEG.xmax = LIMO.data.end;
                    EEG.times = timevect;
                    EEG.trials = 1;
                    EEG.chanlocs = LIMO.data.chanlocs;
                    EEG.nbchan = size(EEG.data,1);
                    pop_topoplot(EEG);
                    tmp_name = sprintf('Plotted_data_partial_coef_%g',b);
                    assignin('base',[tmp_name],squeeze(Partial_coef(:,:,b,1)))
                end
            elseif strcmp(FileName(1:4),'con_')
                toplot = squeeze(con(:,:,4));
                EEG.setname = ['Contrast ',[FileName(5:end-4)],' -- T values'];
                pop_topoplot(EEG); assignin('base',EEG.setname,EEG.data);
            elseif strcmp(FileName(1:4),'ess_')
                toplot = squeeze(ess(:,:,end-1));
                EEG.setname = ['Contrast ',[FileName(5:end-4)],' -- F values'];
                pop_topoplot(EEG); assignin('base',EEG.setname,EEG.data);
            else
                disp('file not supported');
            end
            
            
        case{3}
            
            %--------------------------
            % Time course / Power
            %--------------------------
            
            % which variable(s) to plot
            % ----------------------
            input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2));
            regressor = inputdlg(input_title,'Plotting option'); if isempty(regressor); return; end
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
            extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Adjusted','Adjusted');
            if isempty(extra)
                return;
            elseif strcmp(extra,'Original')
                if regressor == size(LIMO.design.X,2)
                    errordlg('you can''t plot adjusted mean for original data'); return
                end
            end
            
            % which electrode to plot
            % ------------------
            electrode = inputdlg('which electrode to plot','Plotting option');
            try
                electrode = eval(cell2mat(electrode));
                if size(electrode) > 1
                    errordlg('invalid electrode number'); return
                elseif electrode > size(LIMO.data.chanlocs,2)
                    errordlg('invalid electrode number'); return
                end
            catch ME
                load R2; [v,e] = max(R2(:,:,1)); [v,c]=max(v); electrode = e(c);
                % [v,electrode] = max(max(mean(Yr,3),[],2)); % plot at the max of Yr
            end
            
            % timing info
            if LIMO.analysis_flag  == 1
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
            elseif LIMO.analysis_flag  == 2
                freqvect=LIMO.data.freqlist;
            elseif strcmp(LIMO.Analysis,'Time-Frequency') 
                timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                freqvect = LIMO.data.freqlist;
            end
            
            % down to business
            % ----------------------
            probs = [p/2; 1-p/2];
            z = norminv(probs);
            
            if strcmp(extra,'Original')
                load Yr;
                if sum(regressor <= categorical) == length(regressor) % for categorical variables
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        data = squeeze(Yr(electrode,:,index{i}));
                        average(i,:) = mean(data,2);
                        se = (std(data') ./ sqrt(numel(index{i})));
                        ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Yr,2));
                    end
                    mytitle = sprintf('Original ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        continuous(i,:,:) = Yr(electrode,:,sorting_values);
                        mytitle{i} = sprintf('Original single trials \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                    end
                end
                clear Yr
            elseif strcmp(extra,'Modelled')
                load Betas; Betas = squeeze(Betas(electrode,:,:));
                Yh = (LIMO.design.X*Betas')'; % modelled data
                load Yr; R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                
                if sum(regressor <= categorical) == length(regressor) % for categorical variables
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        data = squeeze(Yh(:,index{i}));
                        average(i,:) = mean(data,2);
                        var   = diag(((R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')) / LIMO.model.model_df(2));
                        CI = sqrt(var/size(index{i},1))*z';
                        ci(i,:,:) = (repmat(mean(data,2),1,2)+CI)';
                    end
                    mytitle = sprintf('Modelled ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        continuous(i,:,:) = Yh(:,sorting_values);
                        mytitle{i} = sprintf('Modelled single trials \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                    end
                end
            else % Adjusted
                all = [1:size(LIMO.design.X,2)-1]; all(regressor)=[];
                load Yr; Yr = squeeze(Yr(electrode,:,:));
                load Betas; Betas = squeeze(Betas(electrode,:,:));
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
                    mytitle = sprintf('Adjusted ERP at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                else % continuous variable
                    for i=1:length(regressor)
                        index{i} = find(LIMO.design.X(:,regressor(i)));
                        [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                        continuous(i,:,:) = Ya(:,sorting_values);
                        mytitle{i} = sprintf('Adjusted single trials \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                    end
                end
            end
            
            % make the figure(s)
            % ------------------
            figure;set(gcf,'Color','w')
            if sum(regressor <= categorical) == length(regressor)
                for i=1:size(average,1)
                    if LIMO.analysis_flag  == 1
                        plot(timevect,average(i,:),'LineWidth',1.5); hold on
                    elseif LIMO.analysis_flag  == 2
                        plot(freqvect,average(i,:),'LineWidth',1.5); hold on
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
                            name = sprintf('Condition_effect_%g.mat',i); load(name);
                            if isfield(LIMO,'cache') && ...
                                    strcmp(LIMO.cache.fig.name,name) && ...
                                    LIMO.cache.fig.MCC == MCC && ...
                                    LIMO.cache.fig.threshold == p
                                sig = single(LIMO.cache.fig.mask(electrode,:)); sig(find(sig==0)) = NaN;
                            else
                                [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                                sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                            end
                            h = axis;
                            if LIMO.analysis_flag  == 1
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
                                name = sprintf('Interaction_effect_%g.mat',i); load(name);
                                if isfield(LIMO,'cache') && ...
                                        strcmp(LIMO.cache.fig.name,name) && ...
                                        LIMO.cache.fig.MCC == MCC && ...
                                        LIMO.cache.fig.threshold == p
                                    sig = single(LIMO.cache.fig.mask(electrode,:)); sig(find(sig==0)) = NaN;
                                else
                                    [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                                    sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                                end
                                h = axis;
                                if LIMO.analysis_flag  == 1
                                    plot(timevect,sig.*h(3),'r*','LineWidth',2)
                                elseif LIMO.analysis_flag  == 2
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
                title(mytitle,'FontSize',19); drawnow;
                assignin('base','Plotted_data', average)
                v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
                set(gca,'FontSize',14);
                if LIMO.analysis_flag  == 1
                    xlabel('Time in ms','FontSize',16)
                    ylabel('Amplitude in {\mu}V','FontSize',16)
                elseif LIMO.analysis_flag  == 2
                    xlabel('Freq in Hz','FontSize',16)
                    ylabel('Power spectral density','FontSize',16);
                end
                
                LIMO.cache.ERPplot = average;
                LIMO.cache.ERPplot.ci = ci;
                LIMO.cache.ERPplot.title = mytitle;
                
            else
                for i=1:size(continuous,1)
                    if i > 1; figure;set(gcf,'Color','w'); end
                    index = find(~isnan(squeeze(continuous(i,1,:))));
                    
                    if strcmp(LIMO.Analysis,'Time') 
                        surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Time in ms','FontSize',16)
                        zlabel('Amplitude in {\mu}V','FontSize',16)
                    elseif strcmp(LIMO.Analysis,'Frequency') 
                        surf(index,freqvect,squeeze(continuous(i,:,index)));shading interp
                        ylabel('Frequency in Hz','FontSize',16)
                        zlabel('Spectral power in {\mu}V^2/Hz','FontSize',16)
                    end
                    % --
                    axis tight; title(mytitle{i},'FontSize',14); drawnow;
                    xlabel('Sorted trials','FontSize',16)
                    try
                        set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                    end
                end
                
                LIMO.cache.ERPplot = squeeze(continuous(i,:,index));
                LIMO.cache.ERPplot.title = mytitle{i};
            end
            
    end % closes switch
    
    
    
    
    
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % -------------------                               LEVEL 2
    % -------------------                            GROUP EFFECTs
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
elseif LIMO.Level == 2
    
    if ~strncmp(FileName,'LIMO',4) % in all cases but ERP plot for GLM
        
        % if previously plotted, recover from the cache
        if isfield(LIMO,'cache') && ...
                strcmp(LIMO.cache.fig.name, FileName) && ...
                LIMO.cache.fig.MCC == MCC && ...
                LIMO.cache.fig.threshold == p
            
            mask = LIMO.cache.fig.mask;
            if isempty(mask)
                return
            elseif sum(mask(:)) == 0
                warndlg('  no values under threshold  ','no significant effect','modal');
                toplot = []; return
            else
                M = LIMO.cache.fig.pval;
                mytitle = LIMO.cache.fig.title;
                toplot = LIMO.cache.fig.stats;
                update_cache = 0;
                assignin('base','p_values',M)
                assignin('base','mask',mask)
            end
            
        else % compute and plot
            
            [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,choice,[]);
            if isempty(mask)
                return
            elseif sum(mask(:)) == 0
                warndlg('  no values under threshold  ','no significant effect','modal');
                toplot = []; return
            else
                assignin('base','p_values',M)
                assignin('base','mask',mask)
            end
            
            if LIMO.analysis_flag == 3
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
            update_cache = 1;
        end
        
        % ------------------------------
        %      Image and topoplot
        % ----------------------------
        if Type == 1 || Type == 2
            
            if Type == 1 && LIMO.analysis_flag ~= 3
                %--------------------------
                % imagesc of the results
                %--------------------------
                % what to plot
                scale = toplot.*mask;
                v = max(scale(:)); [e,f]=find(scale==v);
                if min(scale(:))<0
                    scale(scale==0)=min(scale(:))+(min(scale(:))/10);
                else
                    scale(scale==0)=NaN;
                end
                
                % do the figure
                figure; set(gcf,'Color','w');
                
                % cache the results for next time
                if update_cache == 1
                    LIMO.cache.fig.name       = FileName;
                    LIMO.cache.fig.MCC        = MCC;
                    LIMO.cache.fig.stats      = toplot;
                    LIMO.cache.fig.threshold  = p;
                    LIMO.cache.fig.pval       = M;
                    LIMO.cache.fig.mask       = mask;
                    LIMO.cache.fig.title      = mytitle;
                    save LIMO LIMO
                end
                
                ax(1) = subplot(3,3,[1 2 4 5 7 8]);
                if strcmp(LIMO.Analysis,'Time') 
                    timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(toplot,2));
                    ratio = (LIMO.data.end*1000 - LIMO.data.start*1000) / size(toplot,2);
                    if LIMO.data.start < 0
                        frame_zeros = round(abs(LIMO.data.start*1000) / ratio);
                    end
                    imagesc(timevect,1:size(toplot,1),scale);
                    
                else
                    freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                    imagesc(freqvect,1:size(M,1),scale);
                end
                title(mytitle,'FontSize',16);
                color_images_(scale,LIMO);
                set(gca,'layer','top');
                
                if size(toplot,1)>1
                    ax(2) = subplot(3,3,6);
                    if length(f)>0
                        f=f(1);
                    end
                    if LIMO.analysis_flag == 1
                        topoplot(toplot(:,f),LIMO.data.chanlocs);
                        title(['topoplot @' num2str(timevect(f)) 'ms'],'FontSize',12)
                    else
                        topoplot(toplot(:,f),LIMO.data.chanlocs);
                        title(['topoplot @' num2str(freqvect(f)) 'Hz'],'FontSize',12)
                    end
                end
                
                ax(3) = subplot(3,3,9);
                if length(e)>0
                    e=e(1);
                end
                
                if LIMO.analysis_flag == 1
                    plot(timevect,toplot(e,:),'LineWidth',3); grid on; axis tight
                    if size(toplot,1)>1
                        mytitle = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
                    else
                        mytitle = sprintf('time course of \n optimized electrode');
                    end
                elseif LIMO.analysis_flag == 2
                    plot(freqvect,toplot(e,:),'LineWidth',3); grid on; axis tight
                    if size(toplot,1)>1
                        mytitle = sprintf('frequencies @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
                    else
                        mytitle = sprintf('frequencies of \n optimized electrode');
                    end
                end
                title(mytitle,'FontSize',12)
                
                if size(toplot,1)>1
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
                            % topoplot
                            subplot(3,3,6,'replace');
                            frame = frame_zeros + round(x / ratio);
                            topoplot(toplot(:,frame),LIMO.data.chanlocs);
                            if LIMO.analysis_flag == 1
                                title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',12)
                            else
                                title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
                            end
                            % ERP/Power
                            subplot(3,3,9,'replace'); y = round(y);
                            if LIMO.analysis_flag == 1
                                plot(timevect,toplot(y,:),'LineWidth',3); grid on, axis tight
                                if size(toplot,1)>1
                                    mytitle = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
                                else
                                    mytitle = sprintf('time course of \n optimized electrode');
                                end
                            else
                                plot(freqvect,toplot(y,:),'LineWidth',3); grid on, axis tight
                                if size(toplot,1)>1
                                    mytitle = sprintf('frequencies @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
                                else
                                    mytitle = sprintf('frequencies of \n optimized electrode');
                                end
                            end
                            title(mytitle,'FontSize',12)
                            clear x y button
                        end
                    end
                end
                
            elseif Type == 2 && LIMO.analysis_flag ~= 3
                %--------------------------
                % topoplot
                %--------------------------
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
                    if LIMO.analysis_flag == 1
                        EEG.xmin = LIMO.data.start;
                        EEG.xmax = LIMO.data.end;
                        EEG.times =  (LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000); % in sec;
                    else
                        EEG.xmin = LIMO.data.freqlist(1);
                        EEG.xmax = LIMO.data.freqlist(end);
                        EEG.times = LIMO.data.freqlist;
                    end
                    pop_topoplot(EEG);
                    assignin('base','Plotted_data',EEG.data)
                end
            end
            
            
        elseif Type == 3 && LIMO.analysis_flag ~= 3
            
            %--------------------------
            % ERP
            %--------------------------
            
            if strncmp(FileName,'one_sample',10) || strncmp(FileName,'two_samples',11) || strncmp(FileName,'paired_samples',14) || ...
                    strncmp(FileName,'con_',4) || strncmp(FileName,'ess_',4);
                % ------------------------------------------------------------------------------------------------------------
                % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
                % H0 file dim = (electrodes,frames,[t, p],nboot)
                
                if strncmp(FileName,'one_sample',10); [e,f,d]=size(one_sample); data = one_sample;
                elseif strncmp(FileName,'two_samples',11); [e,f,d]=size(two_samples); data = two_samples;
                elseif  strncmp(FileName,'paired_samples',14);  [e,f,d]=size(paired_samples); data = paired_samples;
                elseif  strncmp(FileName,'con',4);  [e,f,d]=size(con); data = con;
                elseif  strncmp(FileName,'ess',4);  [e,f,d]=size(ess); data = ess;
                end
                
                if e > 1
                    electrode = inputdlg('which electrode to plot','Plotting option');
                    if isempty(electrode) || strcmp(cell2mat(electrode),'')
                        [v,e] = max(data(:,:,4)); [v,c]=max(v); electrode = e(c);
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
                
                trimci(:,:,2) = squeeze(data(electrode,:,1));
                if strncmp(FileName,'ess',4)
                    start_at = max(strfind(FileName,'_'))+1;
                    C = LIMO.contrast{eval(FileName(start_at:end-4))}.C;
                    df = rank(C); % rank of the relevant contrast
                    trimci(:,:,1) = squeeze(trimci(:,:,2))-finv(p./2*size(C,1),df,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
                    trimci(:,:,3) = squeeze(trimci(:,:,2))+finv(p./2*size(C,1),df,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
                else
                    trimci(:,:,1) = squeeze(trimci(:,:,2))-tinv(p./2,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
                    trimci(:,:,3) = squeeze(trimci(:,:,2))+tinv(p./2,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
                end
                
                figure;
                set(gcf,'Color','w')
                if LIMO.analysis_flag == 1
                    timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                    plot(timevect,squeeze(trimci(:,:,2)),'LineWidth',3);
                    xlabel('Time in ms','FontSize',14)
                    ylabel('Amplitude (A.U.)','FontSize',14)
                else
                    freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                    plot(timevect,squeeze(trimci(:,:,2)),'LineWidth',3);
                    xlabel('Frequency in Hz','FontSize',14)
                    ylabel('Spectral Power (A.U.)','FontSize',14)
                end
                
                fillhandle = patch([timevect fliplr(timevect)], [squeeze(trimci(:,:,1)),fliplr(squeeze(trimci(:,:,3)))], [1 0 0]);
                set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
                grid on; box on; axis tight
                sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                h = axis;  hold on; plot(timevect,sig.*h(3),'r.','MarkerSize',20)
                set(gca,'FontSize',14,'layer','top');
                if e>1
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
                    regressor = inputdlg(input_title,'Plotting option'); if isempty(regressor); return; end
                    try regressor = sort(eval(cell2mat(regressor)));
                        if max(regressor) > size(LIMO.design.X,2); errordlg('invalid regressor number'); end
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
                
                % which electrode to plot
                % ------------------
                if isempty(LIMO.design.electrode)
                    electrode = inputdlg('which electrode to plot','Plotting option');
                    try
                        electrode = eval(cell2mat(electrode));
                        if size(electrode) > 1
                            errordlg('invalid electrode number'); return
                        elseif electrode > size(LIMO.data.chanlocs,2)
                            errordlg('invalid electrode number'); return
                        end
                    catch ME
                        load R2; [v,e] = max(R2(:,:,1)); [v,c]=max(v); electrode = e(c);
                        % [v,electrode] = max(max(mean(Yr,3),[],2)); % plot at the max of Yr
                    end
                else
                    electrode =1;  % for optimized electrode analyses
                end
                
                
                % down to business
                % ----------------------
                probs = [p/2; 1-p/2];
                z = norminv(probs);
                
                if strcmp(extra,'Original')
                    load Yr;
                    if sum(regressor <= categorical) == length(regressor) % for categorical variables
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            data = squeeze(Yr(electrode,:,index{i}));
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
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                            continuous(i,:,:) = Yr(electrode,:,sorting_values);
                            if isempty(LIMO.design.electrode)
                                mytitle{i} = sprintf('Original subjects'' parameters \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                            else
                                mytitle{i} = sprintf('Original subjects'' parameters, \n sorted by regressor %g at optimized electrode', regressor(i));
                            end
                        end
                    end
                    clear Yr
                elseif strcmp(extra,'Modelled')
                    load Betas; Betas = squeeze(Betas(electrode,:,:));
                    Yh = (LIMO.design.X*Betas')'; % modelled data
                    load Yr; R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
                    
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
                    else % continuous variable
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
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
                    load Yr; Yr = squeeze(Yr(electrode,:,:));
                    load Betas; Betas = squeeze(Betas(electrode,:,:));
                    confounds = (LIMO.design.X(:,all)*Betas(:,all)')';
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
                    else % continuous variable
                        for i=1:length(regressor)
                            index{i} = find(LIMO.design.X(:,regressor(i)));
                            [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
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
                        if LIMO.analysis_flag == 1
                            timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
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
                                if LIMO.analysis_flag == 1
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
                        if LIMO.analysis_flag == 1
                            timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                            surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
                            ylabel('Time in ms','FontSize',16)
                            zlabel('Amplitude (A.U.)','FontSize',16)
                        else
                            freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                            surf(index,freqvect,squeeze(continuous(i,:,index)));shading interp
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
                        electrode = inputdlg('which electrode to plot','Plotting option');
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
                            avg(time,:) = nanmean(C*squeeze(Data(time_or_freq,:,:))',2);
                            S(time,:,:) = cov(squeeze(Data(time_or_freq,:,:)));
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
                            S(time,:,:) = cov(SD); % winsorized covariance
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
                    % c = avg + 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time,:,:))*C'))); % uses Bonferoni inequality
                    % b = avg - 2*(finv(p./(2*size(C,1)),df,dfe).*(sqrt(C*squeeze(S(time,:,:))*C')));
                    c = avg + tinv(p./(2*size(C,1)),dfe).*(sqrt(C*squeeze(S(time,:,:))*C'));
                    b = avg - tinv(p./(2*size(C,1)),dfe).*(sqrt(C*squeeze(S(time,:,:))*C'));
                    
                    if LIMO.analysis_flag == 1
                        timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in ms
                        plot(timevect,avg,'LineWidth',3);
                        fillhandle = patch([timevect fliplr(timevect)], [c',fliplr(b')], [1 0 0]);
                    elseif LIMO.analysis_flag == 2
                        freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
                        plot(freqvect,avg,'LineWidth',3);
                        fillhandle = patch([freqvect fliplr(freqvect)], [c',fliplr(b')], [1 0 0]);
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
                        electrode = inputdlg('which electrode to plot','Plotting option');
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
                    if LIMO.analysis_flag == 1
                        timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
                        plot(timevect,Effect,'LineWidth',3);
                        fillhandle = patch([timevect fliplr(timevect)], [c',fliplr(b')], [1 0 0]);
                    elseif LIMO.analysis_flag == 2
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
                        electrode = inputdlg('which electrode to plot','Plotting option');
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
                    for gp = 1:LIMO.design.nb_conditions
                        index = find(LIMO.data.Cat==gp);
                        dfe = size(Data,2)-length(index)+1;
                        if gp==1
                            if LIMO.analysis_flag == 1
                                plot(timevect,avg(gp,:),'LineWidth',3);
                            elseif LIMO.analysis_flag == 2
                                plot(freqvect,avg(gp,:),'LineWidth',3);
                            end
                            colorOrder = get(gca, 'ColorOrder');
                        else
                            if LIMO.analysis_flag == 1
                                plot(timevect,avg(gp,:),'Color',colorOrder(gp,:),'LineWidth',3);
                            else
                                plot(freqvect,avg(gp,:),'Color',colorOrder(gp,:),'LineWidth',3);
                            end
                        end
                        c = avg(gp,:,:) + tinv(p./(2*size(C,1)),dfe).*(sqrt(C*squeeze(S(gp,:,:,:))*C'));
                        b = avg(gp,:,:) - tinv(p./(2*size(C,1)),dfe).*(sqrt(C*squeeze(S(gp,:,:,:))*C'));
                        if LIMO.analysis_flag == 1
                            fillhandle = patch([timevect fliplr(timevect)], [c,fliplr(b)], colorOrder(gp,:));
                            xlabel('Time in ms','FontSize',14)
                            ylabel('Amplitude (A.U.)','FontSize',14)
                        elseif LIMO.analysis_flag == 2
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
                    if LIMO.analysis_flag == 1
                        plot(timevect,sig.*h(3),'r.','MarkerSize',20)
                    elseif LIMO.analysis_flag == 2
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
    end
    
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
            EEG.xmin = LIMO.data.start;
            EEG.xmax = LIMO.data.end;
            EEG.times =  (LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000); % in sec;
            EEG.trials = 1;
            EEG.chanlocs = LIMO.data.chanlocs;
            EEG.nbchan = size(EEG.data,1);
            pop_topoplot(EEG);
            assignin('base','Plotted_data',EEG.data)
        end
        
    elseif Type == 3
        disp('no ERP Plots to Lateraliztion'); % we could but nobody asked for it
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

cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
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
    for i = 1 : length(LIMO.data.chanlocs)
        try
            label_electrodes{i} = LIMO.data.expected_chanlocs(i).labels;
        catch ME
            label_electrodes{i} = LIMO.data.chanlocs(i).labels;
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
timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
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


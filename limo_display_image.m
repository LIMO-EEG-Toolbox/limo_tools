function limo_display_image(FileName,PathName,p,MCC,LIMO,flag)

% This function displays images with a intensity plotted as function of
% time or frequency (x) and electrodes (y) - for ERSP it precomputes what
% needs to be plotted and call limo_display_image_tf
%
% FORMAT: limo_display_results(FileName,PathName,p,MCC,LIMO,flag)
%
% INPUTS:
%   Filename  = Name of the file to image
%   PathName  = Path of the file to image
%   p         = threshold p value e.g. 0.05
%   MCC       = Multiple Comparison technique
%               1=None, 2= Cluster, 3=TFCE, 4=T max
%   LIMO      = LIMO structure
%   flag      = indicates to allow surfing the figure (1) or not (0)
%
% Although the function is mainly intented to be used via the GUI, figures
% can be generated automatically:
% imo_display_results(1,'R2.mat',pwd,0.05,5,LIMO,0);
% will load the R2.mat file from the current directory, and plot all
% electrodes/time frames F values thresholded using tfce at alpha 0.05
%
% the function comes from limo_display_results (deprecated)
% Cyril Pernet v1 August 2015
% see also limo_stat_values topoplot
% ----------------------------------
%  Copyright (C) LIMO Team 2015

cd(PathName)
load (FileName);
if nargin <= 5
    flag = 1;
end

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
    MCC =2;
end

% -------------------------------------------------------------------------
% -------------------      LEVEL 1     ------------------------------------
% -------------------  SINGLE SUBJECT  ------------------------------------
% -------------------------------------------------------------------------
if LIMO.Level == 1
    
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
        
        % do the figure
        % ---------------
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
                figure; set(gcf,'Color','w','InvertHardCopy','off');
                
                % imagesc
                ax(1) = subplot(3,3,[1 2 4 5 7 8]);
                if strcmp(LIMO.Analysis,'Time')
                    try timevect = LIMO.data.timevect; catch timevect = []; end
                    if size(timevect,2) == 1; timevect = timevect'; end
                    if size(timevect,2) ~= size(toplot,2);
                        timevect = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
                    end
                    ratio =  (timevect(end)-timevect(1)) / length(timevect); % this the diff in 'size' between consecutive frames
                    if LIMO.data.start < 0
                        frame_zeros = find(timevect == 0);
                        if isempty(frame_zeros)
                            frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
                        end
                    else
                        frame_zeros = 1;
                    end
                    scale = toplot.*mask; scale(scale==0)=NaN;
                    imagesc(timevect,1:size(toplot,1),scale);
                    
                elseif strcmp(LIMO.Analysis,'Frequency')
                    freqvect = LIMO.data.freqlist;
                    if size(freqvect,2) == 1; freqvect = freqvect'; end
                    if size(freqvect,2) ~= size(toplot,2)
                        freqvect = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
                    end
                    frame_zeros = 1;
                    ratio =  (freqvect(end)-freqvect(1)) / length(freqvect);
                    scale = toplot.*mask; scale(scale==0)=NaN;
                    imagesc(freqvect,1:size(toplot,1),scale);
                end
                
                v = max(toplot(:)); [e,f]=find(toplot==v);
                try
                    caxis([min(min(scale)), max(max(scale))]);
                catch caxiserror
                end
                color_images_(scale,LIMO);
                title(mytitle,'Fontsize',14)
                
                if length(e)>1 % happen if we have multiple times the exact same max values
                    e = e(1); f = f(1); % then we take the 1st (usually an artefact but allows to see it)
                end
                
                % ERP plot at best electrode
                ax(3) = subplot(3,3,9);
                
                if strcmp(LIMO.Analysis,'Time')
                    if strcmp(LIMO.Type,'Components')
                        mytitle2 = sprintf('time course @ \n component %g', e);
                    elseif strcmp(LIMO.Type,'Channels')
                        mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
                    end
                    plot(timevect,toplot(e,:),'LineWidth',3); grid on; axis tight
                elseif strcmp(LIMO.Analysis,'Frequency')
                    if strcmp(LIMO.Type,'Components')
                        mytitle2 = sprintf('power spectra @ \n component %g', e);
                    elseif strcmp(LIMO.Type,'Channels')
                        mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
                    end
                    plot(freqvect,toplot(e,:),'LineWidth',3); grid on; axis tight
                end
                title(mytitle2,'FontSize',12)
                
                % topoplot at max time
                ax(2) = subplot(3,3,6);
                if strcmp(LIMO.Type,'Components')
                    EEG=pop_loadset([LIMO.data.data_dir LIMO.data.data]);
                    opt = {'maplimits','absmax','electrodes','off','verbose','off'};
                    topoplot(toplot(:,f),EEG.chanlocs,opt{:});
                    
                else
                    chans = LIMO.data.chanlocs;
                    topoplot(toplot(:,f),chans,'maplimits','maxmin','verbose','off');
                    if strcmp(LIMO.Analysis,'Time')
                        title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
                        set(gca,'XTickLabel', timevect);
                    elseif strcmp(LIMO.Analysis,'Frequency')
                        title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
                        set(gca,'XTickLabel', LIMO.data.freqlist);
                    end
                end
                
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
                                if strcmp(LIMO.Type,'Components')
                                    mytitle2 = sprintf('time course @ \n component %g', y);
                                elseif strcmp(LIMO.Type,'Channels')
                                    mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                                end
                                title(mytitle2,'FontSize',12);
                                
                            elseif strcmp(LIMO.Analysis,'Frequency')
                                subplot(3,3,6,'replace');
                                topoplot(toplot(:,frame),LIMO.data.chanlocs);
                                title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
                                
                                subplot(3,3,9,'replace');
                                plot(freqvect,toplot(y,:),'LineWidth',3); grid on; axis tight
                                if strcmp(LIMO.Type,'Components')
                                    mytitle2 = sprintf('power spectra @ \n component %g', y);
                                elseif strcmp(LIMO.Type,'Channels')
                                    mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                                end
                                title(mytitle2,'FontSize',12);
                            end
                        end
                    end
                end
                
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
        
        % plots
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
    
    
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % -------------------                               LEVEL 2
    % -------------------                            GROUP EFFECTS
    % ------------------------------------------------------------------------------------------------------------------------------------------
    % ------------------------------------------------------------------------------------------------------------------------------------------
    
    
    
elseif LIMO.Level == 2
    
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
        data_cached = 0;
    end
    
    % ------------
    %      Image
    % -----------
    
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
        
        ax(1) = subplot(3,3,[1 2 4 5 7 8]);
        if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'ITC')
            timevect = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
            imagesc(timevect,1:size(toplot,1),scale);
        else
            freqvect=linspace(LIMO.data.freqlist(1),LIMO.data.freqlist(end),size(toplot,2));
            imagesc(freqvect,1:size(M,1),scale);
        end
        title(mytitle,'FontSize',16);
        color_images_(scale,LIMO);
        set(gca,'layer','top');
        ratio = (LIMO.data.end - LIMO.data.start) / size(toplot,2);
        if LIMO.data.start < 0
            frame_zeros = round(abs(LIMO.data.start / ratio));
        else
            frame_zeros = 1;
        end
        
        if size(toplot,1)>1
            ax(2) = subplot(3,3,6);
            if length(f)>0
                f=f(1);
            end
            if strcmp(LIMO.Analysis,'Time')
                topoplot(toplot(:,f),LIMO.data.chanlocs);
                title(['topoplot @' num2str(timevect(f)) 'ms'],'FontSize',12)
            elseif strcmp(LIMO.Analysis,'Frequency')
                topoplot(toplot(:,f),LIMO.data.chanlocs);
                title(['topoplot @' num2str(freqvect(f)) 'Hz'],'FontSize',12)
            end
        end
        
        ax(3) = subplot(3,3,9);
        if length(e)>0
            e=e(1);
        end
        
        if strcmp(LIMO.Analysis,'Time')
            plot(timevect,toplot(e,:),'LineWidth',3); grid on; axis tight
            if size(toplot,1)>1
                mytitle = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,LIMO.data.chanlocs(e).urchan);
            else
                mytitle = sprintf('time course of \n optimized electrode');
            end
        elseif strcmp(LIMO.Analysis,'Frequency')
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
                    if strcmp(LIMO.Analysis,'Time')
                        title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',12)
                    elseif strcmp(LIMO.Analysis,'Frequency')
                        title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
                    end
                    % ERP/Power
                    subplot(3,3,9,'replace'); y = round(y);
                    if strcmp(LIMO.Analysis,'Time')
                        plot(timevect,toplot(y,:),'LineWidth',3); grid on, axis tight
                        if size(toplot,1)>1
                            mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,LIMO.data.chanlocs(y).urchan);
                        else
                            mytitle2 = sprintf('time course of \n optimized electrode');
                        end
                    elseif strcmp(LIMO.Analysis,'Frequency')
                        plot(freqvect,toplot(y,:),'LineWidth',3); grid on, axis tight
                        if size(toplot,1)>1
                            mytitle2 = sprintf('frequencies @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,LIMO.data.chanlocs(y).urchan);
                        else
                            mytitle2 = sprintf('frequencies of \n optimized electrode');
                        end
                    end
                    title(mytitle2,'FontSize',12)
                    clear x y button
                end
            end
        end
        
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


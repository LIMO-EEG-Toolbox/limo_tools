function limo_display_image(LIMO,toplot,mask,mytitle)

% This function displays images with a intensity plotted as function of
% time or frequency (x) and electrodes (y) - for ERSP it precomputes what
% needs to be plotted and call limo_display_image_tf
%
% FORMAT: limo_display_image(FileName,PathName,p,MCC,LIMO,flag)
%
% INPUTS:
%   LIMO.mat  = Name of the file to image
%   toplot    = 2D matrix to plot (typically t/F values)
%   mask      = areas for which to show data (to show all mask = ones(size(topolot))
%   MCC       = Multiple Comparison technique
%               1=None, 2= Cluster, 3=TFCE, 4=T max
%
% the function originates from previous version of limo_display_results
% Cyril Pernet v2 January 2016
% ----------------------------------
%  Copyright (C) LIMO Team 2016

figure; set(gcf,'Color','w','InvertHardCopy','off');

% get some informations for the plots
v = max(toplot(:)); [e,f]=find(toplot==v);
if length(e)>1 % happen if we have multiple times the exact same max values
    e = e(1); f = f(1); % then we take the 1st (usually an artefact but allows to see it)
end

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
    
elseif strcmp(LIMO.Analysis,'Frequency')
    freqvect = LIMO.data.freqlist;
    if size(freqvect,2) == 1; freqvect = freqvect'; end
    if size(freqvect,2) ~= size(toplot,2)
        freqvect = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
    end
    frame_zeros = 1;
    ratio =  (freqvect(end)-freqvect(1)) / length(freqvect);
    scale = toplot.*mask; scale(scale==0)=NaN;
end

% ERP plot at best electrode
ax(3) = subplot(3,3,9);

if strcmp(LIMO.Analysis,'Time')
    if isfield(LIMO,'Type')
        if strcmp(LIMO.Type,'Components')
            mytitle2 = sprintf('time course @ \n component %g', e);
        elseif strcmp(LIMO.Type,'Channels')
            mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
        end
    else
        try
            mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
        catch
            mytitle2 = sprintf('time course @ y=%g', e);
        end
    end
    plot(timevect,toplot(e,:),'LineWidth',3); grid on; axis tight
    
elseif strcmp(LIMO.Analysis,'Frequency')
    if isfield(LIMO,'Type')
        if strcmp(LIMO.Type,'Components')
            mytitle2 = sprintf('power spectra @ \n component %g', e);
        elseif strcmp(LIMO.Type,'Channels')
            mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
        end
    else
        try
            mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
        catch
            mytitle2 = sprintf('power spectra @ y=%g', e);
        end
    end
    plot(freqvect,toplot(e,:),'LineWidth',3); grid on; axis tight
end
title(mytitle2,'FontSize',12)

% topoplot at max time
ax(2) = subplot(3,3,6);
if isfield(LIMO,'Type')
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

% images toplot
ax(1) = subplot(3,3,[1 2 4 5 7 8]);
if strcmp(LIMO.Analysis,'Time')
    imagesc(timevect,1:size(toplot,1),scale);
elseif strcmp(LIMO.Analysis,'Frequency')
    imagesc(freqvect,1:size(toplot,1),scale);
end

try
    caxis([min(min(scale)), max(max(scale))]);
catch caxiserror
end
title(mytitle,'Fontsize',14)
color_images_(scale,LIMO);

% update with mouse clicks
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
            if isfield(LIMO,'Type')
                if strcmp(LIMO.Type,'Components')
                    mytitle2 = sprintf('time course @ \n component %g', y);
                elseif strcmp(LIMO.Type,'Channels')
                    mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                end
            else
                try
                    mytitle2 = sprintf('time course @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                catch
                    mytitle2 = sprintf('time course @ \n y=%g)', y);
                end
            end
            title(mytitle2,'FontSize',12);
            
        elseif strcmp(LIMO.Analysis,'Frequency')
            subplot(3,3,6,'replace');
            topoplot(toplot(:,frame),LIMO.data.chanlocs);
            title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
            
            subplot(3,3,9,'replace');
            plot(freqvect,toplot(y,:),'LineWidth',3); grid on; axis tight
            if isfield(LIMO,'Type')
                if strcmp(LIMO.Type,'Components')
                    mytitle2 = sprintf('power spectra @ \n component %g', y);
                elseif strcmp(LIMO.Type,'Channels')
                    mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                end
            else
                try
                    mytitle2 = sprintf('power spectra @ \n electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                catch
                    mytitle2 = sprintf('power spectra @ \n y=%g)', y);
                end
                title(mytitle2,'FontSize',12);
            end
        end
    end
    color_images_(scale,LIMO);
end
    
    %% color map
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    function color_images_(scale,LIMO)
    
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];colormap(cc);
    % limo_colormap
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
end


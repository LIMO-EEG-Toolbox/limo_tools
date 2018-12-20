function LIMO_display_image(LIMO,toplot,mask,mytitle,dynamic)

% This function displays images with a intensity plotted as function of
% time or frequency (x) and electrodes (y) - for ERSP it precomputes what
% needs to be plotted and call LIMO_display_image_tf
%
% FORMAT: LIMO_display_image(LIMO,toplot,mask,mytitle,dynamic)
%
% INPUTS:
%   LIMO.mat  = Name of the file to image
%   toplot    = 2D matrix to plot (typically t/F values)
%   mask      = areas for which to show data (to show all mask = ones(size(topolot))
%   mytitle   = title to show
%   dynamic   = set to 0 for no interaction (default is 1)
%
% the function originates from previous version of LIMO_display_results
% Cyril Pernet v2 January 2016
% ----------------------------------
%  Copyright (C) LIMO Team 2016

if nargin == 4
    dynamic = 1;
end

%% get some informations for the plots

v = max(toplot(:));      % from the 2D data to plot, find max
[e,f]=find(toplot==v);   % which channel and time/frequency frame
if length(e)>1           % if we have multiple times the exact same max values
    e = e(1); f = f(1);  % then take the 1st (usually an artefact but allows to see it)
end

%% what do we plot? 

scale = toplot.*mask;  % the data masked by the mask (tpically of significance)
scale(scale==0)=NaN;   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             ERP            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(LIMO.Analysis,'Time')
    try timevect = LIMO.data.timevect; catch timevect = []; end       % precomputed time info
    if size(timevect,2) == 1; timevect = timevect'; end
    if size(timevect,2) ~= size(toplot,2)
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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Spectrum            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(LIMO.Analysis,'Frequency')
    freqvect = LIMO.data.freqlist;
    if size(freqvect,2) == 1; freqvect = freqvect'; end
    if size(freqvect,2) ~= size(toplot,2)
        freqvect = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
    end
    frame_zeros = 1;
    ratio =  (freqvect(end)-freqvect(1)) / length(freqvect);
end


%% make the main figure
% ----------------------------------------------------------------------------
figure; set(gcf,'Color','w','InvertHardCopy','off');

% course plot at best electrode
ax(3) = subplot(3,3,9);
if ~isfield(LIMO.data, 'chanlocs') || isfield(LIMO.data,'expected_chanlocs')
    LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
end

if size(toplot,1) == 1
    plot(timevect,toplot); grid on; ylabel('stat value'); axis tight
    if isfield(LIMO,'Type')
        if strcmp(LIMO.Type,'Components')
            mytitle2 = 'Average component';
        elseif strcmp(LIMO.Type,'Channels')
            mytitle2 = 'Average electrode';
        end
    else
        mytitle2 = 'Average electrode';
    end
else
    if strcmp(LIMO.Analysis,'Time')
        if isfield(LIMO,'Type')
            if strcmp(LIMO.Type,'Components')
                mytitle2 = sprintf('time course @ \n component %g', e);
            elseif strcmp(LIMO.Type,'Channels')
                label = LIMO.data.chanlocs(e).labels;
                mytitle2 = sprintf('time course @ \n electrode %s (%g)', label,e);
            end
        else
            try
                label = LIMO.data.chanlocs(e).labels;
                mytitle2 = sprintf('time course @ \n electrode %s (%g)', label.labels,e);
            catch
                mytitle2 = sprintf('time course @ y=%g', e);
            end
        end
        plot(timevect,toplot(e,:),'LineWidth',3); grid on; axis tight
        
    elseif strcmp(LIMO.Analysis,'Frequency')
        if isfield(LIMO,'Type')
            if strcmp(LIMO.Type,'Components')
                mytitle2 = sprintf('power spectrum @ \n component %g', e);
            elseif strcmp(LIMO.Type,'Channels')
                mytitle2 = sprintf('power spectrum @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
            end
        else
            try
                mytitle2 = sprintf('power spectrum @ \n electrode %s (%g)', LIMO.data.chanlocs(e).labels,e);
            catch
                mytitle2 = sprintf('power spectrum @ y=%g', e);
            end
        end
        plot(freqvect,toplot(e,:),'LineWidth',3); grid on; axis tight
    end
end
title(mytitle2,'FontSize',12)

% topoplot at max time
% ---------------------
if size(toplot,1) ~= 1
    if isempty(findstr(LIMO.design.name, ['one ' LIMO.Type(1:end-1)])) && ~isempty(LIMO.data.chanlocs)
        
        ax(2) = subplot(3,3,6);
        chans = LIMO.data.chanlocs;
        opt = {'maplimits','maxmin','verbose','off'};
        
        if isfield(LIMO,'Type')
            if strcmp(LIMO.Type,'Components')
                opt = {'maplimits','absmax','electrodes','off','verbose','off'};
                topoplot(toplot(:,f),chans,opt{:});
            else
                topoplot(toplot(:,f),chans,opt{:});
            end
            
            if size(toplot,2) == 1
                title('Topoplot','FontSize',12)
            else
                if strcmp(LIMO.Analysis,'Time')
                    title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
                    set(gca,'XTickLabel', timevect);
                elseif strcmp(LIMO.Analysis,'Frequency')
                    title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
                    set(gca,'XTickLabel', LIMO.data.freqlist);
                end
            end
            
        elseif ~isempty(chans)
            topoplot(toplot(:,f),chans,opt{:});
            if size(toplot,2) == 1
                title('Topoplot','FontSize',12)
            else
                if strcmp(LIMO.Analysis,'Time')
                    title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
                    set(gca,'XTickLabel', timevect);
                elseif strcmp(LIMO.Analysis,'Frequency')
                    title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
                    set(gca,'XTickLabel', LIMO.data.freqlist);
                end
            end
        end
    end
end

% images toplot
% -------------------------------
ax(1) = subplot(3,3,[1 2 4 5 7 8]);
if strcmp(LIMO.Analysis,'Time')
    imagesc(timevect,1:size(toplot,1),scale);
elseif strcmp(LIMO.Analysis,'Frequency')
    imagesc(freqvect,1:size(toplot,1),scale);
end

try
    maxval = max(abs(max(scale(:))),abs(min(scale(:))));
    if max(scale(:)) < 0
        caxis([-maxval 0])
    elseif min(scale(:)) > 0 
        caxis([0 maxval])
    else
        caxis([-maxval maxval])
    end
    % caxis([min(min(scale)), max(max(scale))]);
catch caxiserror
end
title(mytitle,'Fontsize',12)
cc = limo_color_images(scale,LIMO);

% ------------------------
% update with mouse clicks
% ------------------------
if dynamic == 1
    if size(toplot,1) > 1
        update = 0;
        while update ==0
            try
                [x,y,button]=ginput(1);
            catch
                update =1; break
            end
            
            if button > 1
                update = 1; % right click to come out of the dynamic figure
            end
            
            clickedAx = gca;
            if clickedAx ~=ax(1)
                disp('right click to exit')
            else
                % topoplot at new time or freq
                frame = frame_zeros + round(x / ratio);
                if frame<=0; frame = 1; end
                if frame>=size(toplot,2); frame=size(toplot,2); end
                
                % course plot at best electrode and 
                y = round(y);
                if size(toplot,1)> 1 && y>size(toplot,1)
                    y = size(toplot,1);
                elseif size(toplot,1)> 1 && y<1
                    y = 1;
                end
            
                
                if strcmp(LIMO.Analysis,'Time') 
                    
                    if isempty(findstr(LIMO.design.name, ['one ' LIMO.Type(1:end-1)])) && ~isempty(LIMO.data.chanlocs)
                        subplot(3,3,6,'replace');
                        if size(toplot,2) == 1
                            topoplot(toplot(:,1),chans,opt{:});
                        else
                            topoplot(toplot(:,frame),chans,opt{:});
                        end
                        if size(toplot,2) == 1
                            title('Topoplot','FontSize',12)
                        else
                            title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',12)
                        end
                    end
                    
                    subplot(3,3,9,'replace');
                    if size(toplot,2) == 1
                        bar(toplot(y,1)); grid on; axis([0 2 0 max(toplot(:))+0.2]); ylabel('stat value')
                        if isfield(LIMO,'Type')
                            if strcmp(LIMO.Type,'Components')
                                mytitle2 = sprintf('component %g', y);
                            elseif strcmp(LIMO.Type,'Channels')
                                mytitle2 = sprintf('Electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                            end
                        else
                            mytitle2 = sprintf('Electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                        end
                    else
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
                    end
                    title(mytitle2,'FontSize',12);
                    
                elseif strcmp(LIMO.Analysis,'Frequency')
                    
                    if isempty(findstr(LIMO.design.name, ['one ' LIMO.Type(1:end-1)])) && ~isempty(LIMO.data.chanlocs)
                        subplot(3,3,6,'replace');
                        topoplot(toplot(:,frame),LIMO.data.chanlocs);
                        if size(toplot,2) == 1
                            title('Topoplot','FontSize',12)
                        else
                            title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
                        end
                    end
                    
                    subplot(3,3,9,'replace');
                    if size(toplot,2) == 1
                        bar(toplot(e,1)); grid on; axis([0 2 0 max(toplot(:))+0.2]);
                        if isfield(LIMO,'Type')
                            if strcmp(LIMO.Type,'Components')
                                mytitle2 = sprintf('component %g', y);
                            elseif strcmp(LIMO.Type,'Channels')
                                mytitle2 = sprintf('Electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                            end
                        else
                            mytitle2 = sprintf('Electrode %s (%g)', LIMO.data.chanlocs(y).labels,y);
                        end
                    else
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
                        end
                        title(mytitle2,'FontSize',12);
                    end
                end
                
            end
            
            colormap(cc); 
            try
                p_values = evalin('base','p_values');
                if ~isnan(p_values(round(y),frame))
                    fprintf('Stat value: %g, p_value %g \n',toplot(round(y),frame),p_values(round(y),frame));
                end
            end
        end
    end
end
end



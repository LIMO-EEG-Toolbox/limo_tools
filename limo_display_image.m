function limo_display_image(LIMO,toplot,mask,mytitle,dynamic)

% This function displays images with a intensity plotted as function of
% time or frequency (x) and electrodes (y) - for ERSP it precomputes what
% needs to be plotted and call LIMO_display_image_tf
%
% FORMAT: limo_display_image(LIMO,toplot,mask,mytitle,dynamic)
%
% INPUTS:
%   LIMO.mat  = Name of the file to image
%   toplot    = 2D matrix to plot (typically t/F values)
%   mask      = areas for which to show data (to show all mask = ones(size(topolot))
%   mytitle   = title to show
%   dynamic   = set to 0 for no interaction (default is 1)
%
% The colour scales are from https://github.com/CPernet/brain_colours
% using linear luminance across the range with cool for negative and 
% hot for positive maps and the divergent BWR scale for negative and positive
% maps. Note that masked values are always gray.
%
% Reference: Pernet & Madan (2019). Data visualization for inference in
% tomographic brain imaging. 
% https://onlinelibrary.wiley.com/doi/full/10.1111/ejn.14430
%
% ----------------------------------
%  Copyright (C) LIMO Team 2019

if nargin == 4
    dynamic = 1;
end

%% get some informations for the plots

% what do we plot?  the data (toplot) masked (tpically of significance)
scale           = toplot.*single(mask>0);  
scale(scale==0) = NaN;   
cc              = limo_color_images(scale); % get a color map commensurate to that

v = max(scale(:));       % from the 2D data to plot, find max
[e,f]=find(scale==v);    % which channel and time/frequency frame
if length(e)>1           % if we have multiple times the exact same max values
    e = e(1); f = f(1);  % then take the 1st (usually an artefact but allows to see it)
end

% for each cluster, get start/end/max value
% if unthresholded, uncorrected, tfce or max = mask is made up of ones
n_cluster     = max(mask(:));
cluster_start = NaN(1,n_cluster); % start of each cluster
cluster_end   = NaN(1,n_cluster); % end of each cluster
cluster_maxv  = NaN(1,n_cluster); % max value for each cluster
cluster_maxe  = NaN(1,n_cluster); % channel location of the max value of each cluster
cluster_maxf  = NaN(1,n_cluster); % frame location of the max value of each cluster

for c=1:n_cluster
    tmp                               = toplot.*(mask==c);
    sigframes                         = sum(tmp,1);
    cluster_start(c)                  = find(sigframes,1,'first');
    cluster_end(c)                    = find(sigframes,1,'last');
    V                                 = max(tmp(:));
    cluster_maxv(c)                   = V(1);
    [cluster_maxe(c),cluster_maxf(c)] = ind2sub(size(tmp),find(tmp==V(1)));
end


%% get frame information 
if strcmpi(LIMO.Analysis,'Time')
    if isfield(LIMO.data,'timevect')
        timevect = LIMO.data.timevect;
        if size(timevect,2) == 1; timevect = timevect'; end
    else
        timevect = [];
    end

    if size(timevect,2) ~= size(toplot,2)
        timevect           = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
        LIMO.data.timevect =  timevect;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
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
    
elseif strcmpi(LIMO.Analysis,'Frequency')
    if isfield(LIMO.data,'freqlist')
        freqvect = LIMO.data.freqlist;
        if size(freqvect,2) == 1; freqvect = freqvect'; end
    else
        freqvect = [];
    end
    
    if size(freqvect,2) ~= size(toplot,2)
        freqvect           = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
        LIMO.data.freqlist = freqvect;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    end
    
    ratio =  (freqvect(end)-freqvect(1)) / length(freqvect);
    if LIMO.data.start < 0
        frame_zeros = find(freqvect == 0);
        if isempty(frame_zeros)
            frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
        end
    else
        frame_zeros = 1;
    end
    
elseif strcmpi(LIMO.Analysis,'Time-Frequency')
    if isfield(LIMO.data,'tf_times')
        timevect = LIMO.data.tf_times;
        if size(timevect,2) == 1; timevect = timevect'; end
    else
        timevect = [];
    end
    
    if size(timevect,2) ~= size(toplot,2)
        timevect           = linspace(LIMO.data.start,LIMO.data.end,size(toplot,2));
        LIMO.data.tf_times =  timevect;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
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
    
    if isfield(LIMO.data,'tf_freqs')
        freqvect = LIMO.data.tf_freqs;
        if size(freqvect,2) == 1; freqvect = freqvect'; end
    else
        freqvect = [];
    end
    
    if size(freqvect,2) ~= size(toplot,1)
        freqvect           = linspace(LIMO.data.lowf,LIMO.data.highf,size(toplot,1));
        LIMO.data.tf_freqs =  freqvect;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    end
    
else
    error('LIMO.Analysis unspecfied')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the main figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; set(gcf,'Color','w','InvertHardCopy','off');

% course plot at best electrode
ax(3) = subplot(3,3,9);
if ~isfield(LIMO.data, 'chanlocs') || isfield(LIMO.data,'expected_chanlocs')
    LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
end

if size(toplot,1) == 1
    if strcmpi(LIMO.Analysis,'Time')
        plot(timevect,toplot,'LineWidth',3); 
    elseif strcmpi(LIMO.Analysis,'Frequency')
        plot(freqvect,toplot,'LineWidth',3);
    end
    grid on; ylabel('stat value'); axis tight
    
    if isfield(LIMO,'Type')
        if strcmpi(LIMO.Type,'Components')
            mytitle2 = 'Average component';
        elseif strcmpi(LIMO.Type,'Channels')
            mytitle2 = 'Average channel';
        end
    else
        mytitle2 = 'Average channel';
    end
else  
    if strcmpi(LIMO.Analysis,'Time')
        plot(timevect,toplot(e,:),'LineWidth',3); 
    elseif strcmpi(LIMO.Analysis,'Frequency')
        plot(freqvect,toplot(e,:),'LineWidth',3);
   elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        plot(timevect,toplot(e,:),'LineWidth',3);
        mytitle2 = sprintf('stat values @ %g Hz', freqvect(e));
    end
    grid on; axis tight
    
    if ~strcmpi(LIMO.Analysis,'Time-Frequency')
        if isfield(LIMO,'Type')
            if strcmpi(LIMO.Type,'Components')
                mytitle2 = sprintf('stat values @ \n component %g', e);
            elseif strcmpi(LIMO.Type,'Channels')
                label = LIMO.data.chanlocs(e).labels;
                mytitle2 = sprintf('stat values @ \n channel %s (%g)', label,e);
            end
        else
            try
                label = LIMO.data.chanlocs(e).labels;
                mytitle2 = sprintf('stat values @ \n channel %s (%g)', label.labels,e);
            catch
                mytitle2 = sprintf('stat values @ y=%g', e);
            end
        end
    end
end
title(mytitle2,'FontSize',12)

% topoplot at max time
% ---------------------
if size(toplot,1) ~= 1 && ~strcmpi(LIMO.Analysis,'Time-Frequency')
    
    ax(2) = subplot(3,3,6);
    opt   = {'maplimits','maxmin','verbose','off','colormap', limo_color_images(toplot)};
    
    if isfield(LIMO,'Type')
        if strcmpi(LIMO.Type,'Components')
            opt = {'maplimits','absmax','electrodes','off','verbose','off','colormap', limo_color_images(toplot)};
            topoplot(toplot(:,f),LIMO.data.chanlocs,opt{:});
        else
            topoplot(toplot(:,f),LIMO.data.chanlocs,opt{:});
        end
        
        if size(toplot,2) == 1
            title('Topoplot','FontSize',12)
        else
            if strcmpi(LIMO.Analysis,'Time')
                title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
                set(gca,'XTickLabel', timevect);
            elseif strcmpi(LIMO.Analysis,'Frequency')
                title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
                set(gca,'XTickLabel', LIMO.data.freqlist);
            end
        end
        
    elseif ~isempty(LIMO.data.chanlocs)
        topoplot(toplot(:,f),LIMO.data.chanlocs,opt{:});
        if size(toplot,2) == 1
            title('Topoplot','FontSize',12)
        else
            if strcmpi(LIMO.Analysis,'Time')
                title(['topoplot @ ' num2str(round(timevect(f))) 'ms'],'FontSize',12)
                set(gca,'XTickLabel', timevect);
            elseif strcmpi(LIMO.Analysis,'Frequency')
                title(['topoplot @' num2str(round(freqvect(f))) 'Hz'],'FontSize',12);
                set(gca,'XTickLabel', LIMO.data.freqlist);
            end
        end
    end
end

% images toplot
% -------------------------------
ax(1) = subplot(3,3,[1 2 4 5 7 8]);
if strcmpi(LIMO.Analysis,'Time')
    imagesc(timevect,1:size(toplot,1),scale);
    colormap(gca, cc);
elseif strcmpi(LIMO.Analysis,'Frequency')
    imagesc(freqvect,1:size(toplot,1),scale);
    colormap(gca, cc);
elseif strcmpi(LIMO.Analysis,'Time-Frequency')
    imagesc(timevect,freqvect,scale);
    colormap(gca, cc);
end
set_imgaxes(LIMO,scale);
title(mytitle,'Fontsize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% return cluster info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning off
if contains(mytitle,'cluster')
    for c=1:n_cluster
        if strcmpi(LIMO.Analysis,'Time')
        fprintf('cluster %g starts at %gms ends at %gms, max %g @ %gms channel %s \n', c, ...
            timevect(cluster_start(c)),timevect(cluster_end(c)), cluster_maxv(c), timevect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
        elseif strcmpi(LIMO.Analysis,'Frequency')
        fprintf('cluster %g starts at %gHz ends at %gHz, max %g @ %gHz channel %s \n', c, ...
            freqvect(cluster_start(c)),freqvect(cluster_end(c)), cluster_maxv(c), freqvect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
        elseif strcmpi(LIMO.Analysis,'Time-Frequency')
            f1 = scale(:,cluster_start(c)); f2 = scale(:,cluster_end(c)); f3 = scale(:,cluster_maxf(c));
            fprintf('cluster %g at %gms %gHz, ends at %gms %gHz, max %g @ %gms %gHz \n', c, ...
                timevect(cluster_start(c)),freqvect(find(f1==max(f1))), ...
                timevect(cluster_end(c)), freqvect(find(f2==max(f2))), ...
                cluster_maxv(c), timevect(cluster_maxf(c)), freqvect(find(f3==max(f3))));
        end
    end
else % no clusters (we have make one )
    if strcmpi(LIMO.Analysis,'Time')
        fprintf('1st significant frame at %gms, last signifiant frame at %gms, max %g @ %gms channel %s \n', ...
            timevect(cluster_start(c)),timevect(cluster_end(c)), cluster_maxv(c), timevect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
    elseif strcmpi(LIMO.Analysis,'Frequency')
        fprintf('1st significant frame at %gHz, last signifiant frame at %gHz, max %g @ %gHz channel %s \n', ...
            freqvect(cluster_start(c)),freqvect(cluster_end(c)), cluster_maxv(c), freqvect(cluster_maxf(c)), LIMO.data.chanlocs(cluster_maxe(c)).labels);
    elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        f1 = scale(:,cluster_start(c)); f2 = scale(:,cluster_end(c)); f3 = scale(:,cluster_maxf(c));
        fprintf('1st significant frame at %gms %gHz, last signifiant frame at %gms %gHz, max %g @ %gms %gHz \n', ...
            timevect(cluster_start(c)),freqvect(find(f1==max(f1))), ...
            timevect(cluster_end(c)), freqvect(find(f2==max(f2))), ...
            cluster_maxv(c), timevect(cluster_maxf(c)), freqvect(find(f3==max(f3))));
    end
end
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update with mouse clicks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dynamic == 1
    if size(toplot,1) > 1
        update = 0;
        while update ==0
            try
                [x,y,button]=ginput(1);
            catch
                break
            end
            
            if button > 1
                update = 1; % right click to come out of the dynamic figure
            end
            
            clickedAx = gca;
            if clickedAx ~=ax(1)
                disp('right click to exit')
            else
                % topoplot at new time or freq
                % because x axis is already with the correct value
                frame = frame_zeros + round(x / ratio);
                if frame<=0; frame = 1; end
                if frame>=size(toplot,2); frame=size(toplot,2); end
                
                if strcmpi(LIMO.Analysis,'Time-Frequency')
                    % course plot at selected frequency
                    y          = round(y); 
                    [~,yframe] = min(abs(freqvect - y));
                    if size(toplot,1)> 1 && yframe>size(toplot,1)
                        yframe = size(toplot,1);
                        y      = freqvect(yframe);
                    elseif size(toplot,1)> 1 && y<1
                        yframe = 1;
                        y      = freqvect(yframe);
                    end
                else
                    % course plot at selected  channel
                    y = round(y);
                    if size(toplot,1)> 1 && y>size(toplot,1)
                        y = size(toplot,1);
                    elseif size(toplot,1)> 1 && y<1
                        y = 1;
                    end
                end
                
                if strcmpi(LIMO.Analysis,'Time')
                    if ~contains(LIMO.design.name, ['one ' LIMO.Type(1:end-1)]) && ~isempty(LIMO.data.chanlocs)
                        subplot(3,3,6,'replace');
                        if size(toplot,2) == 1
                            topoplot(toplot(:,1),LIMO.data.chanlocs,opt{:});
                            title('topoplot','FontSize',12)
                        else
                            topoplot(toplot(:,frame),LIMO.data.chanlocs,opt{:});
                            title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',12)
                        end
                    end
                    
                elseif strcmpi(LIMO.Analysis,'Frequency')
                    if ~contains(LIMO.design.name, ['one ' LIMO.Type(1:end-1)]) && ~isempty(LIMO.data.chanlocs)
                        subplot(3,3,6,'replace');
                        if size(toplot,2) == 1
                            topoplot(toplot(:,1),LIMO.data.chanlocs,opt{:});
                            title('topoplot','FontSize',12)
                        else
                            topoplot(toplot(:,frame),LIMO.data.chanlocs,opt{:});
                            title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
                        end
                    end
                end
                
                subplot(3,3,9,'replace');
                if size(toplot,2) == 1
                    bar(toplot(y,1)); grid on; axis([0 2 0 max(toplot(:))+0.2]); ylabel('stat value')
                    if isfield(LIMO,'Type')
                        if strcmpi(LIMO.Type,'Components')
                            mytitle2 = sprintf('component %g', y);
                        elseif strcmpi(LIMO.Type,'Channels')
                            mytitle2 = sprintf('channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
                        end
                    else
                        mytitle2 = sprintf('channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
                    end
                else
                    if strcmpi(LIMO.Analysis,'Time')
                        plot(timevect,toplot(y,:),'LineWidth',3);
                    elseif strcmpi(LIMO.Analysis,'Frequency')
                        plot(freqvect,toplot(y,:),'LineWidth',3);
                    elseif strcmpi(LIMO.Analysis,'Time-Frequency')
                        plot(timevect,toplot(yframe,:),'LineWidth',3);
                        mytitle2 = sprintf('stat values @ %g Hz', y);
                    end
                    grid on; axis tight
                    
                    if ~strcmpi(LIMO.Analysis,'Time-Frequency')
                        if isfield(LIMO,'Type')
                            if strcmpi(LIMO.Type,'Components')
                                mytitle2 = sprintf('stat values @ \n component %g', y);
                            elseif strcmpi(LIMO.Type,'Channels')
                                mytitle2 = sprintf('stat values @ \n channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
                            end
                        else
                            try
                                mytitle2 = sprintf('stat values @ \n channel %s (%g)', LIMO.data.chanlocs(y).labels,y);
                            catch
                                mytitle2 = sprintf('stat values @ \n y=%g)', y);
                            end
                        end
                    end
                end
                title(mytitle2,'FontSize',12);

                subplot(3,3,[1 2 4 5 7 8]); colormap(gca, cc);
                
                try
                    p_values = evalin('base','p_values');
                    if strcmpi(LIMO.Analysis,'Time-Frequency')
                        if ~isnan(p_values(yframe,frame))
                            fprintf('Stat value: %g, p_value %g \n',toplot(yframe,frame),p_values(yframe,frame));
                        end
                    else
                        if ~isnan(p_values(round(y),frame))
                            fprintf('Stat value: %g, p_value %g \n',toplot(round(y),frame),p_values(round(y),frame));
                        end
                    end
                catch pvalerror
                    fprintf('couldn''t figure the p value?? %s \n',pvalerror.message)
                end
            end
        end
    end
end
end

%% set axes and labels 
% -------------------------------------------------------------------------
function set_imgaxes(LIMO,scale)

img_prop = get(gca);
set(gca,'LineWidth',2)

% ----- X --------
if strcmpi(LIMO.Analysis,'Time') || strcmpi(LIMO.Analysis,'Time-Frequency')
    xlabel('Time in ms','FontSize',10)
elseif strcmpi(LIMO.Analysis,'Frequency')
    xlabel('Frequency in Hz','FontSize',10)
end
 
% ----- Y --------
if strcmpi(LIMO.Analysis,'Time-Frequency')
    ylabel('Frequency in Hz','FontSize',10)
else
    if strcmpi(LIMO.Type,'Components')
        if size(scale,1) == 1
            ylabel('Optimized component','FontSize',10);
        else
            ylabel('Components','FontSize',10);
        end
    else
        if size(scale,1) == 1
            ylabel('Optimized channel','FontSize',10);
        else
            ylabel('Channels','FontSize',10);
        end
    end
    
    if isfield(LIMO.data, 'chanlocs')
        Ylabels = arrayfun(@(x)(x.labels), LIMO.data.chanlocs, 'UniformOutput', false);
    else
        Ylabels = arrayfun(@(x)(x.labels), LIMO.data.expected_chanlocs, 'UniformOutput', false);
    end
    newticks = round(linspace(1,length(Ylabels),length(img_prop.YTick)*2));
    Ylabels  = Ylabels(newticks);
    if size(scale,1) == 1
        set(gca,'YTick',1);
    else
        set(gca,'YTick',newticks);
        set(gca,'YTickLabel', Ylabels);
    end
end

% ----- Colormap --------
try
    maxval = max(abs(max(scale(:))),abs(min(scale(:))));
    if max(scale(:)) < 0
        caxis([-maxval 0])
    elseif min(scale(:)) > 0 
        caxis([0 maxval])
    else
        caxis([-maxval maxval])
    end
catch caxiserror
    fprintf('axis issue: %s\n',caxiserror.message)
end

end


% Callback function for limo_display_image

function limo_display_image_callback(fig,obj,x,y)

udat = get(fig, 'userdata');
LIMO = udat.LIMO;

if ~isempty(obj)
    if isequal(get(obj, 'Type'), 'axes')
        tmppos = get(obj, 'currentpoint');
    else
        tmppos = get(obj.Parent, 'currentpoint');
    end
    x = tmppos(1,1);
    y = tmppos(1,2);
end

% topoplot at new time or freq
% because x axis is already with the correct value

frame = udat.frame_zeros + round(x / udat.ratio);

if strcmpi(LIMO.Analysis,'Time-Frequency')
    % course plot at selected frequency
    y          = round(y);
    [~,yframe] = min(abs(freqvect - y));
    if size(udat.toplot,1)> 1 && yframe>size(udat.toplot,1)
        yframe = size(udat.toplot,1);
        y      = freqvect(yframe);
    elseif size(udat.toplot,1)> 1 && y<1
        yframe = 1;
        y      = freqvect(yframe);
    end
else
    % course plot at selected  channel
    y = round(y);
    if size(udat.toplot,1)> 1 && y>size(udat.toplot,1)
        y = size(udat.toplot,1);
    elseif size(udat.toplot,1)> 1 && y<1
        y = 1;
    end
end

if strcmpi(LIMO.Analysis,'Time')
    if ~contains(LIMO.design.name, ['one ' LIMO.Type(1:end-1)]) && ~isempty(LIMO.data.chanlocs)
        subplot(3,3,3,'replace');
        if size(udat.toplot,2) == 1
            topoplot(udat.toplot(:,1),LIMO.data.chanlocs,udat.opt{:});
            title('topoplot','FontSize',12)
        else
            topoplot(udat.toplot(:,frame),LIMO.data.chanlocs,udat.opt{:});
            title(['topoplot @ ' num2str(round(x)) 'ms'],'FontSize',12)
        end
    end

elseif strcmpi(LIMO.Analysis,'Frequency')
    if ~contains(LIMO.design.name, ['one ' LIMO.Type(1:end-1)]) && ~isempty(LIMO.data.chanlocs)
        subplot(3,3,3,'replace');
        if size(udat.toplot,2) == 1
            topoplot(udat.toplot(:,1),LIMO.data.chanlocs,udat.opt{:});
            title('topoplot','FontSize',12)
        else
            topoplot(udat.toplot(:,frame),LIMO.data.chanlocs,udat.opt{:});
            title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
        end
    end
end

ax = subplot(3,3,6,'replace');
pos = get(ax, 'position');
set(ax, 'position', [ pos(1) pos(2)+0.07 pos(3) pos(4)-0.1])
if size(udat.toplot,2) == 1
    bar(udat.toplot(y,1)); grid on; axis([0 2 0 max(udat.toplot(:))+0.2]); ylabel('stat value')
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
        plot(udat.timevect,udat.toplot(y,:),'LineWidth',3);
    elseif strcmpi(LIMO.Analysis,'Frequency')
        plot(udat.freqvect,udat.toplot(y,:),'LineWidth',3);
    elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        plot(udat.timevect,udat.toplot(yframe,:),'LineWidth',3);
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

subplot(3,3,[1 2 4 5 7 8],'replace'); 
cla;
h_im = [];
if strcmpi(LIMO.Analysis,'Time')
    h_im = imagesc(udat.timevect,1:size(udat.toplot,1),udat.scale);
elseif strcmpi(LIMO.Analysis,'Frequency')
    h_im = imagesc(udat.freqvect,1:size(udat.toplot,1),udat.scale);
elseif strcmpi(LIMO.Analysis,'Time-Frequency')
    h_im = imagesc(udat.timevect,udat.freqvect,udat.scale);
end
colormap(gca, udat.cc);
set_imgaxes(LIMO,udat.scale);
title(udat.mytitle,'Fontsize',12)

% circled region
% -----------
hold on;
if strcmpi(LIMO.Analysis,'Time')
    h = plot(x, y, 'o');
elseif strcmpi(LIMO.Analysis,'Frequency')
    h = plot(x, y, 'o');
elseif strcmpi(LIMO.Analysis,'Time-Frequency')
    % ?
end
set(h,'MarkerSize', 20, 'color', 'k', 'LineWidth', 3);

% stats text
% -----------
strStat = '';
try
    p_values = evalin('base','p_values');
    if strcmpi(LIMO.Analysis,'Time-Frequency')
        if ~isnan(p_values(yframe,frame))
            strStat = sprintf('Stat value: %g, p value: %g',udat.toplot(yframe,frame),p_values(yframe,frame));
        else
            strStat = sprintf('Stat value: %g, p value: NaN?',udat.toplot(yframe,frame));
        end
    else
        if ~isnan(p_values(round(y),frame))
            strStat = sprintf('Stat value: %g, p value: %g',udat.toplot(round(y),frame),p_values(round(y),frame));
        else
            strStat = sprintf('Stat value: %g, p value: NaN?',udat.toplot(round(y),frame));
        end
    end
catch pvalerror
    fprintf('couldn''t figure the stats values?? %s \n',pvalerror.message)
end
fprintf('%s\n', strStat)
ax = subplot(3,3,9,'replace');
pos = get(ax, 'position');
set(ax, 'position', [ pos(1) pos(2)+0.08 pos(3) pos(4)-0.04])
if ~isempty(strStat)
    commaPos = find(strStat == ',');
    strStat1 = strStat(1:commaPos-1);
    strStat2 = strStat(commaPos+2:end);
    text(1,2.3,strStat2, 'Interpreter','none', 'fontweight', 'bold');
    text(1,3,strStat1, 'Interpreter','none', 'fontweight', 'bold');
end
if strcmpi(LIMO.Analysis,'Frequency')
    helptxt = ['Select frequency/electrode' 10 'by cliking on the image' ];
else
    helptxt = ['Select time/electrode' 10 'by cliking on the image' ];
end
text(1,1, helptxt, 'Interpreter','none');
xlim([1 10]);
ylim([-0.5 2.5]);
axis off;

% interactivity
% -----------
set(h_im , 'ButtonDownFcn', 'limo_display_image_callback(gcbf, gcbo)')

if ~isempty(findobj(gcf, 'tag', 'pval'))
    return
else
    udat = get(gcf, 'userdata');
    if ~isempty(udat.params)
        cb = [ 'gcbf2 = gcbf; uDat = get(gcbf2, ''userdata'');' ...
            'pvalTmp = str2num(get(findobj(gcbf2, ''tag'', ''pval''), ''string''));' ...
            'mccTmp  = get(findobj(gcbf2, ''tag'', ''mcc''), ''value'');' ...
            'limo_display_results(uDat.params.Type, uDat.params.FileName, uDat.params.PathName, pvalTmp, mccTmp, uDat.params.LIMO);' ...
            'clear uDat pvalTmp mccTmp; close(gcbf2);' ];
        ui_p = uipanel('Title', 'Masking/statistics', 'BackgroundColor', [.66 .76 1], 'position', [0.675 0.03 0.3 0.18 ]);
        uicontrol(ui_p, 'style', 'text', 'string', 'p<','unit', 'normalized','position', [0.1 0.6 0.2 0.3], 'backgroundcolor', [.66 .76 1]); %, [0.73 0.12 0.05 0.05],
        uicontrol(ui_p, 'style', 'edit', 'string', num2str(udat.params.p),'tag','pval','unit', 'normalized','position', [0.3 0.6 0.3 0.3], 'callback', cb); %, [0.78 0.12 0.08 0.05])
        options = { 'Uncorrected threshold' 'Cluster correction' 'TFCE correction' 'MAX correction' };
        uicontrol(ui_p, 'style', 'popupmenu', 'string', options, 'value', udat.params.MCC, 'tag', 'mcc', 'unit', 'normalized','position',  [0.05 0.25 0.9 0.2], 'callback', cb); %[0.68 0.05 0.25 0.05])
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
    newticks = unique(newticks);
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


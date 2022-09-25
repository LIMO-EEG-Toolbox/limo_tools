% Callback function for limo_display_image

function limo_display_image_callback(fig,obj)

udat = get(fig, 'userdata');
LIMO = udat.LIMO;

if ~isempty(obj)
    y              =         get(findobj(fig, 'tag', 'channel'), 'value');
    x              = str2num(get(findobj(fig, 'tag', 'time'), 'string'));
    udat.timerange = str2num(get(findobj(fig, 'tag', 'timerange'), 'string'));
    udat.colorlim  = str2num(get(findobj(fig, 'tag', 'colorlim'), 'string'));
    if length(udat.colorlim ) ~= 2 ,udat.colorlim  = []; end
    if length(udat.timerange) ~= 2, udat.timerange = []; end

    if isequal(get(obj.Parent, 'Type'), 'axes')
        tmppos = get(obj.Parent, 'currentpoint');
        x = tmppos(1,1);
        y = tmppos(1,2);
    elseif strcmpi(get(obj, 'tag'), 'channel')
        y = get(obj, 'value');
    end

    % save
    udat.x = x;
    udat.y = y;
end

% ----- Colormap --------
maxval = max(abs(max(udat.scale(:))),abs(min(udat.scale(:))));
if isempty(udat.colorlim) 
%     || min(abs(udat.colorlim)) < maxval
%     if ~isempty(udat.colorlim) && min(abs(udat.colorlim)) < maxval
%         fprintf(2, 'Color limits changed so clipped data is displayed properly');
%     end
    if max(udat.scale(:)) < 0
        udat.colorlim = [-maxval 0];
    elseif min(udat.scale(:)) > 0 
        udat.colorlim = [0 maxval];
    else
        udat.colorlim = [-maxval maxval];
    end
end

% topoplot at new time or freq
% because x axis is already with the correct value
if isempty(udat.timerange)
    udat.timerange = udat.timevect([1 end]);
end
x = udat.x;
y = udat.y;
frame = udat.frame_zeros + round(x / udat.ratio);

if isempty(udat.timerange)
    frameRange = 1:length(udat.timevect);
else
    [~,frameBeg] = min(abs(udat.timevect-udat.timerange(1)));
    [~,frameEnd] = min(abs(udat.timevect-udat.timerange(2)));
    frameRange = frameBeg:frameEnd;
end

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

% topoplot 
% --------
if length(LIMO.data.chanlocs) > 2
    ax = subplot(3,4,3,'replace');
    if strcmpi(LIMO.Analysis,'Time')
        if ~contains(LIMO.design.name, ['one ' LIMO.Type(1:end-1)]) && ~isempty(LIMO.data.chanlocs)
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
            if size(udat.toplot,2) == 1
                topoplot(udat.toplot(:,1),LIMO.data.chanlocs,udat.opt{:});
                title('topoplot','FontSize',12)
            else
                topoplot(udat.toplot(:,frame),LIMO.data.chanlocs,udat.opt{:});
                title(['topoplot @ ' num2str(round(x)) 'Hz'],'FontSize',12)
            end
        end
    end
    if isempty(udat.colorlim)
        udat.colorlim = clim;
    end
    clim(udat.colorlim);
    %axcopy(ax);
end

% curve 
% -----
ax = subplot(3,4,7,'replace');
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
        plot(udat.timevect(frameRange),udat.toplot(y,frameRange),'LineWidth',3);
    elseif strcmpi(LIMO.Analysis,'Frequency')
        plot(udat.freqvect(frameRange),udat.toplot(y,frameRange),'LineWidth',3);
    elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        plot(udat.timevect(frameRange),udat.toplot(yframe,frameRange),'LineWidth',3);
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
if ~isempty(udat.colorlim)
    ylim(udat.colorlim);
end
title(mytitle2,'FontSize',12);
%axcopy(ax);

% image 
% -----
subplot(3,4,[1 2 5 6 9 10],'replace'); 
cla;
h_im = [];
if strcmpi(LIMO.Analysis,'Time')
    h_im = imagesc(udat.timevect(frameRange),1:size(udat.toplot,1),udat.scale(:,frameRange));
elseif strcmpi(LIMO.Analysis,'Frequency')
    h_im = imagesc(udat.freqvect(frameRange),1:size(udat.toplot,1),udat.scale(:,frameRange));
elseif strcmpi(LIMO.Analysis,'Time-Frequency')
    h_im = imagesc(udat.timevect(frameRange),udat.freqvect,udat.scale(:,frameRange));
end
colormap(gca, udat.cc);
set_imgaxes(LIMO,udat.scale);
if ~isempty(udat.colorlim)
    clim(udat.colorlim);
end
title(udat.title,'Fontsize',12)

% circled region
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
ax = subplot(3,4,11,'replace');
if ~isempty(strStat)
    commaPos = find(strStat == ',');
    strStat1 = strStat(1:commaPos-1);
    strStat2 = strStat(commaPos+2:end);
    text(0,2.3,strStat2, 'Interpreter','none', 'fontweight', 'bold');
    text(0,3,strStat1, 'Interpreter','none', 'fontweight', 'bold');
end
if strcmpi(LIMO.Analysis,'Frequency')
    helptxt = ['Select freq./electrode' 10 'by cliking on the image' ];
else
    helptxt = ['Select time/electrode' 10 'by cliking on the image' ];
end
text(0,1, helptxt, 'Interpreter','none');
xlim([1 10]);
ylim([0 2.5]);
axis off;

% interactivity
% -----------
cb_redraw     = 'limo_display_image_callback(gcbf, gcbo)';
set(h_im , 'ButtonDownFcn', cb_redraw, 'tag', 'image')
set(fig, 'userdata', udat);

climStr = [num2str(udat.colorlim(1),2) ' ' num2str(udat.colorlim(2),2)];
if ~isempty(findobj(gcf, 'tag', 'pval'))
    set(findobj(gcf, 'tag', 'channel')  , 'value', y);
    set(findobj(gcf, 'tag', 'time')     , 'string', num2str(round(x)));
    set(findobj(gcf, 'tag', 'colorlim') , 'string', climStr);
    return
else
    udat = get(gcf, 'userdata');
    if ~isempty(udat.params)
        cb = [ ...
            'gcbf2 = gcbf; uDat = get(gcbf2, ''userdata'');' ...
            'nopTmp    = get(findobj(gcbf2, ''tag'', ''nomask''), ''value'');' ...
            'pvalTmp   = str2num(get(findobj(gcbf2, ''tag'', ''pval''), ''string''));' ...
            'mccTmp    = get(findobj(gcbf2, ''tag'', ''mcc''), ''value'');' ...
            'imgTmp    = findobj(gcbf2, ''tag'', ''image'');' ...
            'plot3type = get(findobj(gcbf2, ''tag'', ''model''), ''value'');' ...
            'regressor = get(findobj(gcbf2, ''tag'', ''regressor''), ''value'');' ...
            'if nopTmp, pvalTmp = 1; else pvalTmp = 0.05; end;' ...
            'set(imgTmp, ''CData'', get(imgTmp, ''CData'')*0);' ...
            'tmpRes = limo_display_results(''type'', uDat.params.Type, ''filename'', uDat.params.FileName, ''pathname'', uDat.params.PathName, ''p'', pvalTmp, ''MCC'', mccTmp, ''LIMO'', uDat.params.LIMO, ''fig'', gcbf, ''plot3type'', plot3type, ''regressor'', regressor-1);' ...
            'clear uDat nopTmp pvalTmp mccTmp imgTmp tmpRes plot3type regressor;' ];
%            'tmpRes = limo_display_results(''type'', uDat.params.Type, ''filename'', uDat.params.FileName, ''pathname'', uDat.params.PathName, ''p'', pvalTmp, ''MCC'', mccTmp, ''LIMO'', uDat.params.LIMO, ''fig'', gcbf);' ...
        
        % plot what?
        add = 0.03;
        posPannels = { [0.74 0.550+add 0.23 0.4*0.9 ] ...
                       [0.74 0.345+add 0.23 0.2*0.9 ] ...
                       [0.74 0.050+add 0.23 0.3*0.9 ] };
        ui_w = uipanel('Title', 'Plot what?'        , 'BackgroundColor', [.66 .76 1], 'position', posPannels{1});
        ui_l = uipanel('Title', 'Limits'            , 'BackgroundColor', [.66 .76 1], 'position', posPannels{2});
        ui_m = uipanel('Title', 'Masking/statistics', 'BackgroundColor', [.66 .76 1], 'position', posPannels{3});

        % data selection
        add = -0.12;
        pos = { [0.05  0.75+add 0.9 0.3] ...
                [0.05  0.50+add 0.9 0.3] ...
                [0.05  0.25+add 0.9 0.3] [0.45  0.32+add 0.52 0.25] ...
                [0.05  0.00+add 0.9 0.3] [0.72  0.15+add 0.25 0.17] };
        popupData = { 'Statistics' };
        if isfield(LIMO.data, 'chanlocs') popupValues  = { LIMO.data.chanlocs.labels }; else popupValues  = { LIMO.data.expected_chanlocs.labels }; end 
        if isfield(LIMO.design, 'labels') popupData   = { popupData{:} LIMO.design.labels.description }; end 
        %if ~isempty(popupData) && contains(popupData{end}, 'onstant') popupData(end) = []; end
        if isempty(udat.regressor) udat.regressor = 0; end
        if isempty(udat.plot3type) udat.plot3type = 1; end
        popupModel  = { 'Original data' 'Modeled data' 'Ajusted data' };
        enableModel = fastif(udat.regressor == 0, 'off', 'on');
        uicontrol(ui_w, 'style', 'popupmenu',  'string', popupData,    'unit', 'normalized','position', pos{1}, 'tag','regressor','callback', cb, 'value', udat.regressor+1); 
        uicontrol(ui_w, 'style', 'popupmenu',  'string', popupModel ,  'unit', 'normalized','position', pos{2}, 'tag','model','callback', cb, 'value', udat.plot3type, 'enable', enableModel); 
        uicontrol(ui_w, 'style', 'text',       'string', 'Channel',    'unit', 'normalized','position', pos{3}, 'backgroundcolor', [.66 .76 1], 'horizontalalignment', 'left');
        uicontrol(ui_w, 'style', 'popupmenu',  'string', popupValues,  'unit', 'normalized','position', pos{4}, 'tag','channel','callback', cb_redraw, 'value', y); 
        uicontrol(ui_w, 'style', 'text',       'string', 'Time (millisec)' , 'unit', 'normalized','position', pos{5}, 'backgroundcolor', [.66 .76 1], 'horizontalalignment', 'left');
        uicontrol(ui_w, 'style', 'edit',       'string', num2str(round(x)),  'unit', 'normalized','position', pos{6}, 'tag','time','callback', cb_redraw); 

        % limits
        pos = { [0.05 0.5  0.6 0.4] [0.55 0.55  0.42 0.35] ...
                [0.05 0.1  0.6 0.4] [0.55 0.15  0.42 0.35] };
        tlimStr = [num2str(round(udat.timerange(1))) ' ' num2str(round(udat.timerange(end)))];
        uicontrol(ui_l, 'style', 'text',       'string', 'Time range' ,'unit', 'normalized','position', pos{1}, 'backgroundcolor', [.66 .76 1], 'horizontalalignment', 'left');
        uicontrol(ui_l, 'style', 'edit',       'string', tlimStr,      'unit', 'normalized','position', pos{2}, 'tag','timerange','callback', cb_redraw); 
        uicontrol(ui_l, 'style', 'text',       'string', 'Color lim.' ,'unit', 'normalized','position', pos{3}, 'backgroundcolor', [.66 .76 1], 'horizontalalignment', 'left');
        uicontrol(ui_l, 'style', 'edit',       'string', climStr,      'unit', 'normalized','position', pos{4}, 'tag','colorlim','callback', cb_redraw); 

        % masking
        pos = { [0.1  0.7 0.9 0.3] ...
                [0.15 0.35  0.2 0.3] [0.3 0.4 0.3 0.25] ...
                [0.05 0.1 0.9 0.2] };
        pEnable = fastif(udat.params.p == 1, 'off', 'on');
        uicontrol(ui_m, 'style', 'checkbox', 'string', 'No masking','unit', 'normalized','position', pos{1}, 'backgroundcolor', [.66 .76 1], 'tag', 'nomask', 'callback', cb, 'value', udat.params.p == 1); %, [0.73 0.12 0.05 0.05],
        uicontrol(ui_m, 'style', 'text', 'string', 'p<','unit', 'normalized','position', pos{2}, 'backgroundcolor', [.66 .76 1], 'enable', pEnable, 'horizontalalignment', 'left'); %, [0.73 0.12 0.05 0.05],
        uicontrol(ui_m, 'style', 'edit', 'string', num2str(udat.params.p),'tag','pval','unit', 'normalized','position', pos{3}, 'callback', cb, 'enable', pEnable); %, [0.78 0.12 0.08 0.05])
        options = { 'Uncorrected threshold' 'Cluster correction' 'TFCE correction' 'MAX correction' };
        uicontrol(ui_m, 'style', 'popupmenu', 'string', options, 'value', udat.params.MCC, 'tag', 'mcc', 'unit', 'normalized','position',  pos{4}, 'callback', cb, 'enable', pEnable); %[0.68 0.05 0.25 0.05])
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



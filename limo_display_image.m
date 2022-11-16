function res = limo_display_image(LIMO,toplot,mask,mytitle,params,varargin)

% This function displays images with a intensity plotted as function of
% time or frequency (x) and electrodes (y) - for ERSP it precomputes what
% needs to be plotted and call LIMO_display_image_tf
%
% FORMAT: limo_display_image(LIMO,toplot,mask,title,dynamic)
%
% INPUTS:
%   LIMO.mat  = Name of the file to image
%   toplot    = 2D matrix to plot (typically t/F values)
%   mask      = areas for which to show data (to show all mask = ones(size(topolot))
%   title     = title to show
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

if ~ischar(LIMO)
    options = { 'LIMO', LIMO, 'toplot', toplot, 'mask', mask, 'title', mytitle, 'params', params, varargin{:} };
else 
    options = { LIMO,toplot,mask,mytitle,params varargin{:} };
end

try
    if ~isempty( options )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    end
catch
    error('limo_display_image() error: calling convention {''key'', value, ... } error'); return;
end

try g.LIMO;      catch, g.LIMO     = [];  end % No default values
try g.toplot;    catch, g.toplot   = [];  end % No default values
try g.mask;      catch, g.mask     = [];  end % interactive figure
try g.title;     catch, g.title    = '';  end % interactive figure
try g.params;    catch, g.params   = [];  end % interactive figure
try g.fig;       catch, g.fig      = [];  end % Existing figure

LIMO    = g.LIMO;
toplot  = g.toplot;
mask    = g.mask;

%% get some informations for the plots

if sum(toplot(:)) == 0
    error('the image to plot is empty')
end

% what do we plot?  the data (toplot) masked (tpically of significance)
res = '';
scale           = toplot.*single(mask>0);  
scale(scale==0) = NaN;   
tmpplot = toplot;
if any(mask(:) == 0)
     tmpplot(1) = NaN;
end
cc              = limo_color_images(tmpplot); % get a color map commensurate to that

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
freqvect = [];
timevect = [];

for c=1:round(n_cluster)
    tmp                               = toplot.*(mask==c);
    tmp(tmp==Inf)                     = NaN;
    tmp(tmp==-Inf)                    = NaN;
    sigframes                         = sum(tmp,1);
    cluster_start(c)                  = find(sigframes,1,'first');
    cluster_end(c)                    = find(sigframes,1,'last');
    [V,type]                          = max([abs(min(tmp(:))) max(tmp(:))]);
    if type == 2
        cluster_maxv(c)               = V(1);
    else
        V = -V;
        cluster_maxv(c)               = V(1);
    end
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
        if exist(LIMO.dir,'dir')
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
        end
    end
    
    ratio =  abs(timevect(end)-timevect(1)) / length(timevect); % this the diff in 'size' between consecutive frames
    if LIMO.data.start < 0
        frame_zeros = find(timevect == 0);
        if isempty(frame_zeros)
            frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
        end
    else
        frame_zeros = -round(min(timevect)/ ratio);
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
    
    ratio =  abs(freqvect(end)-freqvect(1)) / length(freqvect);
    if LIMO.data.start < 0
        frame_zeros = find(freqvect == 0);
        if isempty(frame_zeros)
            frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
        end
    else
        frame_zeros = -round(min(freqvect)/ ratio);
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
        LIMO.data.tf_times = timevect;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
    end
    
    ratio =  abs(timevect(end)-timevect(1)) / length(timevect); % this the diff in 'size' between consecutive frames
    if LIMO.data.start < 0
        frame_zeros = find(timevect == 0);
        if isempty(frame_zeros)
            frame_zeros = round(abs(LIMO.data.start) / ratio)+1;
        end
    else
        frame_zeros = -round(min(timevect)/ ratio);
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
    error('LIMO.Analysis unspecified')
end

if isempty(g.title)
    if isfield(LIMO.design,'name')
        g.title = LIMO.design.name;
    else
        g.title = ' ';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the main figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(g.fig)
    fig = figure('Color','w','InvertHardCopy','off', 'paperpositionmode', 'auto');
    pos = get(fig, 'position');
    set(gcf, 'position', [pos(1:2) pos(3)*1.3 pos(4)]);

    % figure parameters
    udat.colorlim = [];
    udat.timerange = [];
    udat.y        = e;
    if ~isempty(timevect) udat.x = timevect(f); else udat.x = freqvect(f); end    
else
    fig = g.fig;
    udat = get(fig, 'userdata');
    clf(fig)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% update with mouse clicks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt   = {'maplimits','absmax','verbose','off','colormap', limo_color_images(toplot)};
udat.regressor = g.params.regressor;
udat.plot3type = g.params.plot3type;
udat.opt       = opt;
udat.cc        = cc;
udat.LIMO      = LIMO;
udat.title     = g.title;
udat.scale     = scale;
udat.params    = g.params; % parameters to limo_display_results
udat.ratio     = ratio;
udat.toplot    = toplot;
udat.timevect  = timevect;
udat.freqvect  = freqvect;
udat.frame_zeros = frame_zeros;

set(fig, 'userdata', udat);
limo_display_image_callback(fig, []);
res = true;

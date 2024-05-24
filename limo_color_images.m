function out = limo_color_images(varargin)

%% limo color map
%
% FORMAT out = limo_color_images(in,'colormapname')
%
% INPUT in - the 2D matrix to data (if in contains NaNs, as for instance a 
%            thresholded statstical map, the last value our out is gray)
%          - a single value
%          - colormapname is optional, choose between 'BWR' (blue-white-red, the default)
%                                                     'BGY' (blue-gray-yellow)
%                                                     'Hot' (black-red-yellow/white)
%
% OUTPUT if nargout==0, set current figure colormap to 'out' 
%        out is an RGB colormap  
%        if in is a singleton and <= 11, colours are taken from ColorBrewer.org
%           http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
%        if in is a singleton and > 12 colurs are sampled sample at equidistance from the xrain colour scale        
%
% References:
% see Pernet & Madan 2020 Eur J Neurosci. 51, 695–705.
% https://onlinelibrary.wiley.com/doi/abs/10.1111/ejn.14430
% ------------------------------
%  Copyright (C) LIMO Team 2024

%% defaults
color_path = [fileparts(which('limo_eeg')) filesep ...
    'external' filesep 'brain_colours' filesep];
colormapname = 'BWR';

%% check inputs
if nargin == 0
    help limo_color_images
    return
elseif nargin > 1
    colormapname = 'varargin{2}';
end

%% get the colors
in = varargin{1};
if all(size(in)==1)
    if in <= 11
        
        % colorbrewer = [166  206 227
        colorbrewer = [31 120 180
            178 223 138
            51  160 44
            251 154 153
            227 26  28
            253 191 111
            255 127 0
            202 178 214
            106 61  154
            255 255 153
            177 89  40];
        colorbrewer = colorbrewer/255;
        out = colorbrewer(1:floor(11/in):11,:);        
    else
        xrain = load([color_path 'x_rain.mat']); 
        xrain = xrain.x_rain;      
        out = xrain(1:floor(256/in):256,:);
        out = out(1:in,:);
    end

else   
    if strcmpi(colormapname,'BWR')
        out = load([color_path 'diverging_bwr.mat']); 
        out = out.diverging_bwr;
    elseif strcmpi(colormapname,'BGY')
        out = load([color_path 'diverging_bgy.mat']); 
        out = out.diverging_bgy;
    elseif strcmpi(colormapname,'hot')
        out = load([color_path 'hot.mat']); 
        out = out.hot;
    else 
        error('unsupported color map')
    end

    % half of the map if signed data in
    if min(in(:)) >= 0
        out = out(128:end,:);
    elseif max(in(:)) <= 0
        out = out(1:128,:);
    end
   
    % set NaNs to gray 
    if sum(isnan(in(:))) ~= 0
        out(1,:)=[.9 .9 .9]; 
    end
end

if nargout == 0
    colormap(out);
end


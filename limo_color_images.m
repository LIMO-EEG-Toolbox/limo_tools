function cc = limo_color_images(map)

%% limo color map
% 
% FORMAT cm = limo_color_images(map)
%
% INPUT map if the data to show (typically stat map masked by NaN)
%
% OUTPUT cc is your colormap
%
% References:
% see Pernet & Madan 2020 Eur J Neurosci. 51, 695–705.
% https://onlinelibrary.wiley.com/doi/abs/10.1111/ejn.14430
% ------------------------------
%  Copyright (C) LIMO Team 2019

% defaults
color_path = [fileparts(which('limo_eeg')) filesep 'external' filesep 'color_maps' filesep];

if min(map(:)) >= 0
    cc = load([color_path 'NIH_fire.mat']); cc = cc.lutmap2;
    cc = cc(floor(length(cc)/2):end,:);
elseif max(map(:)) <= 0
    cc = load([color_path 'NIH_cool.mat']); cc = cc.lutmap2;
else
    cc = load([color_path 'diverging_bwr.mat']); cc = cc.dmap;
    % cc = flipud(cc(1:ceil(length(cc)/2),:));
end

if sum(isnan(map(:))) ~= 0
    cc(1,:)=[.9 .9 .9]; % set NaNs to gray
end
colormap(cc);


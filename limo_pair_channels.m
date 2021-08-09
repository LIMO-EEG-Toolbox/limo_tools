function channels = limo_pair_channels(varargin)

% Routine to pair left vs right channels using the channel location information
%
% FORMAT channels = limo_pair_channels(chanlocs,'pairs',labels,'figure','on');
%
% INPUT chanlocs is the channel location information as per EEGLAB
%       'figure' can be 'on' or 'off' allowing to vizualize pairing
%       'pairs' indicate the pairs of electrode to match using a n*2 cell array of labels
%               !!! it is expected the 1 column to be left channels and 2 column right channels
%
% OUTPUT channels is the n*2 matrix of paired channels
%
% Example: LIMO = load('LIMO.mat');
%          labels = arrayfun(@(x) x.labels, LIMO.LIMO.data.chanlocs, 'UniformOutput', false);
%          table{1} = labels(1:27)'; table{2} = labels([34 35 36 39:46 49:64])';
%          channels = limo_pair_channels(LIMO.LIMO.data.chanlocs,'pairs',table,'figure','on');
%
% Cyril Pernet v1 26-03-2021
% ------------------------------
%  Copyright (C) LIMO Team 2021

chanlocs   = varargin{1};
if ischar(chanlocs)
    tmp = load(chanlocs);
    if isfield(tmp,'LIMO')
        chanlocs = tmp.LIMO.data.chanlocs;
    else
        chanlocs = tmp.(cell2mat(fieldnames(tmp)));
    end
end

fig_option = 'off';
if nargin >1
    for n=2:nargin
        if strcmpi(varargin{n},'figure')
            fig_option = varargin{n+1};
        elseif strcmpi(varargin{n},'pairs')
            inputlabels = varargin{n+1};
        end
    end
end

% create the matrix of channels
if exist('inputlabels','var')
    if iscell(inputlabels)
        if size(inputlabels) ~= [1 2]
            error('the labels cell array must be of dimension [1 2]')
        end
        channels = NaN(length(inputlabels{1}),2);
        labels   = arrayfun(@(x) x.labels, chanlocs, 'UniformOutput', false);
        for n=1:length(inputlabels{1})
            channels(n,1)  = find(cellfun(@(x) strcmpi(x,inputlabels{1}(n)),labels));
            channels(n,2)  = find(cellfun(@(x) strcmpi(x,inputlabels{2}(n)),labels));
        end
    else
        error('a n*2 cell array of labels is exected as input')
    end
    
else
    % check LIMO.data.chanlocs ; get theta and radius to be certain of position
    radius = cell2mat(arrayfun(@(x) x.radius, chanlocs, 'UniformOutput', false));
    theta  = round(cell2mat(arrayfun(@(x) x.theta, chanlocs, 'UniformOutput', false)));
    labels = arrayfun(@(x) x.labels, chanlocs, 'UniformOutput', false);
    
    % remove midline channels
    remove = [find(theta==0) find(theta ==180)];
    radius(remove) = [];
    theta(remove)  = [];
    labels(remove) = [];
    
    % pickup pairs
    Rtheta     = round(theta(theta>0));
    Rradius    = radius(find(theta>0));
    Ltheta     = round(theta(theta<0));
    Lradius    = radius(find(theta<0));
    N          = size(Rtheta,2);
    channels = NaN(N,2);
    for n=1:N
        v = Rtheta(n);
        index = find(Rtheta == v); % could be n and some other value(s)
        if length(index) > 1 % few channels at same angle
            position             = find(index == n);
            index2               = index(position);
            R_value              = find(round(theta) == v);
            channels(index2,2)   = R_value(position);
            R_value              = Rradius(index);
            R_value              = R_value(position);
            [~,position]         = min(rem(Lradius(index),R_value)); % get electrode at 'same' radius in left side
            L_value              = find(round(theta) == -v);
            channels(index2,1)   = L_value(position);
        else
            channels(index,2)  = find(theta == v); % match Rtheta to original theta (- central channels)
            if ~isempty(find(theta == -v))
                channels(index,1)  = find(theta == -v); % same angle on the other side
            else
                newindex = find(Ltheta == Ltheta(n)); % check the same indexing on the other side
                channels(index,1) = newindex;
            end
        end
    end
end

% figure
if strcmpi(fig_option,'on')
    figure('Name','Channel locations')
    subplot(1,3,[1 2]); topoplot([], chanlocs,'style','blank','electrodes','labelpoint','chaninfo',chanlocs);
    C = zeros(length(chanlocs));
    for c=1:length(channels)
        C(channels(c,1),channels(c,2)) = 1;
        C(channels(c,2),channels(c,1)) = 1;
    end
    subplot(1,3,3); imagesc(C); colormap gray; axis square
    xlabel('channels'); ylabel('channels'); title('channel pairs')
end


function [expected_chanlocs, channeighbstructmat] = limo_edit_expected_chanlocs(gp_level_file)

% Routine to edit the neighbouring matrix of (expected) chanlocs.
% While this is typically the group level matrix, this works as well for
% single subjects (i.e. loading a LIMO file).
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% if no input - ask for it
if nargin == 0
    [gp_level_file,filepath,sts]=uigetfile('*.mat','select gp level channel file or LIMO.mat');
    if sts ==0
        expected_chanlocs   = [];
        channeighbstructmat = [];
        return
    end
end

%% load file
tmp = load([filepath gp_level_file]);
if isfield(tmp,'LIMO')
    if isfield(tmp.LIMO.data,'channeighbstructmat')
        channeighbstructmat = tmp.LIMO.data.channeighbstructmat;
        if isfield(tmp.LIMO.data,'chanlocs')
            expected_chanlocs   = tmp.LIMO.data.chanlocs;
        elseif isfield(tmp.LIMO.data,'expected_chanlocs')
            expected_chanlocs   = tmp.LIMO.data.expected_chanlocs;
        end
    else
        errordg2('This LIMO fle does not have any neighbourging matrix, compute first then edit');
        return
    end
        
else % a group level file
    if ~isfield(tmp,'expected_chanlocs')
        errordlg2('this file does not have any channel location?')
        return
    elseif ~isfield(tmp,'channeighbstructmat')
        errordlg2('this file does not have any neighbouring matrix?')
        return
    else
        expected_chanlocs   = tmp.expected_chanlocs;
        channeighbstructmat = tmp.channeighbstructmat;
        clear tmp
    end
end

%% show and edit

figure('Name','Channel locations')
topoplot([], expected_chanlocs,'style','blank','electrodes','labelpoint','chaninfo',expected_chanlocs);

figure('Name','Neighbouring matrix')
imagesc(channeighbstructmat); colormap(gray);
for i=length(expected_chanlocs):-1:1
    if isfield(expected_chanlocs,'urchan')
        label{i}= expected_chanlocs(i).urchan;
    else
        label{i}= expected_chanlocs(i).labels;
    end
end
set(gca,'YTick',1:3:length(expected_chanlocs),'YTickLabel', label(1:3:length(expected_chanlocs)))
set(gca,'XTick',2:3:length(expected_chanlocs),'XTickLabel', label(2:3:length(expected_chanlocs)))
axis([1 length(expected_chanlocs) 1 length(expected_chanlocs)]); axis square
title(sprintf('Connectivity matrix between channels \n (click outside the matrix or right click when done)'),'FontSize',14)
cmap = gray; cmap(1,:) = [0.25 0.25 0.25]; colormap(cmap)

% interactive editing
positive = 1;
while positive == 1
    [x, y, button]=ginput(1);
    if any([x y]< 0) || any([x y]> length(channeighbstructmat)) || button ~= 1
        positive = 0;
    else
        
        if channeighbstructmat(round(x),round(y)) == 0
            channeighbstructmat(round(x),round(y)) = 1;
            channeighbstructmat(round(y),round(x)) = 1;
            imagesc(channeighbstructmat); v = 'on';
        else
            channeighbstructmat(round(x),round(y)) = 0;
            channeighbstructmat(round(y),round(x)) = 0;
            imagesc(channeighbstructmat);  v = 'off';
        end
        colormap(cmap);
        set(gca,'YTick',[1:3:length(expected_chanlocs)],'YTickLabel', label(1:3:length(expected_chanlocs)))
        set(gca,'XTick',[2:3:length(expected_chanlocs)],'XTickLabel', label(2:3:length(expected_chanlocs)))
        axis([1 length(expected_chanlocs) 1 length(expected_chanlocs)]); axis square
        title(sprintf('Connectivity matrix between channels \nconnection %g %g %s',round(x),round(y),v),'FontSize',14)
    end
end

%% check where to save
if nargout == 0
    if isfield(tmp,'LIMO')
        LIMO = tmp.LIMO;
        LIMO.data.channeighbstructmat = channeighbstructmat;
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO'); 
        disp('LIMO file neighbouring matrix saved')
    else
        uisave('expected_chanlocs','channeighbstructmat') % save all in one file
    end
end

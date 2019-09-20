function limo_edit_neighbours(filename,expected_chanlocs, channeighbstructmat,neighbours)

% interactive tool to edit neighbouring matrix
%
% FORMAT limo_edit_neighbours(filename,expected_chanlocs, channeighbstructmat,neighbours)
%
% INPUTS filename is the name of the file saved by limo that contains the matrices
%                 expected_chanlocs, channeighbstructmat,neighbours
%        expected_chanlocs, channeighbstructmat, neighbours are outputs
%                 from limo_expected_chanlocs.m
%        if only filename provided, it loads the file expecting to find
%        expected_chanlocs, channeighbstructmat,neighbours
%
% OUTPUT edited_filename saved on the drive
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% deal with inputs
if nargin == 1
    load(filename);
    if ~exist('expected_chanlocs','var') || ...
            ~exist('channeighbstructmat','var') || ...
            ~exist('neighbours','var')
        error('from loading %s some matrices are missing',filename)
    end
elseif nargin ~= 4
    error('1 or 4 arguments in expected')
end

%% plot data
figure('Name','Channels')
topoplot([], expected_chanlocs,'style','blank','electrodes','labelpoint','chaninfo',expected_chanlocs);

% show connectivity matrix
figure('Name','Neighbouring matrix')
imagesc(channeighbstructmat); colormap(gray);
for i=1:length(expected_chanlocs)
    if isfield(expected_chanlocs,'urchan')
        label{i}= expected_chanlocs(i).urchan;
    elseif isfield(expected_chanlocs,'labels')
        label{i}= expected_chanlocs(i).labels;
    end
end
set(gca,'YTick',[1:3:length(expected_chanlocs)],'YTickLabel', label(1:3:length(expected_chanlocs)))
set(gca,'XTick',[2:3:length(expected_chanlocs)],'XTickLabel', label(2:3:length(expected_chanlocs)))
axis([1 length(expected_chanlocs) 1 length(expected_chanlocs)]); axis square
title(sprintf('Connectivity matrix between channels \n (click outside the matrix or right click when done)'),'FontSize',14)
cmap = gray; cmap(1,:) = [0.25 0.25 0.25]; colormap(cmap)

%% interactive editing

positive = 1;
newneighbours          = neighbours;
newchanneighbstructmat = channeighbstructmat;
while positive == 1
    [x, y, button]=ginput(1);
    if any([x y]< 0) || any([x y]> length(newchanneighbstructmat)) || button ~= 1
        positive = 0;
    else
        
        if newchanneighbstructmat(round(x),round(y)) == 0
            fprintf('Pairing channels %s %s\n',neighbours{round(x)}.label,neighbours{round(y)}.label)
            newchanneighbstructmat(round(x),round(y))  = 1;
            newchanneighbstructmat(round(y),round(x))  = 1;
            imagesc(newchanneighbstructmat); v = 'on';
            newneighbours{round(x)}.neighblabel{end+1} = neighbours{round(y)}.label;
            newneighbours{round(x)}.neighblabel        = sort(newneighbours{round(x)}.neighblabel);
            newneighbours{round(y)}.neighblabel{end+1} = neighbours{round(x)}.label;
            newneighbours{round(y)}.neighblabel        = sort(newneighbours{round(y)}.neighblabel);
        else
            fprintf('Unpairing channels %s %s\n',neighbours{round(x)}.label,neighbours{round(y)}.label)
            newchanneighbstructmat(round(x),round(y)) = 0;
            newchanneighbstructmat(round(y),round(x)) = 0;
            imagesc(newchanneighbstructmat);  v = 'off';
            rmfd = find(cellfun(@isempty,strfind(newneighbours{round(x)}.neighblabel,neighbours{round(y)}.label)) == 0);
            newneighbours{round(x)}.neighblabel(rmfd) = [];
            rmfd = find(cellfun(@isempty,strfind(newneighbours{round(y)}.neighblabel,neighbours{round(x)}.label)) == 0);
            newneighbours{round(y)}.neighblabel(rmfd) = [];            
        end
        
        disp('neighbouring:'); disp(newneighbours{round(x)}.neighblabel')
        disp(newneighbours{round(y)}.neighblabel')
       
        colormap(cmap);
        set(gca,'YTick',[1:3:length(expected_chanlocs)],'YTickLabel', label(1:3:length(expected_chanlocs)))
        set(gca,'XTick',[2:3:length(expected_chanlocs)],'XTickLabel', label(2:3:length(expected_chanlocs)))
        axis([1 length(expected_chanlocs) 1 length(expected_chanlocs)]); axis square
        title(sprintf('Connectivity matrix between channels \nconnection %g %g %s',round(x),round(y),v),'FontSize',14)
    end
end

%% save
eqtest = single(eq(newchanneighbstructmat,channeighbstructmat));
if sum(eqtest(:)) ~= numel(channeighbstructmat)
    neighbours              = newneighbours;
    channeighbstructmat     = newchanneighbstructmat;
    [filepath,filename,ext] = fileparts(filename);
    save([filepath filesep 'edited_' filename ext],'expected_chanlocs','channeighbstructmat','neighbours') % save all in one file
    fprintf('neighbouring matrices saved\n');
end
close all

function out = limo_colormap

% create a hot colormap going from gray to yellow
% if no output - applies the map
% if output returns the map to be used with colormap
% cyril pernet 01-05-2014

% gray   = [0.5 0.5 0.5]
% red    = [1 0 0]
% yellow = [1 1 0]

out = [];
A = [[0.7:0.3/31:1]'; ones(32,1)]; % gray(1) to red(1) to yellow(1) 
B = [[0.7:-0.7/31:0]'; [0:1/31:1]']; % gray(2) to red(2) to yellow(2) 
C = [[0.7:-0.7/31:0]'; zeros(32,1)]; % gray(3) to red(3) to yellow(3) 
if nargin == 1
    out = [A B C];
else
    colormap([A B C]);
end
function limo_system_info

% simple routine to figure out what we can do to speed things up 
% either by allowing large arrays and or multithreading using psom
%
% Cyril Pernet May 2014
% -----------------------------
%  Copyright (C) LIMO Team 2014

root = fileparts(which('limo_eeg'));
new = ~exist([root filesep 'system_info.mat'],'file');

if new % fresh install
    [~,system_info.max_elements] = computer;
    system_info.max_cores = feature('numCores');
    if ispc
        how_much_ram = memory;
        number_of_elements = how_much_ram.MaxPossibleArrayBytes / 8; % store data as double
        memory_available = how_much_ram.MemAvailableAllArrays;
    elseif isunix
        [r,w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        memsize = stats(1)/1e6;
        freemem = (stats(3)+stats(end))/1e6;
    elseif ismac
        % [r,w] = unix('top');
    end
else
    load([root filesep 'system_info.mat'])
    % quick check it still the same system
    [r,w] = unix('free | grep Mem');
end

% stuff for psom
% test psom config
opt_pipe.mode = 'background';
psom_config('',opt_pipe);

% use cell2csv to make a psom config file

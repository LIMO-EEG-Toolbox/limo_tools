function limo_check_ppool

% routine to make sure the parallel pool is working
% ------------------------------
%  Copyright (C) LIMO Team 2020

p = gcp('nocreate');
if isempty(p) % i.e. the parallel toolbox is not on
    
    % check how many cores to use
    % ---------------------------
    N = getenv('NUMBER_OF_PROCESSORS'); % logical nu,ber of cores (i.e. count hyperthreading)
    % N = feature('numcores');          % physical number of cores
    if ischar(N)
        N = str2double(N);
    end
    
    % check and set the local profile NumWorkers
    % ----------------------------------
    c            = parcluster;
    c.NumWorkers = N;
    saveProfile(c);
    
    % go
    % --
    parpool(N-1); % logical
end

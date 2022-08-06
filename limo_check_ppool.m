function limo_check_ppool

% routine to make sure the parallel pool is working
% you can set manually N below, so that the pool is fix
% -----------------------------------------------------
%  Copyright (C) LIMO Team 2022

N      = [];
addons = ver;

if any(strcmpi('Parallel Computing Toolbox',arrayfun(@(x) x.Name, addons, "UniformOutput",false)))
    
    if isempty(N)
        p = gcp('nocreate');
        if isempty(p) % i.e. the parallel toolbox is not already on
            
            % check how many cores to use
            % ---------------------------
            N = getenv('NUMBER_OF_PROCESSORS'); % logical number of cores (i.e. count hyperthreading)
            % N = feature('numcores');          % physical number of cores
            if ischar(N)
                N = str2double(N);
            end
            
            % check and set the local profile NumWorkers
            % ----------------------------------
            c            = parcluster;
            c.NumWorkers = N;
            % saveProfile(c);
            
            % go
            % --
            parpool(N-1);
        end
    else
        parpool(N);
    end
    
else
    disp('Parallel toolbox not found - nothing to worry about (except slower computation in some cases)');
end


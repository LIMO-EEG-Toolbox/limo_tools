function limo_check_ppool

% routine to make sure the parallel pool is working
% you can set manually N below, so that the pool is fix
% -----------------------------------------------------
%  Copyright (C) LIMO Team 2022

N      = [];
addons = ver;

if any(strcmpi('Parallel Computing Toolbox',arrayfun(@(x) x.Name, addons, "UniformOutput",false)))
    % the toolbox exists
    if isempty(N)
        p = gcp('nocreate');
        if isempty(p) % i.e. the parallel toolbox is not turned on
            
            % check how many cores to use
            % ---------------------------
            N = getenv('NUMBER_OF_PROCESSORS'); % logical number of cores (i.e. count hyperthreading)
            if isempty(N)
                N = feature('numcores');        % physical number of cores
            elseif ischar(N)
                N = str2double(N);
            end
            
            % check and set the local profile NumWorkers
            % ----------------------------------
            c            = parcluster('local');
            c.NumWorkers = N;
            % saveProfile(c);
            
            % go
            % --
            if N>1
              c.parpool(N-1);
            else
              disp('Parallel toolbox not started - nothing to worry about (except slower computation in some cases)');
              ps = parallel.Settings;
              ps.Pool.AutoCreate = false; % if parfor loop, do not try to start ppool
            end
        end
    else % this is already running
        try
            parpool(N);
        catch errpool
            warning(errpool.identifier,'could not create a parallel pool: %s',errpool.message)
            try poolobj = gcp('nocreate'); delete(poolobj); limo_check_ppool; end % close pool and restart
        end
    end    
else
    disp('Parallel toolbox not found - nothing to worry about (except slower computation in some cases)');
end


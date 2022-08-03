function limo_check_ppool

% routine to make sure the parallel pool is working
% ------------------------------
%  Copyright (C) LIMO Team 2020

if ~exist('gcp')
    disp('Parallel toolbox not found - nothing to worry about (except slower computation in some cases)');
else
    p = gcp('nocreate');
    if isempty(p) % i.e. it's not on
        parpool(getenv('NUMBER_OF_PROCESSORS')-1); % logical
        % parpool(feature('numcores')-1); only physical
    end
end


function limo_check_ppool

% routine to make sure the parallel pool is working
% ------------------------------
%  Copyright (C) LIMO Team 2020

p = gcp('nocreate');
if isempty(p) % i.e. it's not on
    parpool(feature('numcores')-1);
end

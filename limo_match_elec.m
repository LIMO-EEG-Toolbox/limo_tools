function out = limo_match_elec(c_elec,e_elec,a_beg,a_end,data)

% LIMO_MATCH_ELEC - called to match currently available electrodes from 
% one subject with the set of expected electrodes across all subjects.
%
% syntax: out = limo_match_elec(c_elec,e_elec,a_beg,a_end,data)
%
% INPUTS:
%   c_elec = chanlocs structure of current electrodes available for that
%            subject
%   e_elec = chanlocs structure of expected electrodes for all subjects
%   a_beg  = first frame for the analysis
%   a_end  = last frame for the analysis
%   data   = matrix electrodes x frames x beta/con coefficients
%
% OUTPUT:
%   out   = matrix electrodes x frames x beta coefficients in which the
%           number of electrodes is the number of expected electrodes. When some 
%           expected electrodes are not available for a subject, the data are
%           replaced by NaNs
%
% See also LIMO_RANDOM_SELECT.
%
% Guillaume Rousselet, Cyril Pernet v1 24-September-2009
% -----------------------------
%  Copyright (C) LIMO Team 2010


out = NaN(size(e_elec,2),size([a_beg:a_end],2),size(data,3));

for current=1:length(c_elec)
    for expected=1:size(e_elec,2)
        if strcmp(c_elec(current).labels,e_elec(expected).labels)
            out(expected,:,:) = data(current,a_beg:a_end,:);
        end
    end
end



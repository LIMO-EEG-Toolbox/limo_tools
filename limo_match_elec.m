function out = limo_match_elec(c_elec,e_elec,a_beg,a_end,data)

% LIMO_MATCH_ELEC - called to match currently available electrodes from 
% one subject with the set of expected electrodes (usually to be used
% across all subjects at the second level).
%
% syntax: out = limo_match_elec(c_elec,e_elec,a_beg,a_end,data)
%
% INPUTS:
%   c_elec = chanlocs structure of current electrodes available for that
%            subject
%   e_elec = chanlocs structure of expected electrodes for all subjects
%   a_beg  = first frame for the analysis (for Time-Frquency this is a
%            vector with the 1st freq and 1st time)
%   a_end  = last frame for the analysis (for Time-Frquency this is a
%            vector with the 1st freq and 1st time)
%   data   = matrix electrodes * frames * data
%            matrix electrodes * freq * time * data
%
% OUTPUT:
%   out   = matrix electrodes * frames * data or electrodes * freq * time * data
%           in which the number of electrodes is the number of expected electrodes.
%           When some expected electrodes are not available for a subject, the 
%           data are replaced by NaNs
%
% See also LIMO_RANDOM_SELECT.
%
% Guillaume Rousselet, Cyril Pernet v1 24-September-2009
% Cyril Pernet May 2014, updates for time-freqency
% ----------------------------------------------------------
%  Copyright (C) LIMO Team 2014

if numel(size(data)) == 4
    out = NaN(size(e_elec,2),size([a_beg(1):a_end(1)],2),size([a_beg(2):a_end(2)],2),size(data,4));
    
    for current=1:length(c_elec)
        for expected=1:size(e_elec,2)
            if strcmp(c_elec(current).labels,e_elec(expected).labels)
                out(expected,:,:,:) = data(current,a_beg(1):a_end(1),a_beg(2):a_end(2),:);
            end
        end
    end
else
    out = NaN(size(e_elec,2),size([a_beg:a_end],2),size(data,3));
    
    for current=1:length(c_elec)
        for expected=1:size(e_elec,2)
            if strcmp(c_elec(current).labels,e_elec(expected).labels)
                out(expected,:,:) = data(current,a_beg:a_end,:);
            end
        end
    end
end




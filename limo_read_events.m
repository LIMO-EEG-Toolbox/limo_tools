function cat = limo_read_events(varargin)

% simple routine that reads EEG.event.type to create a categorical variable
% that can be used in LIMO EEG
%
% FORMAT limo_read_events
%        limo_read_events(data.set)
%        limo_read_events(data.set,markers)
%
% INPUT <empty> will read the EEG variable already in the workplace as well
%               as all the events associared and create a cat.mat file in
%               the dirrectory associated to the .set 
%       data.set is the full name (with path) of a .set file to read data 
%       from -- EXPECT EPOCHED DATA
%       markers is a cell array with the marker name in each cell, for
%       instance the .set can have 10 different event type but only 5 need
%       to be analyzed - the markers array indicates which ones should be
%       used
%
% OUTPUT <empty> will create and save a cat.mat file that can be used at
%        import or directly in limo_design_matrix(_tf)
%        cat is the variable to be passed in limo_design matrix(_tf)
%
% Cyril Pernet 13-10-2014
% ------------------------------------------
% Copyright (C) LIMO Team 2014

% get input file to get the EEG structure
if nargin == 0
    storein = EEG.filepath;
else
    EEG=pop_loadset(varargin{1});
    [storein,file,ext]=fileparts(varargin{1});
end

markers = 'all';
if nargin == 2
    clear markers
   markers = varargin{2}; 
end

% now get the different events
for n=1:size(EEG.event,2)
    type{n} = EEG.event(n).type;
end
categories = unique(type);
nevents = length(categories);
fprintf('%g event types were detected \n',nevents)

% check that markers and event match
cat = NaN(size(EEG.event,2),1);
if iscell(markers)
    for n=1:size(EEG.event,2)
        for nn = 1:size(markers,2)
            if strcmp(markers{nn},type{n})
                cat(n) = nn;
            end
        end
    end
else
    for n=1:size(EEG.event,2)
        for nn = 1:nevents
            if strcmp(categories{nn},type{n})
                cat(n) = nn;
            end
        end
    end
end
        



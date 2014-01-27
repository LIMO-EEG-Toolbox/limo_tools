function [data] = limo_eeglab2fieldtrip(EEG, fieldbox, transform)

% LIMO_EEGLAB2FIELDTRIP converts an EEGLAB data structure into a FieldTrip data
% structure, which subsequently can be used for dipole fitting and other
% FieldTrip analysis methods.
%
% Use as
%   [data] = eeglab2fieldtrip( EEG, fieldbox )
% or
%   [data] = eeglab2fieldtrip( EEG, fieldbox, transform )
%
% where the inputs are
%   EEG       - [struct] EEGLAB structure
%   fieldbox  - ['preprocessing'|'timelockanalysis'|'componentanalysis'|...
%                'chanloc', 'chanloc_withfid']
%   transform - ['none'|'dipfit'] transform channel locations for DIPFIT
%               using the transformation matrix in the field 'coord_transform'
%               of the dipfit substructure of the EEG structure.
%
% The output is a structure that is compatible with the FieldTrip functions
% PREPROCESSING, TIMELOCKANALYSIS, COMPONENTANALYSIS, or only containing an
% electrode definition.

% Copyright (C) 2004-2006, Robert Oostenveld, FCDC
% Copyright (C) 2004-2006, Arnaud Delorme, SCCN, INC, UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: eeglab2fieldtrip.m,v $
% Revision 1.14  2007/07/07 12:10:00  roboos
% detect matlab version 7.0.1 and 7.0.4 and treat as 7.0
%
% Revision 1.13  2006/12/07 08:37:32  roboos
% compute component activations from only those channels used in the ICA
%
% Revision 1.12  2006/11/23 11:18:42  roboos
% moved the naming of channel labels into the switch (to avoid confusion with componentanalys fieldbox), fixed non-string channel labels (now using sprintf)
%
% Revision 1.11  2006/11/23 10:40:38  roboos
% Arno fixed an important bug how single trials component activation was
% computed (sphering was forgotten). Restructured help
%
% Revision 1.10  2006/11/23 10:15:00  roboos
% Use numeric channel labels if no chanlocs are present. Removed fieldtripchan2eeglab subfunction, since it was not used any more.
%
% Revision 1.9  2006/10/05 09:45:11  roboos
% included matlabversion as subfunction
%
% Revision 1.8  2006/10/05 09:39:56  roboos
% fixes an important bug when you do not use all the channels to run ica (arno),
% included the fixprecision function as subfunction
%
% Revision 1.7  2006/05/29 08:41:22  roboos
% ensure that data is double precision for matlab versions<7
% ensure that labels are returned as a column (thanks to Theresa)
%
% Revision 1.6  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.5  2006/04/13 09:09:08  roboos
% get the precomputed activations if possible, or alternatively compute them
%
% Revision 1.4  2006/03/20 08:20:41  roboos
% added ; to the end of a line (thanks to arno)
%
% Revision 1.3  2006/02/27 10:02:08  roboos
% added optional transform for channel locations
% added fiducials to channel locations
%
% Revision 1.1  2006/01/20 23:43:18  arno
% Initial revision
%
% renamed for integration in LIMO toolbox: GAR, University of Glasgow, June
% 2010

if nargin < 2
  help eeglab2fieldtrip
  return;
end;

% start with an empty data object
data = [];

if strcmpi(fieldbox, 'chanloc_withfid')
  % also insert the non-data channels (i.e. fiducials) in channel structure
  if isfield(EEG.chaninfo, 'nodatchans')
    chanlen = length(EEG.chanlocs);
    fields = fieldnames( EEG.chaninfo.nodatchans );
    for index = 1:length(EEG.chaninfo.nodatchans)
      ind = chanlen+index;
      for f = 1:length( fields )
        EEG.chanlocs = setfield(EEG.chanlocs, { ind }, fields{f}, ...
          getfield( EEG.chaninfo.nodatchans, { index },  fields{f}));
      end
    end
  end
end

% get the electrode positions from the EEG structure: in principle, the number of
% channels can be more or less than the number of channel locations, i.e. not
% every channel has a position, or the potential was not measured on every
% position. This is not supported by EEGLAB, but it is supported by FIELDTRIP.

Nchan           = length(EEG.chanlocs);
data.elec.pnt   = zeros(Nchan,3);
data.elec.label = cell(Nchan,1);
for ind = 1:length( EEG.chanlocs )
  data.elec.label{ind} = EEG.chanlocs(ind).labels;
  if ~isempty(EEG.chanlocs(ind).X)
    % this channel has a position
    data.elec.pnt(ind,1) = EEG.chanlocs(ind).X;
    data.elec.pnt(ind,2) = EEG.chanlocs(ind).Y;
    data.elec.pnt(ind,3) = EEG.chanlocs(ind).Z;
  else
    % this channel does not have a position
    data.elec.pnt(ind,:) = [0 0 0];
  end;
end;

if nargin > 2
  if strcmpi(transform, 'dipfit')
    if ~isempty(EEG.dipfit.coord_transform)
      disp('Transforming electrode coordinates to match head model');
      transfmat = traditional(EEG.dipfit.coord_transform);
      data.elec.pnt = transfmat * [ data.elec.pnt ones(size(data.elec.pnt,1),1) ]';
      data.elec.pnt = data.elec.pnt(1:3,:)';
    else
      disp('Warning: no transformation of electrode coordinates to match head model');
    end;
  end;
end;

switch fieldbox
  case { 'chanloc' 'chanloc_withfid' }
    % these have already been processed above

  case 'preprocessing'
    for index = 1:EEG.trials
      data.trial{index}  = fixprecision(EEG.data(:,:,index));
      data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    end;
    if ~isempty(EEG.chanlocs)
      % use the channel labels from the EEG structure
      data.label = { EEG.chanlocs(1:EEG.nbchan).labels };
    else
      % construct channel names on the fly
      for i=1:EEG.nbchan
        data.label{i} = sprintf('chan%03d', i);
      end
    end
    data.label   = data.label(:);
    data.fsample = EEG.srate;

  case 'timelockanalysis'
    data.avg  = mean(fixprecision(EEG.data), 3);
    data.var  = std(fixprecision(EEG.data), [], 3).^2;
    data.time = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    data.dimord = 'chan_time';
    if ~isempty(EEG.chanlocs)
      % use the channel labels from the EEG structure
      data.label = { EEG.chanlocs(1:EEG.nbchan).labels };
    else
      % construct channel names on the fly
      for i=1:EEG.nbchan
        data.label{i} = sprintf('chan%03d', i);
      end
    end
    data.label   = data.label(:);
    data.fsample = EEG.srate;

  case 'componentanalysis'
    try,
      for index = 1:EEG.trials
        % the trials correspond to the raw data trials, except that they contain the component activations
        if isempty(EEG.icaact)
          data.trial{index}  = fixprecision(EEG.icaweights*EEG.icasphere*EEG.data(EEG.icachansind,:,index));
        else
          data.trial{index}  = EEG.icaact(:,:,index);
        end;
        data.time{index}   = linspace(EEG.xmin, EEG.xmax, EEG.pnts); % should be checked in FIELDTRIP
      end;
    catch, end;
    for comp = 1:size(EEG.icawinv,2)
      % the labels correspond to the component activations that are stored in data.trial
      data.label{comp} = sprintf('ica%03d', comp);
    end
    % get the spatial distribution and electrode positions
    data.topolabel = { EEG.chanlocs(EEG.icachansind).labels };
    data.topo      = fixprecision(EEG.icawinv);
    % ensure that the labels are in colunms
    data.label     = data.label(:);
    data.topolabel = data.topolabel(:);

  case 'freqanalysis'
    error('freqanalysis fieldbox not implemented yet')

  otherwise
    error('unsupported fieldbox')
end

try
  % get the full name of the function
  data.cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  data.cfg.version.name = st(i);
end

% add the version details of this function call to the configuration
data.cfg.version.id   = '$Id: eeglab2fieldtrip.m,v 1.14 2007/07/07 12:10:00 roboos Exp $';
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXPRECISION ensures that the data matrix has a numeric precision that is
% supported by the Matlab version
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x] = fixprecision(x, precision);
if nargin<2
  if matlabversion<7
    % Matlab versions prior to 7.0 cannot perform computations on single
    % precision numbers
    precision = 'double';
  else
    precision = class(x);
  end
end
if ~strcmp(class(x), precision)
  % conversion is needed
  switch precision
    case 'single'
      x = single(x);
    case 'double'
      x = double(x);
    otherwise
      error('unsupported numeric precision');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLABVERSION returns the Matlab version as a number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v] = matlabversion;
s = ver('matlab');
v = s.Version;
if ischar(v) 
  % try converting to a number
  n = str2num(v);
  if isempty(n)
    switch v
    case '6.5.1'
      n = 6.5; % this is accurate enough
    case '7.0.1'
      n = 7.0; % this is accurate enough
    case '7.0.4'
      n = 7.0; % this is accurate enough
    otherwise  
      warning('cannot convert matlab version into a number');
      v = v;
    end
  end
  v = n;
end


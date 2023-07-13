function res = limo_questdlg(varargin)

% if EEGLAB present using its' GUI
% Arnaud Delorme
% -----------------------------------
%  Copyright (C) LIMO Team 2023

if exist('questdlg2','file')
    res = questdlg2(varargin{:});
else
    res = questdlg(varargin{:});
end
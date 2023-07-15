function limo_warndlg(varargin)

% if EEGLAB present using its' GUI
% Arnaud Delorme
% -----------------------------------
%  Copyright (C) LIMO Team 2023

if exist('warndlg2','file')
    warndlg2(varargin{:});
else
    warning(varargin{:});
end
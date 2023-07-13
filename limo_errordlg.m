function limo_errordlg(varargin)

% if EEGLAB present using its' GUI
% Arnaud Delorme
% -----------------------------------
%  Copyright (C) LIMO Team 2023

if exist('errordlg2','file')
    errordlg2(varargin{1});
else
    errordlg(varargin{1});
end

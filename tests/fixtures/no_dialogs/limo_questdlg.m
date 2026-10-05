function varargout = limo_questdlg(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected limo_questdlg call in a fully specified command line analysis.');
end

function varargout = listdlg(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected listdlg call in a fully specified command line analysis.');
end

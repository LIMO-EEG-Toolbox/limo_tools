function varargout = inputdlg(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected inputdlg call in a fully specified command line analysis.');
end

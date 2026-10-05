function varargout = warndlg(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected warndlg call in a fully specified command line analysis.');
end

function varargout = errordlg(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected errordlg call in a fully specified command line analysis.');
end

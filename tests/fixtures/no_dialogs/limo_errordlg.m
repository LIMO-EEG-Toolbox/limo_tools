function varargout = limo_errordlg(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected limo_errordlg call in a fully specified command line analysis.');
end

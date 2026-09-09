function varargout = limo_inputdlg(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected limo_inputdlg call in a fully specified command line analysis.');
end

function varargout = warndlg2(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected warndlg2 call in a fully specified command line analysis.');
end

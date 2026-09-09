function varargout = uigetfile(varargin)
% Fail instead of blocking unattended command line regression tests.
error('LIMO:testUnexpectedDialog', 'Unexpected uigetfile call in a fully specified command line analysis.');
end

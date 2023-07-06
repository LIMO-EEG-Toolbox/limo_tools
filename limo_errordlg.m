function limo_errordlg(varargin)

if exist('errordlg2','file')
    errordlg2(varargin{1});
else
    errordlg(varargin{1});
end
%error(varargin{1})
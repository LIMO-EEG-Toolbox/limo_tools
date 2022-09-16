function limo_errordlg(varargin)

if exist('errordlg2','file')
    errordlg2(varargin{:});
else
    errordlg(varargin{:});
end
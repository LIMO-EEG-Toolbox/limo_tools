function res = limo_inputdlg(varargin)

if exist('warndlg2','file')
    res = inputdlg2(varargin{:});
else
    res = inputdlg(varargin{:});
end
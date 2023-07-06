function limo_warndlg(varargin)

if exist('warndlg2','file')
    warndlg2(varargin{:});
else
    warndlg(varargin{:});
end
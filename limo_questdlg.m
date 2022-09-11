function res = limo_questdlg(varargin)

if exist('questdlg2','file')
    res = questdlg2(varargin{:});
else
    res = questdlg(varargin{:});
end
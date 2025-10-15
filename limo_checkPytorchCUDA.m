function info = limo_checkPytorchCUDA()
% Returns a struct with fields: isAvailable (logical), deviceCount (int),
% deviceNames (cell array of strings)
%
% make sure you matlab/python are communicating as needed
% for instance to make it work on a server using a virtual
% environment you do
% pyenv( ...
%     "Version", "/users/cyrilpernet/conda-envs/venv/bin/python", ...
%     "ExecutionMode", "OutOfProcess" ...
% );
% where the Version you want to use point to your python, and since
% this is a virtual environement, conda activate, then call matlab 
% and run limo tools 
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2025

try
    % Import torch
    pyrun("import torch");
catch ME
    info.error = ['could not import PyTorch, pip install it: ' ME.message];
end

try
    % info.version = string(pyrun("import warnings; warnings.filterwarnings('ignore'); import torch; ver = torch.__version__","ver"));
    info.version = string(pyrun("import torch; ver = str(torch.__version__)", "ver"));
    % Query availability
    isOk = pyrun("ok = torch.cuda.is_available()", "ok");
    info.isAvailable = logical(isOk);
    if info.isAvailable
        cnt = pyrun("cnt = torch.cuda.device_count()", "cnt");
        info.deviceCount = int32(cnt);
        names = cell(1, double(cnt));
        for k = 0:cnt-1
            nm = pyrun(sprintf("nm=torch.cuda.get_device_name(%d)", k), "nm");
            names{k+1} = char(nm);   % convert Python str to MATLAB char
        end
        info.deviceNames = names;
    else
        info.deviceCount = 0;
        info.deviceNames = {};
    end
catch ME
    % Something went wrong (Python/PyTorch not installed, import error, etc.)
    info.isAvailable = false;
    info.deviceCount = 0;
    info.deviceNames = {};
    info.error       = ME.message;
end


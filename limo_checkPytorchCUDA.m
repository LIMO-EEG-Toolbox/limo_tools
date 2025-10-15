function info = limo_checkPytorchCUDA()
% Returns a struct with fields: isAvailable (logical), deviceCount (int),
% deviceNames (cell array of strings)
% ----------------------------------------
%  Copyright (C) LIMO Team 2025

try
    % Import torch
    pyrun("import torch", []);
    % Query availability
    isOk = pyrun("torch.cuda.is_available()", "torch.cuda.is_available");
    info.isAvailable = logical(isOk);
    if info.isAvailable
        cnt = pyrun("torch.cuda.device_count()", "torch.cuda.device_count");
        info.deviceCount = int32(cnt);
        names = cell(1, double(cnt));
        for k = 0:cnt-1
            nm = pyrun(sprintf("torch.cuda.get_device_name(%d)", k), "torch.cuda.get_device_name");
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


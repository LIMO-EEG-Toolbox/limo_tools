function boot_table = limo_boot_table_get(boot_file, var_name, data, nboot)

% Get the bootstrap resampling table for a chunked/resumable LIMO bootstrap.
%   * if no table exists yet -> create one with limo_create_boot_table and save
%   * if a table exists with >= nboot columns -> reuse it as-is (so a resumed
%     run continues the exact same resampling -> identical results)
%   * if a table exists with < nboot columns -> APPEND fresh resampling columns
%     to reach nboot and re-save (this is the "add more bootstraps later"
%     feature; bootstrap iterations are exchangeable so appended columns are
%     valid and leave the earlier ones untouched)
%
% boot_file : full path to the .mat (with or without extension)
% var_name  : variable name stored inside (e.g. 'boot_table','boot_table1')
% data      : [channels x frames x subjects], passed to limo_create_boot_table
% nboot     : desired number of bootstrap columns
% ------------------------------
%  Copyright (C) LIMO Team 2026

if ~endsWith(boot_file,'.mat'), boot_file = [boot_file '.mat']; end

if exist(boot_file,'file')
    S = load(boot_file);
    boot_table = S.(var_name);
    have = 0;
    for c = 1:numel(boot_table)
        if ~isempty(boot_table{c}), have = size(boot_table{c},2); break; end
    end
    if have < nboot
        fprintf('extending boot table %s: %d -> %d columns\n', var_name, have, nboot);
        extra = limo_create_boot_table(data, nboot - have);
        for c = 1:numel(boot_table)
            if ~isempty(boot_table{c}) && c <= numel(extra) && ~isempty(extra{c})
                boot_table{c} = [boot_table{c}, extra{c}];
            end
        end
        tmp.(var_name) = boot_table; save(boot_file, '-struct', 'tmp');
    end
else
    boot_table = limo_create_boot_table(data, nboot);
    tmp.(var_name) = boot_table; save(boot_file, '-struct', 'tmp');
end
end

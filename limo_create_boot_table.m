function boot_table = limo_create_boot_table(data,nboot)

% this function allows building a table of data to resample
% such as almost the same resampling is applied across electrodes.
% At the 2nd level, not all subjects have the same electrodes, and
% missing data are treated as NaN. The boot_table has indexes which are
% common to all electrodes + some other indexes specific to electrodes
% where there some NaNs.
%
% FORMAT: boot_table = limo_create_boot_table(data,nboot)
%
% INPUT: data: a 3D matrix [electrode x time frames x subjects];
%        nboot: the number of bootstraps to do
%
% OUTPUT: boot_table is a cell array with one resampling matrix per
% electrode
%
% Cyril Pernet v1 24-06-2013
% ----------------------------
% Copyright (C) LIMO Team 2010

% check data for NaNs
if size(data,1) == 1
    chdata=data(1,1,:); 
else
    chdata=squeeze(data(:,1,:)); 
end
    
% create boot_table
B=1;
boot_index=zeros(size(data,3),nboot);
while B~=nboot+1
    tmp = randi(size(data,3),size(data,3),1);
    if length(unique(tmp)) >= 2 && min(sum(~isnan(chdata(:,tmp)),2)) > 2; % at least 3 different observations per boot and data collected
        boot_index(:,B) = tmp;
        B=B+1;
    end
end
clear chdata tmp

% loop per electrode, if no nan use boot_index else change it
for electrode = 1:size(data,1)
    tmp = squeeze(data(electrode,:,:)); % 2D
    Y = tmp(:,find(~isnan(tmp(1,:)))); % remove NaNs
    bad_subjects = find(isnan(tmp(1,:)));
    good_subjects = find(~isnan(tmp(1,:)));
    
    if ~isempty(bad_subjects)
        boot_index2 = zeros(size(Y,2),nboot);
        for c=1:nboot
            common = ismember(boot_index(:,c),good_subjects');
            current = boot_index(find(common),c); % keep resampling of good subjets
            % add or remove indices
            add = size(Y,2) - size(current,1);
            if add > 0
                new_boot = [current ; good_subjects(randi(size(good_subjects),add,1))'];
            else
                new_boot = current(1:size(Y,2));
            end
            % change indices values
            tmp_boot = new_boot;
            for i=1:length(bad_subjects)
                new_boot(find(tmp_boot > bad_subjects(i))) = tmp_boot(find(tmp_boot >  bad_subjects(i))) - i; % change range
            end
            % new boot-index
            boot_index2(:,c) = new_boot;
        end
        boot_table{electrode} = boot_index2;
    else
        boot_table{electrode} = boot_index;
    end
end


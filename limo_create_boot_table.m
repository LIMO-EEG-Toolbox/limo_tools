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
%        for repeated measures - call this function inputing 1 measure and
%        apply to all measures the created table
%
% OUTPUT: boot_table is a cell array with one resampling matrix per
% electrode
%
% Cyril Pernet v1 24-06-2013
% v2 simple edit to specify up front Nmin
% --------------------------------------
% Copyright (C) LIMO Team 2015

%% edit default
Nmin = 3; % this is the minimum number of different trials/subjects 
          % if too low, the variance is < 1 and thre stat values will be
          % too high see Pernet et al. 2014

%% start
% check data for NaNs
if size(data,1) == 1
    chdata=data(1,1,:); 
else
    chdata=squeeze(data(:,1,:)); 
end

if sum(sum(isnan(chdata),2)==size(data,ndims(data))) ~=0
    disp('some electrodes are empty (full of NaN) - still making the table but some cells will be empty')
end

if (sum((size(data,3) - sum(isnan(chdata),2))<=3)) ~=0
    disp('some cells have a very low count <=3 ; bootstrapping cannot work - still making the table but some cells will be empty')
end

% create boot_table
B=1;
boot_index=zeros(size(data,3),nboot);
if size(data,3)-1 <= Nmin
    error(['Not enough subjects in dataset - need at least ' num2str(Nmin+2) ' subjects']);
end;
while B~=nboot+1
    tmp = randi(size(data,3),size(data,3),1);
    if length(unique(tmp)) >= Nmin % at least Nmin different observations per boot 
        boot_index(:,B) = tmp;
        B=B+1;
    end
end
clear chdata tmp

% loop per electrode, if no nan use boot_index else change it
if size(data,1) > 1
    array = find(sum(squeeze(isnan(data(:,1,:))),2) < size(data,3)-3);
else
    array = [1];
end;
for e = 1:size(array,1)
    electrode = array(e);
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


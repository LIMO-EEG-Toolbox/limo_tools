function [expected_chanlocs] = limo_avg_expected_chanlocs(PATH_TO_DERIVATIVES,defaults)

dinfo = dir(fullfile(PATH_TO_DERIVATIVES,'sub-*'));
X = [];
Y = [];
Z = [];
for i = 1 : numel( dinfo )
    if i >= 10
        subfolder = 'sub-0%d';
    else
        subfolder = 'sub-00%d';
    end
    glm_info = dir(fullfile(PATH_TO_DERIVATIVES,sprintf(subfolder,i),'eeg','GLM_*'));
    tmp = load(fullfile(PATH_TO_DERIVATIVES,sprintf(subfolder,i),'eeg',glm_info.name,'LIMO.mat'));
    X(:,i) = [tmp.LIMO.data.chanlocs.X];
    Y(:,i) = [tmp.LIMO.data.chanlocs.Y];
    Z(:,i) = [tmp.LIMO.data.chanlocs.Z];
end
expected_chanlocs.label = {tmp.LIMO.data.chanlocs.labels}'; %assumption: all subjects have same labels
expected_chanlocs.chanpos = [];
expected_chanlocs.chanpos(:,1) = mean(X,2);
expected_chanlocs.chanpos(:,2) = mean(Y,2);
expected_chanlocs.chanpos(:,3) = mean(Z,2);

expected_chanlocs = limo_get_ft_chanlocs(expected_chanlocs, defaults);
expected_chanlocs.expected_chanlocs = expected_chanlocs.chanlocs';
expected_chanlocs.pnt = expected_chanlocs.chanpos;
data_neighb = tmp.LIMO.data; %use the last subject to get the data structure
data_neighb.elec = expected_chanlocs;
cfg = [];
cfg.elec = expected_chanlocs;
cfg.neighbourdist = 40; %max distance between neighbours defined in cm in limo_ft_neighbourselection
[~,expected_chanlocs.channeighbstructmat] = limo_ft_neighbourselection(cfg,data_neighb);
end
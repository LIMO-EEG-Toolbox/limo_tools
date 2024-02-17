function average = limo_get_model_data(LIMO, regressor, extra, p, freq_index)

% short-cut from user provided options to limo_dipslay_results
% ----------------------------------------------------------------------
%  Copyright (C) LIMO Team 2024


if length(regressor) ~= 1
    error('This function takes one regressor as input');
end

if strcmpi(extra,'Original')

    Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
    index = find(LIMO.design.X(:,regressor));
    if strcmpi(LIMO.Analysis,'Time-Frequency')
        data     = Yr(channel,freq_index,:,index);
    else
        data     = Yr(channel,:,index);
    end
    average = nanmean(data,3);

elseif strcmpi(extra,'Modelled')
    Betas = load(fullfile(LIMO.dir,'Betas.mat'));
    Betas = Betas.Betas;
    if strcmpi(LIMO.Analysis,'Time-Frequency')
        Betas = squeeze(Betas(channel,freq_index,:,:));
    else
        Betas = Betas(channel,:,:);
    end
    for iChan = size(Betas,1):-1:1
        Yh(iChan,:,:) = (LIMO.design.X*squeeze(Betas(iChan,:,:))')'; % modelled data
    end
    index        = logical(LIMO.design.X(:,regressor));
    data         = Yh(:,:,index);
    average      = nanmean(data,3);
   
else % Adjusted
    allvar = 1:size(LIMO.design.X,2)-1;
    allvar(regressor)=[];
    if strcmpi(LIMO.Analysis,'Time-Frequency')
        Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
        Yr    = squeeze(Yr.Yr(channel,freq_index,:,:));
        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
        Betas = squeeze(Betas.Betas(channel,freq_index,:,:));
    else
        Yr    = load(fullfile(LIMO.dir,'Yr.mat'));
        Yr    = Yr.Yr;
        Betas = load(fullfile(LIMO.dir,'Betas.mat'));
        Betas = Betas.Betas;
    end
    for iChan = size(Betas,1):-1:1
        confounds(iChan,:,:) = (LIMO.design.X*squeeze(Betas(iChan,:,:))')'; % modelled data
    end
    Ya        = Yr - confounds; clear Yr Betas confounds;
    index     = logical(LIMO.design.X(:,regressor));
    data      = Ya(:,:,index);
    average   = nanmean(data,3);
end
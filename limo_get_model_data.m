function average = limo_get_model_data(LIMO, regressor, extra, p, freq_index)

if length(regressor) ~= 1
    error('This function takes one regressor as input');
end

probs = [p/2; 1-p/2];
z = norminv(probs);
categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
channel = ':';

if regressor <= categorical
    if strcmpi(LIMO.Analysis,'Time')
        mytitle = [ extra ' ERP'];
    elseif strcmpi(LIMO.Analysis,'Frequency')
        mytitle = [ extra ' Power Spectrum' ];
    elseif strcmpi(LIMO.Analysis,'Time-Frequency')
        mytitle = [ extra sprintf(' ERSP at %g Hz', frequency) ];
    end
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
%     se      = nanstd(data,0,3) ./ sqrt(numel(index));
%     ci      = bsxfun(@plus, average, se'*z);

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

    Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
    if strcmpi(LIMO.Analysis,'Time-Frequency')
        Yr = squeeze(Yr(:,freq_index,:,:));
    end
    R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));

    index        = logical(LIMO.design.X(:,regressor));
    data         = Yh(:,:,index);
    average      = mean(data,3);

%    R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
%     var          = diag(((R(index,index)*squeeze(Yr(channel,:,index))')'*(R(index,index)*squeeze(Yr(channel,:,index))')) / LIMO.model.model_df(2));
%     CI           = sqrt(var/size(index,1))*z';
%     ci           = (repmat(mean(data,2),1,2)+CI)';
    
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

%     se        = nanstd(data,0,3) ./ sqrt(numel(index));
%     ci        = repmat(average(i,:),2,1) + repmat(se',2,1).*repmat(z,1,size(Ya,1));
end
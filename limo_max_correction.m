function [mask,p_val] = limo_max_correction(M,bootM,p)

% Correction for multiple testing using the maximum stat value
% Since the FWER is te probability to make one error or more, then it is
% controlled by looking at the probability of the biggest stat value
%
% INPUTS: M = matrix of observed values (channel*[]*[]) 
%        (for a single channel/component the format is 1*[]*[])
%         bootM = matrix of F values for data under H0
%        (for a single channel/component the format is 1*[]*[]*MC)
%         p = threshold to apply
% 
% OUTPUTS: maks = a binary matrix of significant/non-significant cells
%          p_vals = a matrix of corrected p-values
%
% Note this function can be used for data or TFCE transformed data alike
%
% Cyril Pernet - outsourced from limo_stat_values
% ---------------------------------------------------------------------
% Copyright (C) LIMO Team 2016

ndim = numel(size(bootM));
if ndim ~=3 || ndim ~= 4
    error('H0 data must be 3D or 4D')
end

% simply collect highest absolute value across the while space for each MC
nboot = size(bootM,ndim);
parfor boot=1:nboot
    if ndim == 3
        data = squeeze(bootM(:,:,boot));
    elseif ndim == 4
        data = squeeze(bootM(:,:,:,boot));
    end
    maxM(boot) = max(data(:)); %
end

% using percentile to treshold the data
U=round((1-p).*nboot);
sortmaxM = sort(maxM); % sort bootstraps
maxF_th = sortmaxM(U); % get threshold for each parameter
mask = squeeze(M) >= maxF_th;

% compute p-values
p_val = NaN(size(M));

if ndim == 3
    for row =1:size(M,1)
        for column=1:size(M,2)
            p_val(row,column) = 1-(sum(squeeze(M(row,column)) >=sortmaxM) / nboot);
        end
    end
elseif ndim == 4
    for channel = 1:size(M,1)
        for freq = 1:size(M,2)
            for time = 1:size(M,3)
                p_val(channel,freq,time) = 1-(sum(squeeze(M(channel,freq,time)) >=sortmaxM) / nboot);
            end
        end
    end
end

% adjust p-values of zeros
p_val(p_val==0) = 1/nboot;

function [mainSub, gpSub, intSub] = limo_boot_repanova_sub(B, centered_data, boot_table, type, factor_levels, C, X, gp_vector, trimmed, nC)

% One bootstrap iteration (all channels) of the repeated-measures ANOVA under
% H0, returned as the per-iteration sub-arrays. Byte-for-byte the body of the
% limo_random_robust repeated-measures bootstrap parfor, so chunked
% accumulation reproduces the non-chunked H0 exactly.
%
% Returns
%   mainSub  [channels x frames x nC x 2]      (repeated-measures effect F/p)
%   gpSub    [channels x frames x 2]    or []   (group effect, type 3/4)
%   intSub   [channels x frames x nC x 2] or [] (interaction with gp, type 3/4)
% ------------------------------
%  Copyright (C) LIMO Team 2026

nchan   = size(centered_data,1);
nframes = size(centered_data,2);
array   = find(~isnan(centered_data(:,1,1,1)));

mainSub = NaN(nchan,nframes,nC,2);
gpSub   = [];
intSub  = [];
if type == 3 || type == 4
    gpSub  = NaN(nchan,nframes,2);
    intSub = NaN(nchan,nframes,nC,2);
end

for e = 1:length(array)
    channel = array(e);
    tmp = squeeze(centered_data(channel,:,boot_table{channel}(:,B),:));
    if size(centered_data,2) == 1
        Y  = ones(1,size(tmp,1),size(tmp,2)); Y(1,:,:) = tmp;
        gp = gp_vector(find(~isnan(Y(1,:,1))),:);
        Y  = Y(:,find(~isnan(Y(1,:,1))),:);
    else
        Y  = tmp(:,find(~isnan(tmp(1,:,1))),:);
        gp = gp_vector(find(~isnan(tmp(1,:,1))));
    end

    if type == 3 || type == 4
        XB = X(find(~isnan(tmp(1,:,1))));
    else
        XB = [];
    end

    if type == 1
        if trimmed, result = limo_robust_rep_anova(Y,gp,factor_levels,C);
        else,       result = limo_rep_anova(Y,gp,factor_levels,C); end
        mainSub(channel,:,1,1) = result.F;
        mainSub(channel,:,1,2) = result.p;
    elseif type == 2
        if trimmed, result = limo_robust_rep_anova(Y,gp,factor_levels,C);
        else,       result = limo_rep_anova(Y,gp,factor_levels,C); end
        mainSub(channel,:,:,1) = result.F';
        mainSub(channel,:,:,2) = result.p';
    elseif type == 3
        if trimmed, result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
        else,       result = limo_rep_anova(Y,gp,factor_levels,C,XB); end
        mainSub(channel,:,1,1) = result.repeated_measure.F;
        mainSub(channel,:,1,2) = result.repeated_measure.p;
        gpSub(channel,:,1)     = result.gp.F;
        gpSub(channel,:,2)     = result.gp.p;
        intSub(channel,:,1,1)  = result.interaction.F;
        intSub(channel,:,1,2)  = result.interaction.p;
    elseif type == 4
        if trimmed, result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
        else,       result = limo_rep_anova(Y,gp,factor_levels,C,XB); end
        mainSub(channel,:,:,1) = result.repeated_measure.F';
        mainSub(channel,:,:,2) = result.repeated_measure.p';
        gpSub(channel,:,1)     = result.gp.F;
        gpSub(channel,:,2)     = result.gp.p;
        intSub(channel,:,:,1)  = result.interaction.F';
        intSub(channel,:,:,2)  = result.interaction.p';
    end
end
end

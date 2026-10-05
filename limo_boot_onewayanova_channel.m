function [Fmat,pmat] = limo_boot_onewayanova_channel(cdata, boot_table, X_full, channel, b_range)

% Per-channel 1-way robust ANOVA bootstrap statistics for the requested
% bootstrap columns b_range, returned as [nframes x numel(b_range)].
% Byte-for-byte the limo_random_robust N-way (1-way robust) ANOVA bootstrap
% loop, so chunked accumulation reproduces the non-chunked H0 exactly.
%
% cdata      condition-centred data [channels x frames x subjects]
% boot_table cell array of resampling indices, one [nsubj x nboot] per channel
% X_full     the full design matrix LIMO.design.X (the intercept column,
%            i.e. the last column, is dropped here as in the original)
% ------------------------------
%  Copyright (C) LIMO Team 2026

nframes = size(cdata,2);
nb      = numel(b_range);
Fmat = NaN(nframes,nb); pmat = NaN(nframes,nb);

index = find(~isnan(squeeze(cdata(channel,1,:))));
X     = X_full(index,1:end-1);
if any(sum(X) == 0)     % a condition has no observation -> skip (leave NaN), as in the original
    return
end

bt = boot_table{channel};
Fc = cell(1,nb); pc = cell(1,nb);
parfor bi = 1:nb
    [Fc{bi}, pc{bi}] = limo_robust_1way_anova(squeeze(cdata(channel,:,bt(:,b_range(bi)))), X, 20);
end
for bi = 1:nb
    Fmat(:,bi) = Fc{bi}(:);
    pmat(:,bi) = pc{bi}(:);
end
end

function [tmat,pmat] = limo_boot_paired_channel(d1c, d2c, boot_table, channel, b_range, trimmed)

% Per-channel paired-samples bootstrap statistics for the requested bootstrap
% columns b_range, returned as [nframes x numel(b_range)]. The per-iteration
% computation is byte-for-byte the limo_random_robust paired t-test loop, so
% chunked accumulation reproduces the non-chunked H0 exactly.
%
% d1c,d2c   centred data1/data2 [channels x frames x subjects]
% boot_table cell array of resampling indices, one [nsubj x nboot] per channel
% channel    channel index to compute
% b_range    bootstrap columns to compute
% trimmed    true -> limo_yuend_ttest (Trimmed Mean); false -> limo_ttest (Mean)
% ------------------------------
%  Copyright (C) LIMO Team 2026

tmp = d1c(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:))));
tmp = d2c(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:))));
bt  = boot_table{channel};
nb  = numel(b_range);
nframes = size(Y1,2);
tcell = cell(1,nb); pcell = cell(1,nb);

if trimmed
    parfor bi = 1:nb
        [tcell{bi},~,~,~,pcell{bi},~,~] = limo_yuend_ttest(Y1(1,:,bt(:,b_range(bi))), Y2(1,:,bt(:,b_range(bi))));
    end
else
    parfor bi = 1:nb
        [~,~,~,~,~,tcell{bi},pcell{bi}] = limo_ttest(1, Y1(1,:,bt(:,b_range(bi))), Y2(1,:,bt(:,b_range(bi))));
    end
end

tmat = zeros(nframes, nb); pmat = zeros(nframes, nb);
for bi = 1:nb
    tmat(:,bi) = tcell{bi}(:);
    pmat(:,bi) = pcell{bi}(:);
end
end

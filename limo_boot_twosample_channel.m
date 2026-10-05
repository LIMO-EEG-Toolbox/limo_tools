function [tmat,pmat] = limo_boot_twosample_channel(d1c, d2c, bt1, bt2, channel, b_range, robust)

% Per-channel two-samples bootstrap statistics for the requested bootstrap
% columns b_range, returned as [nframes x numel(b_range)]. Byte-for-byte the
% limo_random_robust two-samples t-test bootstrap loop. Two boot tables are
% used (one per group), matching the non-chunked code.
%
% d1c,d2c   centred data1/data2 [channels x frames x subjects]
% bt1,bt2   resampling tables (cell per channel) for group 1 and group 2
% robust    true -> limo_yuen_ttest (Trimmed Mean / Welch); false -> limo_ttest (Mean)
% ------------------------------
%  Copyright (C) LIMO Team 2026

tmp = d1c(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:))));
tmp = d2c(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:))));
B1 = bt1{channel}; B2 = bt2{channel};
nb = numel(b_range);
nframes = size(Y1,2);
tcell = cell(1,nb); pcell = cell(1,nb);

if robust
    parfor bi = 1:nb
        [tcell{bi},~,~,~,pcell{bi},~,~] = limo_yuen_ttest(Y1(1,:,B1(:,b_range(bi))), Y2(1,:,B2(:,b_range(bi))));
    end
else
    parfor bi = 1:nb
        [~,~,~,~,~,tcell{bi},pcell{bi}] = limo_ttest(2, Y1(1,:,B1(:,b_range(bi))), Y2(1,:,B2(:,b_range(bi))), .05);
    end
end

tmat = zeros(nframes, nb); pmat = zeros(nframes, nb);
for bi = 1:nb
    tmat(:,bi) = tcell{bi}(:);
    pmat(:,bi) = pcell{bi}(:);
end
end

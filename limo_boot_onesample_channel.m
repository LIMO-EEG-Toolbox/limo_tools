function [tmat,pmat] = limo_boot_onesample_channel(cdata, boot_table, channel, b_range, trimmed)

% Per-channel one-sample bootstrap statistics for the requested bootstrap
% columns b_range, returned as [nframes x numel(b_range)]. Byte-for-byte the
% limo_random_robust one-sample t-test bootstrap loop.
%
% cdata     centred data [channels x frames x subjects]
% trimmed   true -> limo_trimci (Trimmed Mean); false -> limo_ttest (Mean)
% ------------------------------
%  Copyright (C) LIMO Team 2026

tmp = cdata(channel,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
bt  = boot_table{channel};
nb  = numel(b_range);
nframes = size(Y,2);
tcell = cell(1,nb); pcell = cell(1,nb);

if trimmed
    parfor bi = 1:nb
        [tcell{bi},~,~,~,pcell{bi},~,~] = limo_trimci(Y(1,:,bt(:,b_range(bi))));
    end
else
    parfor bi = 1:nb
        [~,~,~,~,~,tcell{bi},pcell{bi}] = limo_ttest(1, Y(1,:,bt(:,b_range(bi))), 0, 5/100);
    end
end

tmat = zeros(nframes, nb); pmat = zeros(nframes, nb);
for bi = 1:nb
    tmat(:,bi) = tcell{bi}(:);
    pmat(:,bi) = pcell{bi}(:);
end
end

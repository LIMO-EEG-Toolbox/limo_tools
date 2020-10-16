function newmask = limo_transprob(mask,th)

% routine to compute the transitional probability of binary mask
% (i.e. cluster of significant values)
%
% FORMAT: errors = limo_prederror(mask,th)
%
% INPUT: mask is a binary  image
%        th the threshold to split the clusters (default 50%)
%
% OUTPUT: mask splitted according to th
%
% Cyril Pernet - October 2020
% ------------------------------
%  Copyright (C) LIMO Team 2020

if nargin == 1
    th = 0.05;  % threshold prob. > 50% prior could be used to weight this
end

% how many elements transition from signitifcant to non-significant / number of significant elements
% --> cluster in shrinking in size
transitional_prob = zeros(1,size(mask,2));
newmask           = zeros(size(mask));
labels            = unique(mask(:));
label_onset       = 0;

for cluster = 1:length(labels)-1
    Nstates                         = nansum(mask==labels(cluster+1),1);
    state_change                    = diff(mask==labels(cluster+1));
    change_prob                     = (sum(state_change==1,1)+sum(state_change==-1,1)) ./ (2*Nstates);
    change_prob(isnan(change_prob)) = 0;
    newlabels                       = bwlabeln(change_prob>=th);
    tmpmask                         = zeros(size(mask));
    for p=find(Nstates)
        if newlabels(p) ~= 0
            tmpmask(find(mask(:,p)),p)  = newlabels(p)+ label_onset;  %#ok<FNDSB>
        end
    end
    newmask                         = newmask + tmpmask;
    transitional_prob               = transitional_prob + change_prob;
    label_onset                     = label_onset + max(newlabels);
end

figure; 
subplot(3,1,1); imagesc(LIMO.data.timevect,1:Nmask,mask); title('original mask');
subplot(3,1,2); plot(LIMO.data.timevect,transitional_prob);axis tight; title('transitional prob.');
subplot(3,1,3); imagesc(LIMO.data.timevect,1:Nmask,newmask); title('new mask');


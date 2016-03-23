[mask, cluster_p] = limo_chancluster_test(ori_f,ori_p,boot_f,boot_p,alphav);

% limo_chancluster_test finds clusters of significant F values, computes the
% sum of F values inside each cluster, and compares that sum to a threshold
% sum of F values expected by chance. This function works for low dimensional
% data - example of data in are ERP or power specrum vectors of F values, or
% ERSP (time*freq). 
%
% FORMAT: [mask, pval] = limo_chancluster_test(ori_f,ori_p,boot_f,boot_p,p)
%
% INPUTS: ori_f: vector or 2D matrix of observed F values 
%         ori_p: vector or 2D matrix of observed P values 
%         boot_f = 2D/3D matrix of f values for data under H0 
%         boot_p = 2D/3D matrix of f p values for data under H0 
%         alpha level, default 0.05
%
% OUTPUTS: mask = logical matrix of significant effects corrected for multiple comparisons by a cluster test
%          pval = corrected p values of significant clusters
%
% See also limo_ecluster_make limo_ecluster_test limo_tfcluster_make limo_tfcluster_test
% Cyril Pernet
% -------------------------------------------------
%  Copyright (C) LIMO Team 2016

if nargin == 4
    alphav = 0.05;
end

if (exist('spm_bwlabel','file')==3) + (exist('bwlabeln','file')==2) == 0
    errordlg('You need either the Image Processing Toolbox or SPM in your path to execute this limo_chancluster_test');
end

ndim = numel(size(ori_f));
if ndim ~=1 || ndim ~= 2
    error('data must be 1D or 2D')
end

if numel(size(ori_f)) ~= numel(size(ori_p)) || numel(size(boot_f)) ~= numel(size(boot_p))
    error('imput data for stat/p values are not of the same size')
end

if numel(size(boot_f)) ~= ndim+1
    error(['H0 data m ust be ' num2str(ndim+1) ' dimensional'])
end


%% compute clusters in time and/or freq domain under H0, take the max
b = size(boot_f,ndim+1);
boot_values = zeros(b,1);

for kk=1:b % bootstrap samples
    
    % get cluster along 1st dim
    if ndim ==1
        if exist('spm_bwlabel','file') ==3
            [L,NUM] = spm_bwlabel(squeeze(boot_p(:,kk))<=alphav,6);
        elseif exist('bwlabeln','file')==2
            [L,NUM] = bwlabeln(squeeze(boot_p(:,kk))<=alphav);
        end
        
        if NUM~=0
            tmp=zeros(1,NUM);
            for C = 1:NUM % compute sum for each cluster
                tmp(C) = sum( squeeze(boot_f(L==C,kk)));
            end
            boot_values(kk) = max(tmp(:)); % save max across clusters
        else
            boot_values(kk) = 0;
        end
        
    elseif ndim == 2
        if exist('spm_bwlabel','file') ==3
            [L,NUM] = spm_bwlabel(squeeze(boot_p(:,:,kk))<=alphav,6);
        elseif exist('bwlabeln','file')==2
            [L,NUM] = bwlabeln(squeeze(boot_p(:,:,kk))<=alphav);
        end
        
        if NUM~=0
            data = squeeze(boot_f(:,:,kk));
            tmp=zeros(1,NUM);
            for C = 1:NUM % compute sum for each cluster
                tmp(C) = sum(data(L==C));
            end
            boot_values(kk) = max(tmp(:)); % save max across clusters
        else
            boot_values(kk) = 0;
        end
    end    
end % bootstrap loop

sortSC = sort(boot_values);
U = round((1-alphav)*b);
th = sortSC(U);

%% compute clusters in time and/or freq domain and compare to the max distribution
if ndim ==1
    if exist('spm_bwlabel','file') ==3
        [L,NUM] = spm_bwlabel(squeeze(ori_p(:,kk))<=alphav,6);
    elseif exist('bwlabeln','file')==2
        [L,NUM] = bwlabeln(squeeze(ori_p(:,kk))<=alphav);
    end
elseif ndim == 2
    if exist('spm_bwlabel','file') ==3
        [L,NUM] = spm_bwlabel(squeeze(ori_p(:,:,kk))<=alphav,6);
    elseif exist('bwlabeln','file')==2
        [L,NUM] = bwlabeln(squeeze(ori_p(:,:,kk))<=alphav);
    end
end

corrected_p = NaN(size(ori_p));
mask = corrected_p;

for C = 1:NUM % compute cluster sums & compare to bootstrap threshold
    corrected_p((L==C)) = 1-sum(sum(ori_f(L==C))>=th)/b;
    if corrected_p((L==C)) < alphav
        mask(L==C)=1; % flag clusters above threshold
    end
end

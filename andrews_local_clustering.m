function [mask,cluster_p] = andrews_local_clustering(M,P,bootM,bootP,LIMO,MCC,p)
% optimised, slightly 'skinnier' local_clustering function, without call field trip functions

% M = 2D matrix of observed F values (note for a single electrode the format is 1*time frames*trials)
% P = 2D matrix of observed p values (note for a single electrode the format is 1*time frames*trials)
% bootM = 3D matrix of F values for data bootstrapped under H0
% bootP = 3D matrix of F values for data bootstrapped under H0
% LIMO = LIMO structure - information requested is LIMO.data.chanlocs and LIMO.data.neighbouring_matrix
% MCC = 2 (spatial-temporal clustering) or 3 (temporal clustering)
% p = threshold to apply (note this applied to create clusters and to
% threshold the cluster map)

if size(M,1) == 1
    MCC = 3;
end
cluster_p = [];
mask = [];

if MCC == 2 
    nboot = size(bootM,3);
    U = round((1-p)*nboot); % bootstrap threshold
    if size(bootM,1)>1 % many electrodes
        minnbchan = 2;
        expected_chanlocs = LIMO.data.chanlocs;
        channeighbstructmat = LIMO.data.neighbouring_matrix;
        boot_maxclustersum=zeros(nboot,1); % compute bootstrap clusters
        parfor boot=1:nboot
            %boot_maxclustersum(boot) = andrews_getclustersum(bootM(:,:,boot),bootP(:,:,boot),channeighbstructmat,minnbchan,p);
            p_under_alpha = bootP(:,:,boot) <= p;
            bootM_f = bootM(:,:,boot);
            [posclusterslabelmat,nposclusters] = andrews_findcluster(p_under_alpha,channeighbstructmat,minnbchan);
            
            if nposclusters~=0
                tmp=zeros(1,nposclusters);
                for C = 1:nposclusters % compute sum for each cluster
                    tmp(C) = sum( bootM_f(posclusterslabelmat==C) );
                end
                boot_maxclustersum(boot) = max(tmp(:)); % save max across clusters
            else
                boot_maxclustersum(boot) = 0;
            end
            
            boot_progress_msg_at = [1 25 50 75 100];
            if any(boot==boot_progress_msg_at*10)
                boot_progress_txt = sprintf('Currently on bootstrap %d of %d ...',boot,nboot);
                disp(boot_progress_txt)
            end
            
        end
        [mask, cluster_p] = limo_cluster_test(M,P,boot_maxclustersum,channeighbstructmat,minnbchan,p);
        
    elseif size(bootM,1)==1 % one electrode
        th = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);
        sigcluster = limo_ecluster_test(squeeze(M),squeeze(P),th,p);
        mask = sigcluster.elec; cluster_p = [];
    end
    
elseif MCC == 3
    nboot = size(bootM,3);
    U = round((1-p)*nboot); % bootstrap threshold
    th = limo_ecluster_make(squeeze(bootM),squeeze(bootP),p);
    sigcluster = limo_ecluster_test(squeeze(M),squeeze(P),th,p);
    mask = sigcluster.elec;
    
end
end


function [mask,p_val] = max_correction(M,bootM,p)
% correction for multiple testing using the max stat value
% note this works for bootstrapped data under H0 and for TFCE
%
% M = 2D matrix of observed values (note for a single electrode the format is 1*time frames*trials)
% bootM = 3D matrix of F values for data bootstrapped under H0
% p = threshold to apply 

nboot = size(bootM,3);
for boot=1:nboot
    data = squeeze(bootM(:,:,boot));
    maxM(boot) = max(data(:)); % collect highest absolute value in space and time for each boot
end

U=round((1-p).*nboot);
sortmaxM = sort(maxM); % sort bootstraps 
maxF_th = sortmaxM(U); % get threshold for each parameter
mask = squeeze(M) >= maxF_th; 
% figure; imagesc(mask)
for row =1:size(M,1)
    for column=1:size(M,2)
        p_val(row,column) = 1-(sum(squeeze(M(row,column)) >=sortmaxM) / nboot);
    end
end 
 
end
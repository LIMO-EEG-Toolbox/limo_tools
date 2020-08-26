function th = limo_ecluster_make(bootf,bootp,alpha)
% function th = limo_ecluster_make(bootf,bootp,alpha)
%
% ECLUSTER_MAKE creates a multivariate, cluster based, statistical threshold based on bootstrapped F tests.
% The function operates along the last dimension.
% For each bootstrap, independently at each electrode:
%       significant F values are clustered in time;
%       the sum of F values inside each cluster is computed;
%       the maximum sum across clusters is saved.
% The maximum sums of clusters of significant F values are then sorted to
% obtain a (1-alpha)% percentile threshold.
% NOTE: for a two-tailed bootstrap t-test, enter a squared matrix of T values:
% th = ecluster_make(boott.^2,bootp,alpha)
%
% INPUTS:
%           BOOTF = 2D or 3D matrix of F values with format electrode x
%               frames x bootstrap resamples or frames x bootstrap resamples;
%           BOOTP = matrix of p values associated with bootf, with same
%               format;
%           ALPHA = type I error rate, default 0.05
%
% OUTPUTS:
%           TH.ELEC = cluster threshold at each electrode
%           TH.MAX = max cluster threshold across electrodes; this is a
%               more conservative way to control for multiple comparisons 
%               than using a spatial-temporal clustering technique.
%
% v1 Guillaume Rousselet, University of Glasgow, August 2010
% Luisa Frei, 4 Nov 2011: fixed bug in all electrode ouput
% edit Marianne Latinus adding spm_bwlabel
% -----------------------------
%  Copyright (C) LIMO Team 2010
%
% See also LIMO_ECLUSTER_TEST

if nargin < 3
    alpha = 0.05;
end

if ndims(bootf)==3 % Ne x Nf **************************

    b = size(bootf,3);
    U = round((1-alpha)*b);
    Ne = size(bootf,1);
    SC = zeros(b,Ne);

    for E = 1:Ne % electrodes

        for kk=1:b % bootstrap samples

            %[L,NUM] = bwlabeln( squeeze(bootf(E,:,kk)) .* (squeeze(bootp(E,:,kk))<=alpha) ); % find clusters
            try
                [L,NUM] = bwlabeln( squeeze(bootf(E,:,kk)) .* (squeeze(bootp(E,:,kk))<=alpha) );
            catch ME
                try
                    [L,NUM] = spm_bwlabel( squeeze(bootf(E,:,kk)) .* (squeeze(bootp(E,:,kk))<=alpha), 6);
                catch ME
                    errordlg('You need either the Image Processing Toolbox or SPM in your path');
                end
            end
            if NUM~=0
                tmp=zeros(1,NUM);
                for C = 1:NUM % compute sum for each cluster
                    tmp(C) = sum( squeeze(bootf(E,L==C,kk)) );
                end
                SC(E,kk) = max(tmp(:)); % save max across clusters
            else
                SC(E,kk) = 0;
            end

        end % bootstrap loop

    end % electrode loop

    sortSC = sort(SC,2);
    th.elec = sortSC(:,U); % threshold at each electrode

    maxSC = max(SC,[],1); % max across electrodes
    sortmaxSC = sort(maxSC);
    th.max = sortmaxSC(U); % threshold of max across electrodes

elseif ndims(bootf)==2 % Nf only, no electrode dimension **************************

    b = size(bootf,2);
    U = round((1-alpha)*b);
    SC = zeros(b,1);

    for kk=1:b % bootstrap samples

        [L,NUM] = bwlabeln( squeeze(bootf(:,kk)) .* (squeeze(bootp(:,kk))<=alpha) ); % find clusters

        if NUM~=0
            tmp=zeros(1,NUM);
            for C = 1:NUM % compute sum for each cluster
                tmp(C) = sum( squeeze(bootf(L==C,kk)) );
            end
            SC(kk) = max(tmp(:)); % save max across clusters
        else
            SC(kk) = 0;
        end

    end % bootstrap loop

    sortSC = sort(SC);
    th.elec = sortSC(U); % threshold at the unique electrode

else

    error('ecluster_make: bootf should be 2D or 3D')

end




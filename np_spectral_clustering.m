function [p_values, cluster_pvalues, cluster_mass, FDR_threshold] = ...
    np_spectral_clustering(spectra1, spectra2, varargin)

% Compute Non-Parametric (i.e. rank based) tests between spectra1 and
% spectra2 at each frequencies - correction for multiple testing is
% obtained by clustering consecutive frequency points - it returns the 
% cluster p values corrected for multiple comparisons using type 1 error
% family-wise error correction based on maximum statistics under the null
% (see Pernet et al., 2015 J Neurosci Methods 250:85-93. 
% doi: 10.1016/j.jneumeth.2014.08.003) -- for comparison it also returns
% the non-paramteric False Discovery Rate threshold
%
% FORMAT [p_values, cluster_pvalues, cluster_mass, mask] = np_spectral_clustering(spectra1, spectra2, 'type', value)
%        [~,cluster_pvalues] = np_spectral_clustering(spectra1, spectra2, 'type', value, 'tail', value, 'threshold', value, 'figure', value)
%
% INPUTS spectra1 and spectra2 are matrices of the same length (row are
%                              observations, columns are frequencies)
%       'type' must be either 'paired' or 'independent' 
%       'tail' must be 'both' (default) | 'right' | 'left'
%       'threshold' is the p_value threshold i.e. 0.05 (default)
%       'figure' is 'new' (default), 'on' (plot in current one) or 'off' 
%       'verbose' is 'on' by default
%
% OUTPUTS p_values is a vector of p values for the Wilcoxon rank sign test
%                  (paired data) or Wilcoxon rank sum test (independent data).
%         cluster_pvalues are the p values of the clusters corrected for
%                         multiple testing based on type 1 FWE control 
%         cluster_mass are the cluster masses (sum of ranks)
%         mask is a N-ary vector indicating significant clusters
%
% see also signrank, ranksum, bwlabeln
%
% Cyril Pernet - September 2023

%% set defaults
Ncores    = [];     % how nany cores to run the bootstrap default N-1
Nboot     = 1000;   % 1000 boostraps by default
tail      = 'both';
threshold = 0.05;
fig       = 'new';
verbose   = 'on';

%% check data in

if ~ismatrix(spectra1) || ~ismatrix(spectra2)
    error('data in are not recognize as matrices')
else
    [n1,p1] = size(spectra1);
    [n2,p2] = size(spectra2);
    if p1 ~= p2
        error('spectra must be of the same length - s1=%g, s2=%g',p1,p2);
    end
end

% update
for arg = 1:2:nargin-2
    if strcmpi(varargin{arg}, 'type')
        type = varargin{arg+1};
        if contains(type,'indep')
            type = 'ranksum';
        elseif contains(type,'pair')
            type = 'signrank';
        end
    elseif strcmpi(varargin{arg}, 'tail')
        tail = varargin{arg+1};
    elseif strcmpi(varargin{arg}, 'threshold')
        threshold = varargin{arg+1};
    elseif strcmpi(varargin{arg}, 'figure')
        fig = varargin{arg+1};
    elseif strcmpi(varargin{arg}, 'verbose')
        verbose = varargin{arg+1};
    else
        warning('unrecognized argument %s',varargin{arg})
    end
end

if  ~any(strcmp(type,{'signrank','ranksum'}))
    error('type argument must be paired or independent')
end

if strcmp(type,'signrank')
    if n1 ~= n2
        error('for paired data, spectra 1 (n=%g) and spectra 2 (n=%g) must have the same number of observations/rows',n1,n2)
    end
end

if threshold > 1
    threshold = threshold/100;
    warning('threshold set to %g',threshold);
end

if  ~any(strcmp(fig,{'new', 'on', 'off'}))
    warning('figure argument mnot recognized - using default ''new''')
end

%% Compute statistics for observed values

for freq = p1:-1:1
    if strcmp(type,'ranksum') % independent data
        [p_values(freq),~,stats] = ranksum(spectra1(:,freq),spectra2(:,freq),...
            'tail',tail); test_rank(freq) = stats.ranksum;
    elseif strcmp(type,'signrank') % paired data
        [p_values(freq),~,stats] = signrank(spectra1(:,freq),spectra2(:,freq),...
            'tail',tail); test_rank(freq) = stats.signedrank;
    end
end

%% Compute statistics under the null

% check parallel computing
p = gcp('nocreate');
if isempty(p)
    c            = parcluster;
    if ~isempty(Ncores)
        c.NumWorkers = Ncores;
    else
        Ncores = c.NumWorkers;
    end
    saveProfile(c);
    try
        parpool(Ncores-1);
    catch flasestart %#ok<NASGU>
        delete(c.Jobs)
        parcluster('local')
        parpool(Ncores-1);
    end
end

% resample under H0
if strcmp(type,'ranksum') % independent data
    % tests the null hypothesis that data in x and y have equal medians
    null_spectra1 = spectra1 - median(spectra1, 1,'omitnan');
    resample1     = randi(n1,n1,Nboot);
    null_spectra2 = spectra2 - median(spectra2, 1,'omitnan');
    resample2     = randi(n2,n2,Nboot);
    parfor b = 1:Nboot
        if strcmp(verbose,'on')
            fprintf('Computing the null %g/%g\n',b,Nboot)
        end
        for freq = p1:-1:1
            [p_null(b,freq),~,stats] = ranksum(null_spectra1(resample1(:,b),freq),...
                null_spectra2(resample2(:,b),freq), 'tail',tail); 
            null_rank(b,freq) = stats.ranksum;
        end
    end
elseif strcmp(type,'signrank') % paired data
    % tests the null hypothesis that x – y has zero median
    % since signrank(x,y) = signrank(x-y), the null is (x-y)-median(x-y)
    null_data = (spectra1 - spectra2) - median((spectra1-spectra2), 1,'omitnan');
    resample = randi(size(null_data,1),size(null_data,1),Nboot);
    parfor b = 1:Nboot
        if strcmp(verbose,'on')
            fprintf('Computing the null %g/%g\n',b,Nboot)
        end
        for freq = p1:-1:1
            [p_null(b,freq),~,stats] = signrank(null_data(resample(:,b),freq));
            null_rank(b,freq) = stats.signedrank;
        end
    end
end

%% Correct using clustering and max

% 1st get the distribution of maxima under H0
for b=1:Nboot 
    [L,NUM] = bwlabeln(squeeze(p_null(b,:))<=threshold); % find clusters
    if NUM~=0
        tmp=zeros(1,NUM);
        for C = 1:NUM % compute sum for each cluster
            tmp(C) = sum(squeeze(null_rank(b,L==C)));
        end
        boot_values(b) = max(tmp(:)); % save max across clusters
    else
        boot_values(b) = 0;
    end
end 

boot_values                     = sort(boot_values);
boot_values(isnan(boot_values)) = [];
U                               = round((1-threshold)*length(boot_values)); % percentile
boot_threshold                  = boot_values(U); % threshold at the unique electrode
% hist(boot_values)

% 2nd threshold observed data
[L,NUM] = bwlabeln(p_values<=threshold); % find clusters

cluster_mass    = zeros(1,NUM);
mask            = zeros(1,p1);
cluster_pvalues = nan(1,p1);
cluster_label   = 1;
for C = 1:NUM % compute cluster sums & compare to bootstrap threshold
    cluster_mass(C) = sum(squeeze(test_rank(L==C)));
    if cluster_mass(C) >= boot_threshold
        mask(L==C)    = cluster_label; % flag clusters above threshold
        cluster_label =  cluster_label+1;
        if ~isempty(boot_values)
            p = sum(cluster_mass(C) >= boot_threshold)/length(boot_values);
            P = min(p,1-p);
            if p ==0
                p = 1/length(boot_values);
            end
            cluster_pvalues(L==C) = p;
        end
    end
end

%% figure
if strcmpi(fig,'new')
    figure('Name','np_spectral_clustering')
end

if ~strcmpi(fig,'off')
    % plot the data and significant frames
    subplot(1,2,1); plot(mean(spectra1,1),'LineWidth',2); grid on
    hold on; plot(mean(spectra2,1),'LineWidth',2); 
    tmp = mask; tmp(mask==0) = NaN; plot(tmp,'r*'); 
    title('Observed data'); axis('tight'); box on

    % plot the null and threshold
    subplot(1,2,2); plot(boot_values,'LineWidth',3); grid on; hold on;     
    plot(min(find(boot_values==boot_threshold)),boot_threshold,'r*','LineWidth',5)
    txt = ['bootstrap threashold ' num2str(boot_threshold) '\rightarrow'];
    text(min(find(boot_values==boot_threshold)),boot_threshold,txt,'FontSize',10,'HorizontalAlignment','right');
    
    if NUM ~=0
        [val,loc]=min(abs(boot_values-max(cluster_mass)));
        plot(loc,max(cluster_mass),'r*','LineWidth',5)
        txt = ['biggest observed cluster mass: ' num2str(max(max(cluster_mass))) '\rightarrow'];
        if loc < Nboot/2
            text(loc,double(max(max(cluster_mass))),txt,'FontSize',10,'HorizontalAlignment','left');
        else
            text(loc,double(max(max(cluster_mass))),txt,'FontSize',10,'HorizontalAlignment','right');
        end
    end

    title('Cluster-mass Maxima under H0','FontSize',12)
    xlabel('sorted bootstrap iterations','FontSize',12); 
    ylabel('Freq.','FontSize',12)
    box on; set(gca,'Layer','Top')
end

%% non-parametric FDR 

p   = sort(p_values);
V   = length(p);
I   = (1:V);
cVN = sum(1./(1:V));
FDR_threshold = p(max(find(p<=I/V*threshold/cVN)));
if isempty(FDR_threshold)
    FDR_threshold = 0; 
end

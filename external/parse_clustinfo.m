% parse_clustinfo() - parse cluster info in STUDY structure
%
% Usage:
%   >>  [clust_stat,SubjClusIC_Matrix] = parse_clustinfo(STUDY, parentclust);
%
% Inputs:
%      STUDY         - studyset structure containing some or all files in ALLEEG
%      parentclust   - String with the name of the parentcluster who
%                      comprise the clusters of interest. If no 'parentclust' is provided,
%                      wil use the first Parentcluster by default.
%    
% Outputs:
%    clust_stat        - Structure with one substructure for each cluster under
%                        the parentcluster. Each substructure have the following info:                                      
%                        cluster's name
%                        And then a paired info of Subject/dataset name, IC's,
%                        1st and 2nd independet variables.
%   SubjClusIC_Matrix  - 3D binary matrix of incidence of Clusters Vs Datasets/Subjects Vs IC's
%
% See also: 
%   std_infocluster
%
% Author: Ramon Martinez-Cancino, SCCN, 2014
%
% Copyright (C) 2014  Ramon Martinez-Cancino, 
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [clust_stat,SubjClusIC_Matrix] = parse_clustinfo(STUDY,parentcluster)

hits_temp = cellfun(@(x)strcmp(x,deblank(parentcluster)),{STUDY.cluster.name});
parent_indx = find(hits_temp);
if isempty(parent_indx)
    frpintf('parse_clustinfo error: Invalid input parentcluster');
    return
end

nsets =  length(STUDY.datasetinfo);
SubjClusIC_Matrix = zeros(length(STUDY.cluster(parent_indx).child),nsets,max(STUDY.cluster(parent_indx).comps));
 
%--------------------------------------------------------------------------

% Getting cls
for i = 1:length({STUDY.cluster.name})
    tmpval = STUDY.cluster(i).parent;
    if isempty(tmpval)
        hits_tmp(i) = 0;
    else
        hits_tmp(i) = strcmp(tmpval{1},deblank(parentcluster));
    end
end
cls = find(hits_tmp); clear hits_tmp tmpval hits_tmp;

for icls = 1:length(cls)
    store_sets = STUDY.cluster(cls(icls)).sets';
    store_sets = store_sets(:);
    
    store_ics = repmat(STUDY.cluster(cls(icls)).comps,1,size(STUDY.cluster(cls(icls)).sets,1));
    
    for i = 1:size(store_sets,1)
        SubjClusIC_Matrix(icls,store_sets(i),store_ics(i)) = 1;
    end
    
    % Creating output structure
    clust_stat.clust(icls).name          = STUDY.cluster(cls(icls)).name;
    clust_stat.clust(icls).subj          = {STUDY.datasetinfo(store_sets).subject};
    clust_stat.clust(icls).datasetinfo   = {STUDY.datasetinfo(store_sets)};
    clust_stat.clust(icls).ics           = store_ics;
end
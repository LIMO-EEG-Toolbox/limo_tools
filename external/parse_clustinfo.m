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

hits_temp = cellfun(@(x)strcmp(x,parentcluster),{STUDY.cluster.name});
parent_indx = find(hits_temp);


nsets =  length(STUDY.datasetinfo);
SubjClusIC_Matrix = zeros(length(STUDY.cluster(parent_indx).child),nsets,max(STUDY.cluster(1).comps));

dsgn = STUDY.currentdesign;

%--------------------------------------------------------------------------
cellval = {STUDY.design(dsgn).cell.value}';

var1_label = STUDY.design(dsgn).variable(1).label;
var2_label = STUDY.design(dsgn).variable(2).label;

var1_uniquevals = STUDY.design(dsgn).variable(1).value;
var2_uniquevals = STUDY.design(dsgn).variable(2).value;

var1_dataval = cellfun(@(x) x(1),cellval);
var2_dataval = cellfun(@(x) x(2),cellval);
%--------------------------------------------------------------------------

ivar1 = length(STUDY.design(dsgn).variable(1).value);
ivar2 = length(STUDY.design(dsgn).variable(2).value);

c1 = 1;
for cls = (parent_indx+1):(parent_indx + length(STUDY.cluster(parent_indx).child))
    c2 = 1;
    for nv1 = 1:ivar1
        for nv2 = 1:ivar2
            
            allinds   = STUDY.cluster(1,cls).allinds{nv1,nv2};
            setinds = STUDY.cluster(1,cls).setinds{nv1,nv2};
            tmpset = cell2mat({(STUDY.design(dsgn).cell(setinds).dataset)})';
            for i = 1:length(tmpset)
                SubjClusIC_Matrix(cls-1,tmpset(i),allinds(i)) = 1;
                store_sets(c2) = tmpset(i);
                store_ics(c2)   = allinds(i);
                c2 = c2+1;
            end
        end
    end
    % Creating output structure
    clust_stat.clust(c1).name     = STUDY.cluster(cls).name;
    clust_stat.clust(c1).subj     = {STUDY.datasetinfo(store_sets).subject};
    clust_stat.clust(c1).ics      = store_ics;
    clust_stat.clust(c1).iv1      = var1_dataval(store_sets)';
    clust_stat.clust(c1).iv2      = var2_dataval(store_sets)';
    
    store_sets = 0;
    store_ics = 0;
    c1 = c1 +1;
end
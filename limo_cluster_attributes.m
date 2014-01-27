function [extend_mask, height_mask, extend, height, H0_extend, H0_height] = limo_cluster_attributes(data,dataH0,neighbouring_matrix,p_value)
% data is 3D [electrode,time frames,[F,p])
% dataH0 is 4D [electrode,time frames,[F,p],resamples)
% neighbouring_matrix is the matrix describing which electrodes are neighbour
% p_value is the cluster forming threshold and min p value of the cluster

% observed data
F_map = squeeze(data(:,:,1));
p_map = squeeze(data(:,:,2));
[clustered_map, num] = limo_ft_findcluster((p_map < p_value), neighbouring_matrix,2);

extend_map = zeros(size(F_map)); % same as cluster map but contains extend value instead
height_map = zeros(size(F_map)); % same as cluster map but contains height value instead
for i=1:num
    extend(i) = sum(clustered_map(:)==i);
    extend_map = extend_map + (clustered_map == i).*extend(i); % not a 'true' sum since there is no overlap
    height(i) = max(max(F_map.*(clustered_map == i)));
    height_map = height_map + (clustered_map == i).*height(i);
end


% under H0 - all we want is the distribution of max
nboot = size(dataH0,4);
index = 1;
for b=1:nboot
    fprintf('processing bootstrap %g \n',b);
    [clustered_map, num] = limo_ft_findcluster((squeeze(dataH0(:,:,2,b)) < p_value), neighbouring_matrix,2);
    if num ~=0
        for i=1:num
            tmp_extend(i) = sum(clustered_map(:)==i);
            tmp_height(i) = max(max(squeeze(dataH0(:,:,1,b)).*(clustered_map == i)));
        end
        H0_extend(index) = max(tmp_extend);
        H0_height(index) = max(tmp_height);
        index = index+1; clear tmp_extend tmp_height
    end
end


% now we can create a mask to threshold the observed data
U = round((1-p_value)*nboot);
H0_extend = sort(H0_extend);
th_extend = H0_extend(U);
extend_mask = extend_map(:,:) > th_extend;
H0_height = sort(H0_height);
th_height = H0_height(U);
height_mask = height_map(:,:) > th_height;


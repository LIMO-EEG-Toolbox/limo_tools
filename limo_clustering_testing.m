function limo_clustering_testing

% routine to check that the limo_clustering functions do what they are
% supposed to -- it simulates a matrix 10 electrodes * 4 freqencies * 6
% time points. Electrodes 2 4 6 8 are neighbours and various amont of 
% signal is introduced there, the rest is 0. Each time one tests if the
% mask of 'significant' values corresponds to the data
% 
% Cyril Pernet 16 May 2014
% -----------------------------
%  Copyright (C) LIMO Team 2014

NullData = zeros(10,4,6);
Neighbour = zeros(10,10);
Neighbour(2,4) = 1; Neighbour(2,6) = 1; Neighbour(2,8) = 1;
Neighbour(4,2) = 1; Neighbour(4,6) = 1; Neighbour(4,8) = 1;
Neighbour(6,2) = 1; Neighbour(6,4) = 1; Neighbour(6,8) = 1;
Neighbour(8,2) = 1; Neighbour(8,4) = 1; Neighbour(8,6) = 1;
figure; imagesc(Neighbour); colormap('gray'); title('Neighbourhood matrix')

elec_display = [1 3 5 7 9 10 2 4 6 8]; % indices to look at data
elec_signal = [2 4 6 8]; % indices to input signal

%% signal = 1 over 2 frequencies and 3 time points
Data = NullData;
Data(elec_signal,2:3,3:5) = 1; % 4 elect 2 freq 3 time
 
% let's do an ERP one one electrode, we call limo_ecluster_test
one_d = squeeze(Data(2,2,:))';
expected_cluster_mass = sum(one_d(:));
th.elec = expected_cluster_mass; % threshold for the cluster
sigcluster = limo_ecluster_test(one_d,(1-one_d),th,0.5);
if sum(sigcluster.elec == one_d) == numel(one_d)
    fprintf('Clustering of 1D data is ok \n')
else
    fprintf('ERROR Clustering of 1D data failed\n')
end

% let's do an ERP over several electrodes, we can call limo_ecluster_test
% or call limo_cluster_test

two_d = squeeze(Data(2,:,:));
expected_cluster_mass = sum(two_d(:));

% let's do an ERSP over several electrodes, we can call limo_ecluster_test
% or call limo_cluster_test (which calls limo_ft_findcluster)

[posclusterslabelmat,nposclusters] = limo_ft_findcluster(Data,Neighbour,2);
test = posclusterslabelmat == Data;
if sum(test(:)) == numel(Data)
    fprintf('Clustering of 3D data is ok \n')
else
    fprintf('ERROR Clustering of 3D data failed\n')
end

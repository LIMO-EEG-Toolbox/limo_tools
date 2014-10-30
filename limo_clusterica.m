function IDX = eeg_clusterica(data,fields)

% data cell array with one cell per measure
% fields
% IDX is the indices of the compoments after clustering (same output format as K-means)

% ERP [ic,time] abs(xcorr)
% Spect [ic, freqs] xcorr
% ERSP [ic, freqs, time] xcorr
% ITC [ic, freqs, time] abs(xcorr)
% dip [x,y,z,xx,yy,zz] squareform(pdist([x y z],'euclidean')) and normalize to 1
%                      abs(dot([xx,yy,zz]1,[xx,yy,zz]2))
%                      average the distance matrix with orientation matrix
% Scalp map [64*64] abs(corr)
% we use abs values to be insensitive to polarity inversion - not for spec and ersp
% because spec and ersp show power, does  not depend on polarity
%


base = pwd;
list = dir;
index = 1; 
for n=3:size(list,1)
    if list(n).isdir == 1 && strncmp(list(n).name,'S',1)
        dataname.erp{index} = [base filesep list(n).name filesep 'Probe.set_single_trials.icaerp'];
        index = index+1;
    end
end

% dataname.erp = cell array of icaerp files
% dataname.ersp = the icaersp files
% dataname.spec = the icaspec files
% dataname.dip = the icadip files


if isfield(dataname,'erp')
    % load data
    for n=1:length(dataname.erp)
        signal{n} = load('-mat',dataname.erp{n});
    end       
    
    
    % compute the similarity matrix
    
    
    % cluster
    
    
end


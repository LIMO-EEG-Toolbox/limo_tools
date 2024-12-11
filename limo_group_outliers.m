function [global_weights, local_weights, errors, outliers] = limo_group_outliers(Beta_files,expected_chanlocs)

% FORMAT
%
% INPUTS LIMO_files the list of LIMO-mat files of subjects to analyze
%        neighbourgh_matrix a binary of wich channels are neighbourgh
%
% OUTPUTS
%
%
%

%% go fetch beta
if ischar(expected_chanlocs)
        LIMO.data = load(expected_chanlocs);
        if isfield(LIMO.data,'expected_chanlocs')
            LIMO.data.chanlocs = LIMO.data.expected_chanlocs;
        end
        if isfield(LIMO.data,'channeighbstructmat')
            LIMO.data = renameStructField(LIMO.data, 'channeighbstructmat', 'neighbouring_matrix');
        end
    else
        data.chanlocs            = expected_chanlocs.expected_chanlocs;
        data.expected_chanlocs   = expected_chanlocs.expected_chanlocs;
        data.neighbouring_matrix = expected_chanlocs.channeighbstructmat;
 end

% load each beta, reading the Beta list
% use the neighbour matrix to mmake sure channels match
% fullfile(LUIMO_files.dir,LIMO_files.name)
% LIMO_files = dir(fileparts(Beta_files),'LIMO*.txt');

% to do cyril
% for subject = 1:N
%     % load LIMO.mat of the subject use LIMO.data.chanloc
%     data(:,:,:,subject) = limo_match_elec(LIMO.data.chanloc,data.expected_chanlocs,a_beg,a_end,betas)
% end



%% Autoencoder -- arguments in: matrix of beta values, neighbourgh_matrix
%             -- argument out: learned matrix of beta values

learned_betas = pyrunfile("NiPyAEoutliers.py", learned_betas, ...
    datain = beta_values, binatry_matrix = neighbourgh_matrix);


%% Error and Weight
% load Yhat, computed learned_Yhat, get error






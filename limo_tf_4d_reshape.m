function reshaped = limo_tf_4d_reshape(reshape_in,forced_dim)

% Simple function to 'stack' elec-freqs-times-N or 'unstack' 
% elec-freqs*times-N matrices
%
% FORMAT reshaped = limo_tf_4d_reshape(reshape_in)
%
% INPUT/OUTPUT 4D elec-freqs-times-N / 3D elec-freqs*times-N (all freq at a given time point)
%              3D elec-freqs*times-N / 4D elec-freqs-times-N
%              
% By default the function looks for the LIMO variable (either already used
% as global, or in the current dir or above) -- the forced_dim variable can
% be used instead of the LIMO structure, forcing a set of dimensions. If
% the input matrix is 4D, forced dim must be 3D, if the input matrix if 3D,
% the forced dim must be 4D 
%
% The reshaping is done explicitly within simple loops to be clear to read.
% Could be vectorised and/or use reshape(), but that makes it easier to get
% lost in dimensions.
%
% ------------------------------------------------------------------------
%  Copyright (C) LIMO Team 2014

% Andrew X Stewart, nov13
% Cyril Pernet, fixed the last dim to be arbitrary + size check, Jan 2014
% add the fored_dim argument -- making the function more generic

current = pwd;
if ~exist('LIMO','var')
    try
        load LIMO
    catch NO_FILE
        try
            cd ..; load LIMO
        catch NO_FILE
            if nargin == 1
                error('no LIMO variable/file found - use forced dim as 2nd argument')
            end
        end
    end
end

% if forced dim, override LIMO.data.size
if nargin == 2
    if numel(forced_dim) == 3
        LIMO.data.size3D = forced_dim;
    elseif numel(forced_dim) == 4
        LIMO.data.size4D = forced_dim;
    else
        error('forced dim dimension issue')
    end
end
cd(current)

% Check the size of input
reshape_size = size(reshape_in);

%%  If we have 4D input, reshape it to be 3D
if numel(reshape_size) == 4
    
    [n_elec, n_freqs, n_times, N] = size(reshape_in);
    n_freq_times = LIMO.data.size3D(2);
    if n_freq_times ~= n_freqs*n_times
        error('dimensions disagreement to reshape freq*time')
    else
        reshaped=nan(n_elec,n_freq_times,N);
    end
    
    %% For each trial, build a concatenated 2D matrix of (elec x (freqs timepoint)), and then populate the 3D Y with this
    for tr = 1:N
        for tm = 1:n_times
            if tm==1
                tf_long_2d = squeeze(reshape_in(:,:,tm,tr));
            else
                tf_long_2d = [tf_long_2d squeeze(reshape_in(:,:,tm,tr))];
            end
        end
        reshaped(:,:,tr) = tf_long_2d;
    end
    
elseif numel(reshape_size) == 3 || 2
    %% Else, if we have 3D input, reshape it to be 4D
    
    [n_elec, n_freq_times, N] = size(reshape_in);
    n_freqs = LIMO.data.size4D(2);
    n_times = LIMO.data.size4D(3);
    if n_freq_times ~= n_freqs*n_times
        error('dimensions disagreement to reshape freq*time')
    else
        reshaped = nan(n_elec,n_freqs, n_times, N);
    end
    
    % For each trial, take the stack of (elec x (freqs timepoint)) and split it into frames of (elec x freqs x times), then populate the 4D Y with this
    
    for tr = 1:N
        for tm = 1:n_times
            this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 2D tf 
                eft_3d(:,:,tm) =reshape_in(:,this_freq_start_index:this_freq_start_index+n_freqs-1,tr);
        end
        reshaped(:,:,:,tr) = eft_3d;
    end
end
    
    
    
    
    
    
                
                
                
                
                
    
    
    
    
    

function reshaped = limo_tf_4d_reshape(reshape_in)

% Simple function to 'stack' loaded elec-freqs-times-trials 4D
% time-frequency (ERSP) data into a 3D matrix
%
% This is done explicitly within simple loops to be clear to read.
% Could be vectorised and/or use reshape(), but that makes it easier to get
% lost in dimensions.
%
% Andrew X Stewart, nov13
% -----------------------------------------------------
%  Copyright (C) LIMO Team 2013

global LIMO
n_freqs = LIMO.data.size4D(2);
n_times = LIMO.data.size4D(3);


% Check the size of input
reshape_size = size(reshape_in);

%%  If we have 4D input, reshape it to be 3D
if numel(reshape_size) == 4
    
    disp('Reshaping 4D tf input...');
    
    [n_elec, n_freqs, n_times, n_trials] = size(reshape_in);
    reshaped=nan(n_elec,(n_freqs*n_times),n_trials);
    
    %% For each trial, build a concatenated 2D matrix of (elec x (freqs timepoint)), and then populate the 3D Y with this
    for tr = 1:n_trials
        
        for tm = 1:n_times
            
            if tm==1
                tf_long_2d = squeeze(reshape_in(:,:,tm,tr));
            else
                tf_long_2d = [tf_long_2d squeeze(reshape_in(:,:,tm,tr))];
            end
        end
        
        reshaped(:,:,tr) = tf_long_2d;
        
    end
    
elseif numel(reshape_size) == 3
    %% Else, if we have 3D input, reshape it to be 4D
    
     disp('Reshaping 3D tf input...');
    
    
    [n_elec, n_freq_times, n_trials] = size(reshape_in);
    reshaped = nan(LIMO.data.size4D);
   
    
    % For each trial, take the stack of (elec x (freqs timepoint)) and split it into frames of (elec x freqs x times), then populate the 4D Y with this
    
    for tr = 1:n_trials
        
        for tm = 1:n_times
            
            this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 2D tf 
            
 
                eft_3d(:,:,tm) =reshape_in(:,this_freq_start_index:this_freq_start_index+n_freqs-1,tr);

      
        end
        reshaped(:,:,:,tr) = eft_3d;
    end
end
    
    
    
    
    
    
                
                
                
                
                
    
    
    
    
    

function reshaped = limo_tf_5d_reshape(reshape_in,varargin)

% Simple function to 'unstack' loaded elec-freqs*times-N-boot 4D
% data into a 5D matrix elec-freqs-times-N-boot
%
% INPUT    4D electrode*[freq*time]*N*boot
% OUTPUT   5D electrode*freq*time*N*boot
%
% Copy of limo_tf_4d_reshape with a bootstrap loop
%
% Cyril Pernet, January 2014
% ------------------------------
%  Copyright (C) LIMO Team 2019

if nargin == 2
    LIMO = varargin{1};
    clear varargin
else
    if ~exist('LIMO','var')
        if exist(fullfile(pwd,'LIMO.mat'),'file')
            LIMO = load(fullfile(pwd,'LIMO.mat')); LIMO = LIMO.LIMO;
        elseif exist(fullfile(fileparts(pwd),'LIMO.mat'),'file')
            LIMO = load(fullfile(fileparts(pwd),'LIMO.mat')); LIMO = LIMO.LIMO;
        else
            global LIMO
            if ~isempty(LIMO)
                warning('No LIMO file file to unstack 5D data accessing global');
            else
                error('no LIMO.mat file found to unstack 5D file')
            end
        end
    end
end

n_freqs = LIMO.data.size4D(2);
n_times = LIMO.data.size4D(3);
[n_elec, n_freqs_times, N, n_boot] = size(reshape_in);
reshaped=nan(n_elec,n_freqs,n_times,N,n_boot);

for b = 1:n_boot
    for tr = 1:N
        for tm = 1:n_times
            this_freq_start_index = tm*n_freqs - n_freqs + 1;  % Set index in the long 3D tf
            eft_3d(:,:,tm) =reshape_in(:,this_freq_start_index:this_freq_start_index+n_freqs-1,tr,b);
        end
        reshaped(:,:,:,tr,b) = eft_3d;
    end
end
    
    
    
    
    
    
                
                
                
                
                
    
    
    
    
    

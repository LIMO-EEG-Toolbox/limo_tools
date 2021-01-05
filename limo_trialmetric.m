function [DT,TP,AC,DA] = limo_trialmetric(data,varargin)

% this function allows to compute 4 trial related metrics
% 1 - amplitude variations over time = std(trial) ; how much change over time
% 2 - total power = sum(abs(trial)^2)/length(trial) ; how much spectral energy [Parseval theorem] 
% 3 - autocorrelation = size of the peak of the autocorrelation function ; smoothness of the signal over time
% 4 - amplitude variations relative to a reference ~std(trial/ref); how much change relative to others
%
% FORMAT [DT,TP,AC,DA] = limo_trialmetric(data)
%        [DT,TP,AC,DA] = limo_trialmetric(data,options)
%        
% INPUT data is a matrix in 3D [channels*frames*trials] or 2D [frames*trials]
%       options are 'key' and 'value' pairs 
%                   'std_time':           'on' (default) or 'off'
%                   'delta_amplitude':    'on' or 'off' (default)
%                   'reference':          reference curve for delta amplitude
%                   'power':              'on' (default) or 'off'
%                   'autocorrelation':    'on' (default) or 'off'
%                   'sampling_frequency': used to smooth and express autocorrelation in ms
%
% OUTPUT [DT,TP,AC,DA] are the four computed metrics
%                      these are vectors, 1 value per trial
%
% Cyril Pernet 05-01-2021
% ----------------------------
% Copyright (C) LIMO Team 2021

%% deal with inputs
if nargin == 0
    help limo_trialmetric
    return
end

% data in
if ~isnumeric(data)
    error('1st argument must be the data matrix')
else
    if numel(size(data)) == 3
        [e,p,n]=size(data);
        Data = data; clear data
    else
        e = 1; % only one channel
        [p,n]=size(data);
        Data(1,:,:) = data; clear data
    end
end

if n == 1
    warning('Possible data size error, only 1 trial found');
end

% options
sampling_freq = 1;
op = struct('std_time','on','delta_amplitude','off','power','on','autocorrelation','on');
for o=1:length(varargin)
    if ischar(varargin{o})
        if strcmpi(varargin{o},'std_time') || strcmpi(varargin{o},'std time')
            op.std_time = varargin{o+1};
        elseif strcmpi(varargin{o},'delta_amplitude') || strcmpi(varargin{o},'delta amplitude')
            op.delta_amplitude = varargin{o+1};
        elseif strcmpi(varargin{o},'reference')
            if isnumeric(varargin{o+1})
                ref_amplitude = varargin{o+1};
                if all(size(ref_amplitude)==[p e])
                    ref_amplitude = ref_amplitude';
                elseif all(size(ref_amplitude,2)~=[e p])
                    warning('reference curve is not commensurate to the data, ''delta_amplitude'' is turned off')
                    op.delta_amplitude = 'off';
                end
            else
                warning('reference curve is not a vector as expected, ''delta_amplitude'' is turned off')
                op.delta_amplitude = 'off';
            end
        elseif strcmpi(varargin{o},'power')
            op.power = varargin{o+1};
        elseif contains(varargin{o},'autocor','IgnoreCase',true)
            op.autocorrelation = varargin{o+1};
        elseif contains(varargin{o},'sampling','IgnoreCase',true)
            if isnumeric(varargin{o+1})
                sampling_freq = varargin{o+1};
            end
        end
    end
end

%% compute
if strcmpi(op.std_time,'on')
    DT = squeeze(std(Data,0,2));
end

if strcmpi(op.power,'on')
    TP = squeeze(sum(abs(Data.^2),2)/p);       
end

if strcmpi(op.delta_amplitude ,'on')
    DA = squeeze(mean(sqrt(((Data-repmat(ref_amplitude,1,1,n)).^2)./p),2));
end

if strcmpi(op.autocorrelation,'on')
    AC = NaN(e,n);
    nfft = 2^nextpow2(2*p-1);
    parfor channel=1:e
        R = ifft(fft(detrend(squeeze(Data(channel,:,:))')',nfft).*conj(fft(detrend(squeeze(Data(channel,:,:))')',nfft)));
        R = R(1:p,:)'; % [R(end-p+2:end,:) ; R(1:p,:)]'; for a full window
        if sampling_freq ~=1
            R = medfilt1(R',1/sampling_freq*1000)'; % smooth to avoid edge effect
        end
        ac = NaN(n,1);
        for v=1:n
            if sum(R(v,:)) == 0
                ac(v) = 0;
            else
                [~,locs]=findpeaks(R(v,:));
                if ~isempty(locs)
                    [~,maxloc]=max(R(v,:));
                    tmp = locs-maxloc;
                    tmp(tmp <= 0) = NaN; % avoid being on itself
                    [~,closestpeak]=min(tmp);
                    closestpeak = locs(closestpeak);
                    [~,minloc]=min(R(v,maxloc:closestpeak));
                    ac(v) = maxloc+minloc;
                else
                    ac(v) = 0;
                end
            end
        end
        if sampling_freq ~=1
            AC(channel,:) = ac.*(1/sampling_freq*1000);
        else
            AC(channel,:) = ac;
        end
    end
end



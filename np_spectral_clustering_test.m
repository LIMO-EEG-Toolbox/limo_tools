function np_spectral_clustering_test(type)

% type can be H0 test that FWE is correct (validation)
% type can be H1 test that the computed stuff is correct

if strcmpi(type,'H0')

    for MC = 1:1000
        fprintf('running MC %g/1000\n',MC) % make null data
        fundamental_frequency = 5;         % Fundamental frequency in Hz
        num_harmonics         = 2;         % Number of harmonics
        duration              = 1;         % Duration of the time series in seconds
        sampling_rate         = 250;       % Sampling rate in Hz
        t = 0:1/sampling_rate:duration;    % Time vector
        signal = randn(20,size(t,2));      % Initialize the signal
        signal = signal + sin(2*pi*fundamental_frequency*t); % Add the fundamental frequency component
        for harmonic = 1:num_harmonics     % Add harmonics
            harmonic_frequency = fundamental_frequency * (harmonic + 1);     % Harmonic frequencies
            signal = signal + sin(2*pi*harmonic_frequency*t)./(1+harmonic);  % Dampen by half
        end
        N = size(signal,2); % Perform FFT
        for n=20:-1:1
            fft_result          = fft(signal(n,:), N);
            power_spectrum(n,:) = abs(fft_result).^2 / N;
        end
        [pvalues1(MC,:), cluster_pvalues1(MC,:)] = np_spectral_clustering(power_spectrum(1:10,1:length(power_spectrum)/2), ...
            power_spectrum(11:20,1:length(power_spectrum)/2),'type','paired','figure','off','verbose','off');
        [pvalues2(MC,:), cluster_pvalues2(MC,:)] = np_spectral_clustering(power_spectrum(1:10,1:length(power_spectrum)/2), ...
            power_spectrum(11:20,1:length(power_spectrum)/2),'type','independent','figure','off','verbose','off');
    end

    % compute the type 1 error rate
    error_singrank =  mean(nansum(cluster_pvalues1,2)~=0);
    error_ranksum  = mean(nansum(cluster_pvalues2,2)~=0);
    fprintf('matlab tests returns an error rate of %g for singrank and %g for ranksum\n', ...
       mean(pvalues1(:)<0.05), mean(pvalues2(:)<0.05))
     fprintf('the FWER error rate is of %g for singrank and %g for ranksum\n', ...
       mean(sum(cluster_pvalues1<0.05,2)), mean(sum(cluster_pvalues2<0.05,2)))
 

elseif strcmpi(type,'H1')
% ---------------------

        signal = randn(20,250); % a signal 
        signal(1:10,50:100)  = signal(1:10,50:100)+5; % +5 for 50 frames
        signal(1:10,150:250) = signal(1:10,150:250)+2.5; % +2.5 for 100 frames
        [~, cluster_pvalues] = np_spectral_clustering(signal(1:10,:), ...
            signal(11:20,:),'type','paired','figure','new','verbose','off');
        corrected_values = unique(cluster_pvalues(:));
        corrected_values(isnan(corrected_values)) = [];
        if length(corrected_values) == 1
            fprintf('signrank same cluster values because the mass is the same\n')
        else
            fprintf('signrank cluster value error because the mass is the same but p values are not\n')
        end

        subplot(1,2,2);
        [~, cluster_pvalues] = np_spectral_clustering(signal(1:10,:), ...
            signal(11:20,:),'type','independent','figure','new','verbose','off');
        corrected_values = unique(cluster_pvalues(:));
        corrected_values(isnan(corrected_values)) = [];
        if length(corrected_values) == 1
            fprintf('ranksum same cluster values because the mass is the same\n')
        else
            fprintf('ranksum cluster value error because the mass is the same but p values are not\n')
        end
end
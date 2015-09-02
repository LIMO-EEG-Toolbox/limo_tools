function [first_frame,last_frame,subj_chanlocs,limo] = match_frames(Paths,limo)


% once we have a series of files to assemble, we need to collect information 
% to match the frames across subjects
%
% FORMAT: [first_frame,last_frame,subj_chanlocs,channeighbstructmat] = limo_match_frames(Paths,limo)
%
% INPUT: Paths the full paths of files to assenble ; use to read the LIMO.mat
%        limo a LIMO.mat structure to be updated 
%
% OUTPUTS: first_frame last _frame returns the beginning and end in terms of
%          indices ; for time-frequency these are vectors with time st then frequency
%          subj_chanlocs is a cell of all subjects channel locations
%          channeighbstructmat is cell of all subject neighbouringhood matrices (if any)
% if the structure limo is set as global, then it is also updated the reflect the
% smallest interval(s) across subjects, which is used for the second leve analysis
%
% Cyril Pernet v2 August 2015
% ---------------------------------------------------------
%  Copyright (C) LIMO Team 2015

current = pwd; ME = [];

disp('match frames between subjects ...')
% check Paths format
if iscell(Paths{1})
    tmp = Paths; clear Paths
    index = 1;
    for gp=1:size(tmp,2)
        for s=1:size(tmp{gp},2)
            Paths{index} = tmp{gp}(s);
            index = index + 1;
        end
    end
end

% now loop loading the LIMO.mat for each subject to collect information
for i=1:size(Paths,2)
    try
        cd (Paths{i});
    catch
        cd (cell2mat(Paths{i}))
    end
    load LIMO;
    
    if i==1
        Analysis = LIMO.Analysis;
    else
        if ~strcmp(LIMO.Analysis,Analysis)
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
        end
    end
    
    sampling_rate(i)          = LIMO.data.sampling_rate;
    subj_chanlocs(i).chanlocs = LIMO.data.chanlocs;
    try
       channeighbstructmat = LIMO.data.channeighbstructmat;
    catch ME
    end
       
    if strcmp(Analysis,'Time-Frequency')
        first_frame(i,1)            = LIMO.data.trim1;
        last_frame(i,1)             = LIMO.data.trim2;
        start(i,1)                  = LIMO.data.start;
        stop(i,1)                   = LIMO.data.end;
        
        first_frame(i,2)            = LIMO.data.trim_low_f;
        last_frame(i,2)             = LIMO.data.trim_high_f;
        start(i,2)                  = LIMO.data.tf_freqs(1);
        stop(i,2)                   = LIMO.data.tf_freqs(end);
        
        tf_times{i}(1,:)            = LIMO.data.tf_times;
        tf_freqs{i}(1,:)            = LIMO.data.tf_freqs;
    else
        first_frame(i)              = LIMO.data.trim1;
        last_frame(i)               = LIMO.data.trim2;
        start(i)                    = LIMO.data.start;
        stop(i)                     = LIMO.data.end;
        
        if strcmp(Analysis,'Frequency')
            freqlist{i}(1,:)        = LIMO.data.freqlist;
        end
    end
end

% quick check things are ok
if ~isempty(ME) && isempty(limo.data.neighbouring_matrix)
    error(sprintf('some subject(s) have a different channel structure \nplease load an expected chanloc when choosing a test'));                          
end

if (sum(sampling_rate == sampling_rate(1))) ~= length(sampling_rate)
    error('data have different sampling rates')
end

% match and return into limo - the temp structure passed as global
limo.Analysis = Analysis;
limo.data.sampling_rate = sampling_rate(1);

% we need 1) to find the highest start in time and freq 2) the lowest end
% in time and freq and 3) match that on freqlist or tf_times/tf_freqs

[v,c] = max(first_frame);
if strcmp(Analysis,'Time-Frequency')
    limo.data.trim1 = v(1);
    limo.data.start = start(c(1),1);
    limo.data.trim_low_f = v(2);
    limo.data.low_f = start(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-limo.data.start));
        tf_times{i} = tf_times{i}(ind:end);
        [~,ind] = min(abs(tf_freqs{i}-limo.data.low_f));
        tf_freqs{i} = tf_freqs{i}(ind:end);
    end
else
    limo.data.trim1 = v;
    limo.data.start = start(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-limo.data.start));
            freqlist{i} = freqlist{i}(ind:end);
        end
    end
end

[v,c] = min(last_frame);
if strcmp(Analysis,'Time-Frequency')
    limo.data.trim2 = v(1);
    limo.data.end = stop(c(1),1);
    limo.data.trim_high_f = v(2);    
    limo.data.high_f = stop(c(2),2);
    % loop to adjust each subject list
    for i=1:size(Paths,2)
        [~,ind] = min(abs(tf_times{i}-limo.data.end));
        tf_times{i} = tf_times{i}(1:ind);
        [~,ind] = min(abs(tf_freqs{i}-limo.data.high_f));
        tf_freqs{i} = tf_freqs{i}(1:ind);
    end
else
    limo.data.trim2 = v;
    limo.data.end = stop(c);
    if strcmp(Analysis,'Frequency')
        % loop to adjust each subject list
        for i=1:size(Paths,2)
            [~,ind] = min(abs(freqlist{i}-limo.data.end));
            freqlist{i} = freqlist{i}(1:ind);
        end
    end
end


% finally match everything
if strcmp(Analysis,'Frequency')
    % check all lists match ; if sampled the same across subject, all same sizes
    try
        freqlist = cell2mat(freqlist');
    catch list_issue
        assignin('base','freqlist',freqlist)
        error('the resolution of frequency lists doesn''t match between subjects, check your workspace for details')
    end
    limo.data.freqlist = mean(freqlist,1);
    limo.data.start    = limo.data.freqlist(1);
    limo.data.end      = limo.data.freqlist(end);
    
elseif strcmp(Analysis,'Time-Frequency')
    % check all lists match
    try
        tf_times = cell2mat(tf_times');
        tf_freqs = cell2mat(tf_freqs');
    catch list_issue
        error('the resolution of time/frequency lists doesn''t match between subjects, check your workspace for details')
    end
    limo.data.tf_times = mean(tf_times,1);
    limo.data.tf_freqs = mean(tf_freqs,1);
    limo.data.start    = limo.data.tf_times(1);
    limo.data.low_f    = limo.data.tf_freqs(1);
    limo.data.end      = limo.data.tf_times(end);
    limo.data.high_f   = limo.data.tf_freqs(end);
end
cd(current)




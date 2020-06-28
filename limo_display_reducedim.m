function [dataout,channel,freq,time] = limo_display_reducedim(datain,LIMO)

% routine to reduce dimentionality of datain using user prompt
% returning dataout to plot
%
% FORMAT dataout = limo_display_reducedim(datain,LIMO)
% 
% INPUT datain 3D/4D data (usually stat file) or filename of data 
%       LIMO LIMO structure or filename of the LIMO file to use
%
% OUTPUT dataout are the data to plot
%        channel is the channel/component number selected
%        time_or_freq is the time frame or frequency selected
%
% Cyril Pernet 
% ------------------------------
%  Copyright (C) LIMO Team 2020

%% Input checks
if nargin == 0
    help limo_display_reducedim
    return
elseif nargin == 1 || nargin > 2
    help limo_display_reducedim
    error('Wrong number of arguments in')
end

if ~isnumeric(datain)
    datain = load(datain);
    datain = datain.(cell2mat(fieldnames(datain)));
end

if ~isstruct(LIMO)
    LIMO = load(LIMO);
    LIMO = LIMO.LIMO;
end    

if ~isfield(LIMO, 'Analysis')
    error('Analysis field missing from LIMO file, can''t workout dimensions to reduce')
end

%% ask user what to do
channel      = [];
if strcmpi(LIMO.Analysis,'Time-Frequency')
    freq   = 1:size(datain,2);
    time   = 1:size(datain,3);
else
    freq   = 1:size(datain,2);
    time   = 1:size(datain,2);
end

% for ERSP choose time or frequency
if strcmpi(LIMO.Analysis,'Time-Frequency')
    restrict = questdlg('Which domain to plot:','ERSP plot option','Time','Frequency','Time');
    
    if strcmpi(restrict,'Frequency')
        time = cell2mat(inputdlg('At which time to plot the spectrum?','Plotting option'));
        if isempty(time)|| strcmp(time,'0')
            v              = max(datain(:,:,:,1),[],2);
            [channel,time] = ind2sub(size(v),find(v == max(v(:))));
        else
            time = eval(time);
        end
        [~,time] = min(abs(LIMO.data.tf_times - time));
        datain = squeeze(datain(:,:,time,:));
    elseif strcmpi(restrict,'Time')
        freq =  cell2mat(inputdlg('At which frequency to plot the time course?','Plotting option'));
        if isempty(freq)|| strcmp(freq,'0')
            v              = max(datain(:,:,:,1),[],3);
            [channel,freq] = ind2sub(size(v),find(v == max(v(:))));
        else
            freq = eval(freq);
        end
        [~,freq] = min(abs(LIMO.data.tf_freqs - freq));
        datain = squeeze(datain(:,freq,:,:));
    end
end

% reduce space to a single channel
% data alsways 3D because space or freq squished above if ERSP
if size(datain,1) == 1
    channel = 1;
elseif isempty(channel)
    if strcmpi(LIMO.Type,'Channels')
        channel = inputdlg('which channel to plot','Plotting option');
    else
        channel = inputdlg('which component to plot','Plotting option');
    end
    
    if isempty(channel) || strcmp(cell2mat(channel),'')
        [v,e]   = max(datain(:,:,1));
        [~,c]   = max(v);
        channel = e(c);
    else
        channel = eval(cell2mat(channel));
        if length(channel) > 1
            error('1 %s only can be plotted',LIMO.Type(1:end-1));
        elseif channel > size(datain,1)
            error('%s number invalid',LIMO.Type(1:end-1));
        end
    end
end
dataout = squeeze(datain(channel,:,:));


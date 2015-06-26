function limo_itc_import_data(setfile,cat,cont,defaults)
global EEGLIMO

EEGLIMO=pop_loadset(setfile);
[root,name,ext] = fileparts(setfile); 
cd(root); mkdir(['ITC_analysis']);
LIMO.dir                    = [pwd filesep 'ITC_analysis'];
LIMO.data.data              = [name ext];
LIMO.data.data_dir          = root;
LIMO.data.sampling_rate     = EEGLIMO.srate;


if isfield(defaults,'chanlocs')
    LIMO.data.chanlocs      = defaults.chanlocs;
else
    LIMO.data.chanlocs      = EEGLIMO.chanlocs;
end

if strcmp(defaults.analysis,'Time') 
    
    % start
    if isempty(defaults.start)
        LIMO.data.start = min(EEGLIMO.times);
        LIMO.data.trim1 = 1;    
    elseif defaults.start < min(EEGLIMO.times)
        error(['The earliest time possible is:',num2str(EEGLIMO.times(1)),'ms']);
        return
    else
        [~,position]=min(abs(EEGLIMO.times - defaults.start));
        LIMO.data.start = EEGLIMO.times(position);
        LIMO.data.trim1 = position;
    end
    
    % end
    if isempty(defaults.end)
        LIMO.data.start = max(EEGLIMO.times);
        LIMO.data.trim2 = length(EEGLIMO.times);    
    elseif defaults.end > max(EEGLIMO.times)
        error(['The latest time possible is:',num2str(EEGLIMO.times(end)),'ms']);
        return
    else
        [~,position]=min(abs(EEGLIMO.times - defaults.end));
        LIMO.data.end = EEGLIMO.times(position);
        LIMO.data.trim2 = position;
    end
    
elseif strcmp(defaults.analysis,'Frequency') 
    
    disp('Frequency import for ITC data')
    
    % start
    if isempty(defaults.lowf)
        LIMO.data.start = EEGLIMO.etc.tf_freqs(1);
        LIMO.data.trim1 = 1;
    elseif defaults.lowf < EEGLIMO.etc.tf_freqs(1)
        error(['The lowest frequency possible is:',num2str(EEGLIMO.etc.tf_freqs(1))]);
        return
    else
        [~,position] = min(abs(EEGLIMO.etc.tf_freqs-defaults.lowf));
        LIMO.data.start = EEGLIMO.etc.tf_freqs(position);
        LIMO.data.trim1 = position; 
    end
    
    % end
    if isempty(defaults.highf)
        LIMO.data.end = EEGLIMO.etc.tf_freqs(end);
        LIMO.data.trim2 = length(EEGLIMO.etc.tf_freqs);
    elseif defaults.highf > EEGLIMO.etc.tf_freqs(end)
        error(['The highest frequency possible is:',num2str(EEGLIMO.etc.tf_freqs(end))]);
        return
    else
        [~,position] = min(abs(EEGLIMO.etc.tf_freqs-defaults.highf));
        LIMO.data.end = EEGLIMO.etc.tf_freqs(position);
        LIMO.data.trim2 = position; 
    end
    
    LIMO.data.freqlist = EEGLIMO.etc.tf_freqs(LIMO.data.trim1:LIMO.data.trim2);


elseif strcmp(defaults.analysis,'Time-Frequency') || strcmp(defaults.analysis,'ITC')
    
    LIMO.data.chanlocs      = EEGLIMO.chanlocs;
    
    LIMO.data.tf_data_filepath = EEGLIMO.etc.tf_path;
    if ~exist(EEGLIMO.etc.tf_path,'file')
        cd(LIMO.data.data_dir)
        [p,f,e]=fileparts(EEGLIMO.etc.tf_path);
        if exist([pwd filesep f e])
            LIMO.data.tf_data_filepath = [pwd filesep f e];
        end
    end
    
    % start
    if isempty(defaults.start)
        LIMO.data.start = min(EEGLIMO.etc.tf_times);
        LIMO.data.trim1 = 1;    
    elseif defaults.start < min(EEGLIMO.etc.tf_times)
        error(['The earliest time possible is:',num2str(EEGLIMO.etc.tf_times(1)),'ms']);
        return
    else
        [~,position]=min(abs(EEGLIMO.etc.tf_times - defaults.start));
        LIMO.data.start = EEGLIMO.etc.tf_times(position);
        LIMO.data.trim1 =  find(EEGLIMO.etc.tf_times == LIMO.data.start);
    end
    
    % end
    if isempty(defaults.end)
        LIMO.data.end = max(EEGLIMO.etc.tf_times);
        LIMO.data.trim2 = length(EEGLIMO.etc.tf_times);    
    elseif defaults.end > max(EEGLIMO.etc.tf_times)
        error(['The latest time possible is:',num2str(EEGLIMO.etc.tf_times(end)),'ms']);
        return
    else
        [~,position]=min(abs(EEGLIMO.etc.tf_times - defaults.end));
        LIMO.data.end = EEGLIMO.etc.tf_times(position);
        LIMO.data.trim2 =  find(EEGLIMO.etc.tf_times == LIMO.data.end);
    end

    LIMO.data.tf_times = EEGLIMO.etc.tf_times(LIMO.data.trim1:LIMO.data.trim2);

    % start
    if isempty(defaults.lowf)
        LIMO.data.lowf = EEGLIMO.etc.tf_freqs(1);
        LIMO.data.trim_low_f = 1;
    elseif defaults.lowf < EEGLIMO.etc.tf_freqs(1)
        error(['The lowest frequency possible is:',num2str(EEGLIMO.etc.tf_freqs(1))]);
        return
    else
        [~,position] = min(abs(EEGLIMO.etc.tf_freqs-defaults.lowf));
        LIMO.data.lowf = EEGLIMO.etc.limo_psd_freqlist(position);
        LIMO.data.trim_low_f = position; 
    end
    
    % end
    if isempty(defaults.highf)
        LIMO.data.highf = EEGLIMO.etc.tf_freqs(end);
        LIMO.data.trim_high_f = length(EEGLIMO.etc.tf_freqs);
    elseif defaults.highf > EEGLIMO.etc.tf_freqs(end)
        error(['The highest frequency possible is:',num2str(EEGLIMO.etc.tf_freqs(end))]);
        return
    else
        [~,position] = min(abs(EEGLIMO.etc.tf_freqs-defaults.highf));
        LIMO.data.hightf = EEGLIMO.etc.tf_freqs(position);
        LIMO.data.trim_high_f = position; 
    end
    
    LIMO.data.tf_freqs = EEGLIMO.etc.tf_freqs(LIMO.data.trim_low_f:LIMO.data.trim_high_f);
end

% if strcmp(cat(end-3:end),'.txt')
%     LIMO.data.Cat = load(cat);
% else strcmp(cat(end-3:end),'.mat')
%     load(cat); LIMO.data.Cat = eval(cat(1:end-4));
% end
% 
% if strcmp(cont(end-3:end),'.txt')
%     LIMO.data.Cont = load(cont);
% else strcmp(cont(end-3:end),'.mat')
%     load(cont); LIMO.data.Cont = eval(cont(1:end-4));
% end

LIMO.Analysis                = defaults.analysis;
LIMO.design.zscore           = defaults.zscore;
LIMO.design.method           = defaults.method;
LIMO.design.type_of_analysis = defaults.type_of_analysis;
LIMO.design.fullfactorial    = defaults.fullfactorial;
LIMO.design.bootstrap        = defaults.bootstrap;
LIMO.design.tfce             = defaults.tfce;
LIMO.Level                   = 1;
cd(LIMO.dir); save LIMO LIMO; cd ..

end

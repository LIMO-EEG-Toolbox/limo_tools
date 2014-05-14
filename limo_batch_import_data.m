function limo_batch_import_data(setfile,cat,cont,defaults)
global EEG

EEG=pop_loadset(setfile);
[root,name,ext] = fileparts(setfile); 
cd(root); mkdir(['GLM_' defaults.analysis]);
LIMO.dir                    = [pwd filesep 'GLM_' defaults.analysis];
LIMO.data.data              = [name ext];
LIMO.data.data_dir          = root;
LIMO.data.sampling_rate     = EEG.srate;
if isfield(defaults,'chanlocs')
    LIMO.data.chanlocs      = defaults.chanlocs;
else
    LIMO.data.chanlocs      = EEG.chanlocs;
end

if strcmp(defaults.analysis,'Time') 
    
    % start
    if isempty(defaults.start)
        LIMO.data.start = min(EEG.times);
        LIMO.data.trim1 = 1;    
    elseif defaults.start < min(EEG.times)
        error(['The earliest time possible is:',num2str(EEG.times(1)),'ms']);
        return
    else
        [~,position]=min(abs(EEG.times - defaults.start));
        LIMO.data.start = EEG.times(position);
        LIMO.data.trim1 = position;
    end
    
    % end
    if isempty(defaults.end)
        LIMO.data.start = max(EEG.times);
        LIMO.data.trim2 = length(EEG.times);    
    elseif defaults.end > max(EEG.times)
        error(['The latest time possible is:',num2str(EEG.times(end)),'ms']);
        return
    else
        [~,position]=min(abs(EEG.times - defaults.end));
        LIMO.data.end = EEG.times(position);
        LIMO.data.trim2 = position;
    end
    
elseif strcmp(defaults.analysis,'Frequency') 
    
    % start
    if isempty(defaults.lowf)
        LIMO.data.start = EEG.etc.limo_psd_freqlist(1);
        LIMO.data.trim1 = 1;
    elseif defaults.lowf < EEG.etc.limo_psd_freqlist(1)
        error(['The lowest frequency possible is:',num2str(EEG.etc.limo_psd_freqlist(1))]);
        return
    else
        [~,position] = min(abs(EEG.etc.limo_psd_freqlist-defaults.lowf));
        LIMO.data.start = EEG.etc.limo_psd_freqlist(position);
        LIMO.data.trim1 = position; 
    end
    
    % end
    if isempty(defaults.highf)
        LIMO.data.end = EEG.etc.limo_psd_freqlist(end);
        LIMO.data.trim2 = 1;
    elseif defaults.highf > EEG.etc.limo_psd_freqlist(end)
        error(['The highest frequency possible is:',num2str(EEG.etc.limo_psd_freqlist(end))]);
        return
    else
        [~,position] = min(abs(EEG.etc.limo_psd_freqlist-defaults.highf));
        LIMO.data.end = EEG.etc.limo_psd_freqlist(position);
        LIMO.data.trim2 = position; 
    end
    
    LIMO.data.freqlist = EEG.etc.limo_psd_freqlist(LIMO.data.trim1:LIMO.data.trim2);


elseif strcmp(defaults.analysis,'Time-Frequency')
    
    LIMO.data.tf_data_filepath = EEG.etc.tf_path;
    if ~exist(EEG.etc.tf_path,'file')
        cd(LIMO.data.data_dir)
        [p,f,e]=fileparts(EEG.etc.tf_path);
        if exist([pwd filesep f e])
            LIMO.data.tf_data_filepath = [pwd filesep f e];
        end
    end
    
    % start
    if isempty(defaults.start)
        LIMO.data.start = min(EEG.etc.tf_times);
        LIMO.data.trim1 = 1;    
    elseif defaults.start < min(EEG.etc.tf_times)
        error(['The earliest time possible is:',num2str(EEG.etc.tf_times(1)),'ms']);
        return
    else
        [~,position]=min(abs(EEG.etc.tf_times - defaults.start));
        LIMO.data.start = EEG.etc.tf_times(position);
        LIMO.data.trim1 =  find(EEG.etc.tf_times == LIMO.data.start);
    end
    
    % end
    if isempty(defaults.end)
        LIMO.data.end = max(EEG.etc.tf_times);
        LIMO.data.trim2 = length(EEG.etc.tf_times);    
    elseif defaults.end > max(EEG.etc.tf_times)
        error(['The latest time possible is:',num2str(EEG.etc.tf_times(end)),'ms']);
        return
    else
        [~,position]=min(abs(EEG.etc.tf_times - defaults.end));
        LIMO.data.end = EEG.etc.tf_times(position);
        LIMO.data.trim2 =  find(EEG.etc.tf_times == LIMO.data.end);
    end

    LIMO.data.tf_times = EEG.etc.tf_times(LIMO.data.trim1:LIMO.data.trim2);

    % start
    if isempty(defaults.lowf)
        LIMO.data.lowf = EEG.etc.tf_freqs(1);
        LIMO.data.trim_low_f = 1;
    elseif defaults.lowf < EEG.etc.tf_freqs(1)
        error(['The lowest frequency possible is:',num2str(EEG.etc.tf_freqs(1))]);
        return
    else
        [~,position] = min(abs(EEG.etc.tf_freqs-defaults.lowf));
        LIMO.data.lowf = EEG.etc.limo_psd_freqlist(position);
        LIMO.data.trim_low_f = position; 
    end
    
    % end
    if isempty(defaults.highf)
        LIMO.data.highf = EEG.etc.tf_freqs(end);
        LIMO.data.trim_high_f = 1;
    elseif defaults.highf > EEG.etc.tf_freqs(end)
        error(['The highest frequency possible is:',num2str(EEG.etc.tf_freqs(end))]);
        return
    else
        [~,position] = min(abs(EEG.etc.tf_freqs-defaults.highf));
        LIMO.data.hightf = EEG.etc.tf_freqs(position);
        LIMO.data.trim_high_f = position; 
    end
    
    LIMO.data.tf_freqs = EEG.etc.tf_freqs(LIMO.data.trim_low_f:LIMO.data.trim_high_f);
end

if strcmp(cat(end-3:end),'.txt')
    LIMO.data.Cat = load(cat);
else strcmp(cat(end-3:end),'.mat')
    load(cat); LIMO.data.Cat = eval(cat(1:end-4));
end

if strcmp(cont(end-3:end),'.txt')
    LIMO.data.Cont = load(cont);
else strcmp(cont(end-3:end),'.mat')
    load(cont); LIMO.data.Cont = eval(cont(1:end-4));
end

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

function limo_checkweight(LIMO_files, expected_chanlocs, varargin)

% general function designed to look at the weights computed for each trial
% at the each channel for each subject
%
% FORMAT limo_CheckWeight(list_of_LIMO.mat,expected_chanlocs)
%        limo_CheckWeight(list_of_LIMO.mat,expected_chanlocs,'CheckBias','on','TestDifference','on','SingleSubjectAnalysis','off','PlotRank','off')
%
% INPUT list_of_LIMO.mat empty [] calls a gui to select a txt file
%                        txt file listing where the LIMO.mat are located
%                        cell array listing where the LIMO.mat are located
%
%       expected_chanlocs the default channel locations for all subject
%
%       options 'PlotRank' is 'on' (default) or 'off'
%                since weights are between 0 and 1, it computes the average response for each decile
%               'TestDifference' is 'on' (default) or 'off' to compute an OLS between the good trials (weights = 1/0.9)
%                                 and the outliers (by reverse engineering the weights to outlier detection)
%               'CheckBias' is 'on' (default) or 'off' to check that the weights are distributed across trials in a uniform manner, 
%                            i.e. that not one conditions is more affected than another which would bias the results, but
%                            also indicate that something is going on in the data
%               'SingleSubjectAnalysis' is 'on' or 'off' (default) to perform ttests subject wise 
%                                         (difference may exist but not consistent across subjects)
%
% OUTPUT creates a folder called 'Weights_checking' with the different
%        results in it
%
% Cyril Pernet January 2021
% -----------------------------
%  Copyright (C) LIMO Team 2021

%% if no input do it all
limo = struct('PlotRank','on','TestDifference','on','CheckBias','on','SingleSubjectAnalysis','off');

if nargin == 0
    [~,~,LIMO_files] = limo_get_files([],'*txt','choose a list of LIMO files');
    if isempty(LIMO_files)
        return
    end
    [to_load,path] = uigetfile2('expected_chanlocs.mat','load chanlocs'); 
    if to_load == 0
        return
    end
    chan                          = load([path to_load]);
    limo.data.chanlocs            = chan.expected_chanlocs;
    limo.data.neighbouring_matrix = chan.channeighbstructmat;
end

%% input checks
if isempty(LIMO_files)
    [~,~,LIMO_files] = limo_get_files([],'*txt','choose a list of LIMO files');
elseif ischar(LIMO_files)
    files = textread(LIMO_files,'%s','delimiter','');  % select a txt file listing all files
    clear LIMO_files
    for f=1:size(files,1)
        LIMO_files{f} = files{f};
    end
end

% those files are there?
for f = length(LIMO_files):-1:1
    if ~exist(LIMO_files{f},'file')
        error([LIMO_files{f} ' doesn''t exist'])
    end
    
    [limo_paths{f},name,ext]=fileparts(LIMO_files{f});
    if ~strcmp([name ext],'LIMO.mat')
        error([LIMO_files{f} ' is not a LIMO.mat file'])
    end
    
    LIMO_sub = load(LIMO_files{f});
    LIMO_sub = LIMO_sub.LIMO;
    if f==length(LIMO_files)
        limo.Analysis = LIMO_sub.Analysis;
        limo.Type     = LIMO_sub.Type;
    else
        if limo.Analysis ~= LIMO_sub.Analysis
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
        end
        
        if limo.Type ~= LIMO_sub.Type
            error('Looks like different type of analyses (Channels/Components) are mixed up')
        end    
    end
end

if ~isfield(limo,'data')
    if nargin >= 2
        chan                          = load(expected_chanlocs);
    else
        [to_load,path]                = uigetfile2('*mat','load chanlocs');
        chan                          = load([path to_load]);
    end
    try
        limo.data.chanlocs                = chan.expected_chanlocs;
        limo.data.neighbouring_matrix     = chan.channeighbstructmat;
    catch different_names
        warndlg('couldn''t read chanlocs: %s\n',different_names.message)
        FN = fieldnames(chan);
        for i=1:length(FN)
            if isstruct(chan.(FN{i}))
                limo.data.chanlocs = chan.(FN{i});
            elseif ismatrix(chan.(FN{i}))
                limo.data.neighbouring_matrix = chan.(FN{i});
            end
        end
    end
end

if ~isempty(varargin)
    for n=1:length(varargin)
        if strcmpi(varargin{n},'plotrank')
            limo.PlotRank = varargin{n+1};
        end
        
        if strcmpi(varargin{n},'testdifference')
            limo.TestDifference = varargin{n+1};
        end
        
        if strcmpi(varargin{n},'checkbias')
            limo.CheckBias = varargin{n+1};
        end
        
        if strcmpi(varargin{n},'singlesubjectanalysis')
            limo.SingleSubjectAnalysis = varargin{n+1};
        end
    end
end

% -----------
% Compute
% ----------
mkdir('Weights_checking'); cd('Weights_checking'); LIMO.dir = pwd;
[~,~,subj_chanlocs,limo] = limo_match_frames(limo_paths,limo);

%% Compute
if strcmpi(limo.Analysis,'Time-Frequency')
    data       = NaN(size(limo.data.chanlocs,2),(limo.data.highf-limo.data.lowf+1),(limo.data.trim2-limo.data.trim1+1),length(LIMO_files),10);
    difference = NaN(size(data,1),size(data,2),size(data,3),size(data,4));
else
    data       = NaN(size(limo.data.chanlocs,2),(limo.data.trim2-limo.data.trim1+1),length(LIMO_files),10);
    difference = NaN(size(data,1),size(data,2),size(data,3));
end

for f=1:length(LIMO_files)
    fprintf('processing data subject %g\n',f)
    LIMO_sub = load(LIMO_files{f});
    LIMO_sub = LIMO_sub.LIMO;
    W{f}     = LIMO_sub.design.weights;
    Yr       = load([LIMO_sub.dir filesep 'Yr.mat']);
    Yr       = Yr.Yr;
    
    %% make averages
    if strcmp(limo.PlotRank,'on')
        if strcmpi(limo.Analysis,'Time-Frequency')
            array = find(~isnan(Yr(:,1,1,1))); % skip empty electrodes
            tmp   = NaN(size(Yr,1),size(Yr,2),size(Yr,3),10);
            for e=1:length(array)
                for freq = 1:size(Yr,2)
                    if strcmpi(LIMO_sub.design.method,'WLS')
                        meanw = squeeze(W{f}(e,freq,:));
                    elseif strcmpi(LIMO_sub.design.method,'IRLS')
                        meanw = mean(squeeze(W{f}(e,freq,:,:)),1); % mean over time
                    end
                    % compute mean per quantile
                    Q = quantile(meanw,10);
                    for w=1:10
                        if w == 1
                            index = meanw <= Q(w);
                        elseif w == 10
                            index = meanw > Q(w);
                        else
                            index = logical((meanw>(Q(w-1))) .* (meanw<=(Q(w))));
                        end
                        tmp(e,:,w) = mean(Yr(e,:,index),3);
                    end
                    data(:,freq,:,f,:) = limo_match_elec(subj_chanlocs(f).chanlocs,chan.expected_chanlocs,1,size(Yr,2),tmp);
                end
            end
        else
            array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
            tmp = NaN(size(Yr,1),size(Yr,2),10);
            for e=1:length(array)
                if strcmpi(LIMO_sub.design.method,'WLS')
                    meanw = squeeze(W{f}(e,:));
                elseif strcmpi(LIMO_sub.design.method,'IRLS')
                    meanw = mean(squeeze(W{f}(e,:,:)),1); % mean ove time
                end
                % compute mean per quantile
                Q = quantile(meanw,10);
                for w=1:10
                    if w == 1
                        index = meanw <= Q(w);
                    elseif w == 10
                        index = meanw > Q(w);
                    else
                        index = logical((meanw>(Q(w-1))) .* (meanw<=(Q(w))));
                    end
                    tmp(e,:,w) = mean(Yr(e,:,index),3);
                end
            end
            data(:,:,f,:) = limo_match_elec(subj_chanlocs(f).chanlocs,chan.expected_chanlocs,1,size(Yr,2),tmp);
        end
        
        % save using the same format as limo_central_tendency, then call
        % limo_add_plots to make the figure - only parameters and subjedcts are
        % reversed because we want to plot per subjects and not per parameters
        if f==length(LIMO_files)
            clear Data; Data.data = data; Data.limo = limo;
            save subjects_weighted_data Data
        end
    end % close rank computation
    
    %% t-test good vs outliers
    if f == 1 && strcmpi(limo.SingleSubjectAnalysis,'on')
        LIMO.design.method = 'Mean';
        LIMO.Level         = 2;
        LIMO.Analysis      = limo.Analysis;
        LIMO.Type          = limo.Type;
    end

    % is there a difference between the outlier trials and the best trials
    if strcmp(limo.TestDifference,'on') || strcmpi(limo.SingleSubjectAnalysis,'on')
        if strcmpi(limo.Analysis,'Time-Frequency')
            [chan,freq,time,tr]=size(Yr);
            Yr = limo_tf_4d_reshape(Yr,[chan freq*time tr]);
        end
        
        tmp = NaN(size(Yr,1),size(Yr,2));
        for e=1:size(Yr,1)
            if strcmpi(LIMO_sub.design.method,'WLS')
                Q      = quantile(W{f}(e,:),10);
                index1 = find(W{f}(e,:) <=Q(1));
                index2 = find(W{f}(e,:) >=Q(10));
            elseif strcmpi(LIMO_sub.design.method,'IRLS')
                Q      = quantile(mean(squeeze(W{f}(e,:,:)),1),10);
                index1 = find(mean(squeeze(W{f}(e,:,:)),1) <=Q(1));
                index2 = find(mean(squeeze(W{f}(e,:,:)),1) >=Q(10));
            end
            
            for expected=1:size(limo.data.chanlocs,2)
                if strcmp(LIMO_sub.data.chanlocs(e).labels,limo.data.chanlocs(expected).labels)
                    difference(expected,:,f) = mean(squeeze(Yr(e,:,index1)),2) - mean(squeeze(Yr(e,:,index2)),2);
                    
                    % rather than mean, we can test differences per subjects on trials
                    if strcmpi(limo.SingleSubjectAnalysis,'on')
                        if f == 1 && e==1 % first subject
                            mkdir('single_subjects'); cd('single_subjects'); 
                        end
                        
                        if e==1 % first channel
                            if contains(LIMO_sub.dir,'sub-')
                                nb = LIMO_sub.dir(strfind(LIMO_sub.dir,'sub-')+4:end);
                                nb = nb(1:min(strfind(nb,filesep))-1);
                                mkdir(['sub-' nb]); cd(['sub-' nb])
                            else
                                mkdir(['sub-' num2str(f)]); cd(['sub-' num2str(f)])
                            end
                        end
                        
                        [two_samples(e,:,1),two_samples(e,:,3),~,sd,~,two_samples(e,:,4),two_samples(e,:,5)] = ...
                            limo_ttest(2, squeeze(Yr(e,:,index1)), squeeze(Yr(e,:,index2)), 0.05);
                        sd = sd.^2; a = sd(1,:)./length(index1); b = sd(1,:)./length(index2);
                        two_samples(e,:,2) = sqrt(a + b);                        
                    end
                end
            end
        end
                               
        if strcmpi(limo.SingleSubjectAnalysis,'on')
            if strcmp(LIMO_sub.Analysis,'Time-Frequency') ||  strcmp(LIMO_sub.Analysis,'ITC')
                two_samples = limo_tf_4d_reshape(two_samples);
            end
            save ('ttest_good_vs_outliers_trials','two_samples', '-v7.3')
            LIMO.dir = pwd; save(fullfile(LIMO.dir,'LIMO.mat'));
            Y1r = Yr(:,:,index1); save Y1r Y1r; clear Y1r
            Y2r = Yr(:,:,index2); save Y2r Y2r; clear Y2r
            % some trial metrics
            [DT1,TP1,AC1] = limo_trialmetric(Yr(:,:,index1),'std_time','on','power','on',...
                'autocorrelation','on','sampling_frequency',LIMO_sub.data.sampling_rate);
            [DT2,TP2,AC2] = limo_trialmetric(Yr(:,:,index2),'std_time','on','power','on',...
                'autocorrelation','on','sampling_frequency',LIMO_sub.data.sampling_rate);
                        
            % finally pick the channel with the most variance
            [~,maxchannels(f)]=max(mean(squeeze(two_samples(:,:,2)),2)); 
            N = max(length(index1),length(index2)); X = NaN(N,8);
            X(1:length(index1),1:4) = [index1' DT1(maxchannels(f),:)',TP1(maxchannels(f),:)',AC1(maxchannels(f),:)'];
            X(1:length(index2),5:8) = [index2' DT2(maxchannels(f),:)',TP2(maxchannels(f),:)',AC2(maxchannels(f),:)'];
            t = table(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6),X(:,7),X(:,8),'VariableNames',...
                {'low_weight_trials', 'time_var_lt', 'power_lt', 'autocorr_lt', ...
                'high_weight_trials','time_var_ht', 'power_ht', 'autocorr_ht'});
            writetable(t,['metrics_maxchannel' num2str(maxchannels(f)) '.csv']);
            SM(f,:) = nanmedian(X(:,[2 3 4 6 7 8]));
            SMT(f,:) = [median(DT1(:)) median(TP1(:)) median(AC1(:)) median(DT2(:)) median(TP2(:)) median(AC2(:))];
            cd .. ;
            
            % where are those channels
            if f==length(LIMO_files)
                t = table(SM(:,1),SM(:,2),SM(:,3),SM(:,4),SM(:,5),SM(:,6),'VariableNames',...
                    {'time_var_lt', 'power_lt', 'autocorr_lt', ...
                    'time_var_ht', 'power_ht', 'autocorr_ht'});
                writetable(t,'mean_metrics_maxchannels.csv');
                
                t = table(SMT(:,1),SMT(:,2),SMT(:,3),SMT(:,4),SMT(:,5),SMT(:,6),'VariableNames',...
                    {'time_var_lt', 'power_lt', 'autocorr_lt', ...
                    'time_var_ht', 'power_ht', 'autocorr_ht'});
                writetable(t,'mean_metrics_allchannels.csv');

                % locations
                data = zeros(1,length(limo.data.chanlocs));
                for s=1:length(LIMO_files)
                    data(maxchannels(s)) = data(maxchannels(s))+1;
                end
                all = unique(data); all(1) = [];
                save maxchannels maxchannels
                
                % create the frequency map
                h = figure('Color','w','NumberTitle','off','Name','limo_tools: best electrode frequency map');
                [~,grid_or_val]= topoplot( data,limo.data.chanlocs,'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on','conv','on');
                freqmap = grid_or_val(end:-1:1,:); % reverse row order to get back of the head at the bottom of the image
                freqmap = freqmap + abs(min(freqmap(:))); % min 0
                imagesc(freqmap,[0 max(freqmap(:))])
                axis tight;axis square;axis off
                colormap(limo_color_images(freqmap));
                title(sprintf('max variance channels from %g to %g subjects',min(all),max(all)));
                try saveas(gcf, 'max_var_channels.fig','fig'); end
                cd ..;
            end
        end
        
        if f==length(LIMO_files)
            mkdir('trial_differences'); cd('trial_differences'); save Yr difference; 
            
            % stats
            disp('Computing t-test between good and bad trials across all conditions')
            LIMO = limo; LIMO.dir = pwd; LIMO.Level = 2; LIMO.design.bootstrap = 1000;
            LIMO.design.electrode = []; LIMO.design.name = 'one sample ttest'; 
            LIMO.design.type_of_analysis = 'Mass-univariate'; LIMO.design.tfce = 0; 
            save LIMO LIMO; limo_random_robust(1,fullfile(LIMO.dir,'Yr.mat'),0,LIMO);
            movefile(fullfile(LIMO.dir,'one_sample_ttest_parameter_0.mat'),...
                fullfile(LIMO.dir,'one_sample_ttest_low_vs_high_weight_trials.mat'));
            movefile(fullfile(LIMO.dir,[filesep 'H0' filesep 'H0_one_sample_ttest_parameter_0.mat']),...
                fullfile(LIMO.dir,[filesep 'H0' filesep 'H0_one_sample_ttest_low_vs_high_weight_trials.mat']));
            cd ..; clear LIMO; 
        end
    end % close test difference
    
    %% Bias analysis
    if strcmp(limo.CheckBias,'on')
        if LIMO_sub.design.nb_conditions == 0
            sprintf('skipping bias test, there are no conditions, subject %g \n',f);
        else
            if LIMO_sub.design.nb_interactions ~= 0
                tmpX = LIMO_sub.design.X(:,sum(LIMO_sub.design.nb_conditions):(size(LIMO_sub.design.X,2)-LIMO_sub.design.nb_continuous));
            elseif LIMO_sub.design.nb_interactions == 0 && length(LIMO_sub.design.nb_conditions)>1
                [tmpX,~] = limo_make_interactions(LIMO_sub.design.X, LIMO_sub.design.nb_conditions);
            else
                tmpX = LIMO_sub.design.X(:,1:LIMO_sub.design.nb_conditions);
            end
            % tmpX is always a matrix with one column per condition, no
            % matter the design - we want to know if the weights are
            % distributed equally
            if strcmpi(limo.Analysis,'Time-Frequency')
                % LIMO.data.4d
                tmp = NaN(size(W{f},1),size(Yr,2),size(tmpX,2));
                for expected=1:size(chan.expected_chanlocs,2)
                    for e=1:size(W{f},1)
                        if strcmp(LIMO_sub.data.chanlocs(e).labels,chan.expected_chanlocs(expected).labels)
                            for f=1:size(Yr,2)
                                for c=1:size(tmpX,2)
                                    tmp(expected,f,c) = mean(W{f}(e,f,logical(tmpX(:,c))));
                                end
                            end
                        end
                    end
                end
           else
                if strcmpi(LIMO_sub.design.method,'WLS')
                    tmp = NaN(size(chan.expected_chanlocs,2),size(tmpX,2));
                elseif strcmpi(LIMO_sub.design.method,'IRLS')
                    tmp = NaN(size(chan.expected_chanlocs,2),size(Yr,2),size(tmpX,2));
                end
                for expected=1:size(chan.expected_chanlocs,2)
                    for e=1:size(W{f},1)
                        if strcmp(LIMO_sub.data.chanlocs(e).labels,chan.expected_chanlocs(expected).labels)
                            for c=1:size(tmpX,2)
                                if strcmpi(LIMO_sub.design.method,'WLS')
                                    tmp(expected,c) = mean(W{f}(e,logical(tmpX(:,c))));
                                elseif strcmpi(LIMO_sub.design.method,'IRLS')
                                    for frame = 1:size(W{f},2)
                                        tmp(expected,frame,c) = mean(W{f}(e,frame,logical(tmpX(:,c))));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        Bias{f} = tmp;
        
        if f==length(LIMO_files)
            if  strcmpi(LIMO_sub.design.method,'WLS') && length(unique(cell2mat(cellfun(@size,Bias,'UniformOutput',false)'))) ~=2 || ...
                    strcmpi(LIMO_sub.design.method,'IRLS') && length(unique(cell2mat(cellfun(@size,Bias,'UniformOutput',false)'))) ~=3
                disp('the computed bias matrices are in the workspace');  assignin('base','Bias',Bias);
                error('It is likely not all subjects have the same number of conditions, conversion from cell to mat impossible - Bias analysis stopped')
            else    
                factor_nb = size(tmpX,2);
                if strcmpi(LIMO_sub.design.method,'WLS')
                    Yr = NaN(length(limo.data.chanlocs),1,f,factor_nb);
                    for s=1:f
                        for e=1:length(limo.data.chanlocs)               
                            Yr(e,1,s,:)=Bias{s}(e,:); % single frame
                        end
                    end
                elseif strcmpi(LIMO_sub.design.method,'IRLS')
                    Yr = NaN(length(limo.data.chanlocs),size(Yr,2),f,factor_nb);
                    for s=1:f
                        for e=1:length(limo.data.chanlocs)               
                            Yr(e,:,s,:)=Bias{s}(e,:,:); 
                        end
                    end
                end
            end
            
            % stats
            disp('Testing for bias across all conditions')
            LIMO = limo; mkdir('Bias testing');
            cd('Bias testing');
            LIMO.dir        = pwd;
            LIMO.Level      = 2;
            LIMO.data_dir   = pwd;
            LIMO.data.data  = LIMO_files;
            LIMO.data.start = 1;
            LIMO.data.end   = 1;
            LIMO.data.trim1 = 0;
            LIMO.data.trim2 = 0;
            LIMO.design.electrode = [];
            LIMO.design.name      = 'Rep_ANOVA';
            LIMO.design.neighbouring_matrix = chan.channeighbstructmat;
            LIMO.design.bootstrap =1000;
            LIMO.design.tfce = 0;
            save LIMO LIMO;
            save Yr Yr; clear Bias
            limo_random_robust(6,Yr,ones(size(Yr,3),1),factor_nb,LIMO,'go','Yes');
            try close('Design matrix'); end
        end
    end % close bias
end

cd ..
disp('analysis done')
disp('Plot central tendency to check weights per subject and decile')
disp('view results ''all'' for outliers and bias')
limo_results


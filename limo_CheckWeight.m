function limo_CheckWeight(LIMO_files, expected_chanlocs, varargin)

% general function designed to look at the weights computed for each trial
% at the each channel for each subject
%
% FORMAT limo_CheckWeight(list_of_LIMO.mat,options)
%
% INPUT list_of_LIMO.mat empty [] calls a gui to select a txt file
%                        txt file listing where the LIMO.mat are located
%                        cell array listing where the LIMO.mat are located
%
%       expected_chanlocs the default channel locations for all subject
%
%       options 'plot rank' on/off since weights are between 0 and 1,
%                it computes the average response for each decile
%               'test difference' on/off compute an OLS between the good
%               trials weights = 1/0.9 and the outliers (by reverse
%               engineering the weights to outlier detection)
%               'check bias' on/off check that the weights are distributed
%               across trials in a uniform manner, i.e. that not one conditions
%               is more affected than another which would bias the results, but
%               also indicate that something is going on in the data
%
% OUTPUT creates a folder called 'Weights_checking' with the different
%        results in it
%
% Cyril Pernet 21-08-2015
% -----------------------------
% Copyright (C) LIMO Team 2015

%% if no input do it all
limo = struct('plotrank','on','testdifference','on','checkbias','on');

if nargin == 0
    [~,~,LIMO_files] = limo_get_files([],'*txt','choose a list of LIMO files');
    if isempty(LIMO_files); return; end
    [to_load,path] = uigetfile2('expected_chanlocs.mat','load chanlocs'); 
    if to_load == 0; return; end
    chan = load([path to_load]);
    limo.data.chanlocs = chan.expected_chanlocs;
    limo.data.neighbouring_matrix =  chan.channeighbstructmat;
end

%% input checks
if isempty(LIMO_files)
    [~,~,LIMO_files] = limo_get_files([],'*txt','choose a list of LIMO files');
elseif ischar(LIMO_files)
    files=textread(LIMO_files,'%s','delimiter','');  % select a txt file listing all files
    clear LIMO_files
    for f=1:size(files,1)
        LIMO_files{f} = files{f};
    end
end

% those files are there?
for f=1:length(LIMO_files)
    if ~exist(LIMO_files{f},'file');
        error([LIMO_files{f} ' doesn''t exist'])
    end
    
    [limo_paths{f},name,ext]=fileparts(LIMO_files{f});
    if ~strcmp([name ext],'LIMO.mat');
        error([LIMO_files{f} ' is not a LIMO.mat file'])
    end
    
    load(LIMO_files{f});
    if f==1
        limo.Analysis = LIMO.Analysis;
        limo.Type = LIMO.Type;
    else
        if limo.Analysis ~= LIMO.Analysis
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
        end
        
        if limo.Type ~= LIMO.Type
            error('Looks like different type of analyses (Channels/Components) are mixed up')
        end    
    end
end

if ~isfield(limo,'data')
    [to_load,path] = uigetfile2('*mat','load chanlocs'); chan = load([path to_load]);
    limo.data.chanlocs = chan.expected_chanlocs;
    limo.data.neighbouring_matrix =  channeighbstructmat;
end

if ~isempty(varargin)
    for n=1:length(varargin)
        if strcmpi(varargin{n},'plot rank')
            limo.plotrank = varargin{n+1};
        end
        
        if strcmpi(varargin{n},'test difference')
            limo.testdifference = varargin{n+1};
        end
        
        if strcmpi(varargin{n},'check bias')
            limo.checkbias = varargin{n+1};
        end
    end
end
   
%% compute
mkdir('Weights_checking'); cd('Weights_checking');
[first_frame,last_frame,subj_chanlocs,limo] = limo_match_frames(limo_paths,limo);
if strcmp(limo.plotrank,'on')
    if strcmpi(limo.Analysis,'Time-Frequency')
        data = NaN(size(chan.expected_chanlocs,2),(limo.data.highf-limo.data.lowf+1),(limo.data.trim2-limo.data.trim1+1),length(LIMO_files),10);
        difference = NaN(size(data,1),size(data,2),size(data,3),size(data,4));
    else
        data = NaN(size(chan.expected_chanlocs,2),(limo.data.trim2-limo.data.trim1+1),length(LIMO_files),10);
        difference = NaN(size(data,1),size(data,2),size(data,3));
    end
end

for f=1:length(LIMO_files)
    fprintf('reading data subject %g\n',f)
    load(LIMO_files{f}); W{f} = LIMO.design.weights;
    load([LIMO.dir filesep 'Yr.mat']);
    
    if strcmp(limo.plotrank,'on')
        if strcmpi(limo.Analysis,'Time-Frequency')
            array = find(~isnan(Yr(:,1,1,1))); % skip empty electrodes
            tmp = NaN(size(Yr,1),size(Yr,2),size(Yr,3),10);
            for e=1:length(array)
                for w=1:10
                    index = logical((W{f}(e,:) >((w/10)-1)) .* (W{f}(e,:) <=(w/10)));
                    tmp(e,:,:,w) = mean(Yr(e,:,:,index),4);
                end
                data(:,:,:,f,:) = limo_match_elec(subj_chanlocs(f).chanlocs,chan.expected_chanlocs,1,size(Yr,2),tmp);
            end
        else
            array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
            tmp = NaN(size(Yr,1),size(Yr,2),10);
            for e=1:length(array)
                for w=1:10
                    index = logical((W{f}(e,:) >((w/10)-1)) .* (W{f}(e,:) <=(w/10)));
                    tmp(e,:,w) = mean(Yr(e,:,index),3);
                end
                data(:,:,f,:) = limo_match_elec(subj_chanlocs(f).chanlocs,chan.expected_chanlocs,1,size(Yr,2),tmp);
            end
        end
        
        % save using the same format as limo_central_tendency, then call
        % limo_add_plots to make the figure - only parameters and subjedcts are
        % reversed because we want to plot per subjects and not per parameters
        if f==length(LIMO_files)
            clear Data; Data.data = data; Data.limo = limo;
            save subjects_weighted_data Data
        end
    end % close rank computation
    
    % is there a difference between the outlier trials and the best trials
    if strcmp(limo.testdifference,'on')
        if strcmpi(limo.Analysis,'Time-Frequency')
            array = find(~isnan(Yr(:,1,1,1))); % skip empty electrodes
            tmp = NaN(size(Res,1),size(Res,2),size(Res,3),10);
        else
            tmp = NaN(size(Yr,1),size(Yr,2));
            for e=1:size(Yr,1)
                out = round(W{f}(e,:) + 0.25)'; index1 = find(out == 0);
                index2 = find(W{f}(e,:) >=0.9);
                for expected=1:size(chan.expected_chanlocs,2)
                    if strcmp(LIMO.data.chanlocs(e).labels,chan.expected_chanlocs(expected).labels)
                        difference(expected,:,f) = mean(Yr(e,:,index1),3) - mean(Yr(e,:,index2),3);
                    end
                end
            end
        end
        
        if f==length(LIMO_files)
            clear Data; Data.data = difference; Data.limo = limo;
            mkdir('trial_differences'); cd('trial_differences'); 
            save Yr Data; clear difference
            
            % stats
            disp('Computing t-test between good and bad trials across all conditions')
            one_sample = NaN(size(Data.data,1), size(Data.data,2), 5);
            [one_sample(:,:,4),one_sample(:,:,1),~,one_sample(:,:,2),one_sample(:,:,5),~,one_sample(:,:,3)]=limo_trimci(Data.data);
            save ('one_sample_ttest_outliers','one_sample', '-v7.3');
            
            nboot = 1000;
            H0_one_sample = NaN(size(Data.data,1), size(Data.data,2),2,nboot); % stores T and p values for each boot under H0
            centered_data = Data.data - repmat(limo_trimmed_mean(Data.data),[1 1 size(Data.data,3)]);
            boot_table = limo_create_boot_table(Data.data,nboot);
            for electrode = 1:size(Data.data,1)
                fprintf('bootstrapping electrode %g\n',electrode)
                tmp = centered_data(electrode,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                parfor b=1:nboot
                    [t{b},~,~,~,p{b},~,~]=limo_trimci(Y(1,:,boot_table{electrode}(:,b)));
                end
                
                for b=1:nboot
                    H0_one_sample(electrode,:,1,b) = t{b};
                    H0_one_sample(electrode,:,2,b) = p{b};
                end
            end
            mkdir('H0'); save (['H0', filesep, 'H0_one_sample_ttest_outliers'],'H0_one_sample','-v7.3');
            LIMO = limo; LIMO.Level = 2; LIMO.design.bootstrap = 1000;
            LIMO.design.electrode = []; LIMO.design.name = 'one sample ttest'; save LIMO LIMO
            cd ..; clear LIMO;  load(LIMO_files{f});
        end
    end % close test difference
    
    if strcmp(limo.checkbias,'on')
        if LIMO.design.nb_conditions == 0
            sprintf('skipping bias test, there are no conditions, subject %g \n',f);
        else
            if LIMO.design.nb_interactions ~= 0
                tmpX = LIMO.design.X(:,sum(LIMO.design.nb_conditions):(size(LIMO.design.X,2)-LIMO.design.nb_continuous));
            elseif LIMO.design.nb_interactions == 0 && length(LIMO.design.nb_conditions)>1
                [tmpX,~] = limo_make_interactions(LIMO.design.X, LIMO.design.nb_conditions);
            else
                tmpX = LIMO.design.X(:,1:LIMO.design.nb_conditions);
            end
            % tmpX is always a matrix with one column per condition, no
            % matter the design - we want to know if the weights are
            % distributed equally
            if strcmpi(limo.Analysis,'Time-Frequency')
                % LIMO.data.4d
                tmp = NaN(size(W{f},1),size(Res,2),size(tmpX,2));
            else
                tmp = NaN(size(chan.expected_chanlocs,2),size(tmpX,2));
                for expected=1:size(chan.expected_chanlocs,2)
                    for e=1:size(W{f},1)
                        if strcmp(LIMO.data.chanlocs(e).labels,chan.expected_chanlocs(expected).labels)
                            for c=1:size(tmpX,2)
                                tmp(expected,c) = mean(W{f}(e,logical(tmpX(:,c))));
                            end
                        end
                    end
                end
            end
        end
        Bias{f} = tmp;
        
        if f==length(LIMO_files)
            try 
                Bias = cell2mat(Bias); 
                Bias = reshape(Bias,size(chan.expected_chanlocs,2),size(Bias,2)/f,f);
            catch SizeIssue
               disp('the computed bias matrices are in the workspace');  assignin('base','Bias',Bias); 
               error('It is likely not all subjects have the same number of conditions, conversion from cell to mat impossible')
            end
                
            % stats
            disp('Testing for bias across all conditions')
            LIMO = limo; mkdir('Bias testing'); 
            cd('Bias testing'); LIMO.data_dir = pwd;
            LIMO.data.data = LIMO_files; LIMO.data.start = 1;
            LIMO.data.end = 1; LIMO.data.trim1 = 0; LIMO.data.trim2 = 0; 
            % LIMO.design.electrode = chan.expected_chanlocs;
            LIMO.design.electrode = []; LIMO.design.name = 'Rep_ANOVA';
            LIMO.design.neighbouring_matrix = chan.channeighbstructmat;
            LIMO.Level = 2; LIMO.design.bootstrap =1000; save LIMO LIMO; 
            
            factor_nb = size(Bias,2);
            Yr = NaN(size(Bias,1),1,size(Bias,3),factor_nb);
            for e=1:size(Bias,1)
                Yr(e,1,:,:)=squeeze(Bias(e,:,:))';
            end
            save Yr Yr; clear Bias
            limo_random_robust(6,Yr,ones(size(Yr,3),1),factor_nb,LIMO,1000,0);
            try close('Design matrix'); end
        end
    end % close bias
end
cd ..
disp('analysis done')
disp('Plot central tendency to check weights per subject and decile')
disp('view results ''all'' for outliers and bias')



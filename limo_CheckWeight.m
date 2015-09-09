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
    [to_load,path] = uigetfile2('expected_chanlocs.mat','load chanlocs'); chan = load([path to_load]);
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
    else
        if limo.Analysis ~= LIMO.Analysis
            error('Looks like different type of analyses (Time/Freq/Time-Freq) are mixed up')
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
    
    if strcmp(limo.plotrank,'on')
        load([LIMO.dir filesep 'Yr.mat']);
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
    end
    
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
            save subjects_outlier_difference Data; clear difference
            % stats
            disp('Conputing t-test being good and ad trials')
            one_sample = NaN(size(Data.data,1), size(Data.data,2), 5);
            [one_sample(:,:,1),one_sample(:,:,3),~,sd,n,one_sample(:,:,4),one_sample(:,:,5)] = limo_ttest(1,Data.data,0,0.5);
            one_sample(:,:,2) = sd./sqrt(n); save ('one_sample_outliers','one_sample', '-v7.3');
            H0_one_sample=limo_ttest_permute(Data.data,1000); mkdir('H0'); 
            save (['H0', filesep, 'H0_one_sample_outliers'],'H0_one_sample','-v7.3');
            LIMO = limo; LIMO.Level = 2; LIMO.design.bootstrap = 1000; 
            LIMO.design.electrode = []; save LIMO LIMO
        end
    end
end

%$ bias
% needs to match electrode for W
% if LIMO.design.nb_conditions == 0
%     sprintf('skipping difference test, there are no conditions, subject %g \n',f);
% else
%     if LIMO.design.nb_interactions ~= 0
%         tmpX = LIMO.design.X(:,sum(LIMO.design.nb_conditions):(size(LIMO.design.X,2)-LIMO.design.nb_continuous));
%     elseif LIMO.design.nb_interactions == 0 && length(LIMO.design.nb_conditions)>1
%         [tmpX,~] = limo_make_interactions(LIMO.design.X, LIMO.design.nb_conditions);
%     else
%         tmpX = LIMO.design.X(:,1:LIMO.design.nb_conditions);
%     end
%     tmpX = [tmpX ones(size(tmpX,1),1)];
%     B(:,:,f) = pinv(tmpX)*W{f}';
% end




function filesout = limo_contrast_sessions(varargin)

% routine to compute a contrast between two sessions
% unlike limo_contrast.m we do not assume a common error - this
% is used to compare parameters from identical models but from
% different sessions, which is the same as a 2 samples t-test on 
% weighted data (if WLS or IRLS used)
%
% FORMAT fileout = limo_contrast_session(LIMO1, LIMO2, pair)
%
% INPUTS LIMO1 and LIMO2 are LIMO file names from two GLM
%                       of the same dimensions - comparisons are performed
%                       columns-wise (ie from each regressor in LIMO.design.X)
%        pair is the string indicating the session pair e.g. '12'
%
% OUTPUT filesout is a cell array of the contrast file names creates
%
% Exemple
% filesout = limo_contrast_sessions('F:\sub-1\sess-1\design1_GLM_Channels_Time_OLS\LIMO.mat',...
%       'F:\sub-1\sess-2\design1_GLM_Channels_Time_OLS\LIMO.mat', '12')
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2021

if nargin<2
    error('at least two LIMO files in are expected')
    
elseif nargin<=3
    LIMO1 = varargin{1};
    if ~isstruct(LIMO1)
        if iscell(LIMO1); LIMO1 = cell2mat(LIMO1); end
        if ischar(LIMO1)
            LIMO1 = load(LIMO1);
            LIMO1 = LIMO1.LIMO;
        else
            error('can''t load LIMO1 file')
        end
    end
    
    LIMO2 = varargin{2};
    if ~isstruct(LIMO2)
        if iscell(LIMO2); LIMO2 = cell2mat(LIMO2); end
        if ischar(LIMO2)
            LIMO2 = load(LIMO2);
            LIMO2 = LIMO2.LIMO;
        else
            error('can''t load LIMO2 file')
        end
    end
    
    if nargin ==3
        pair = varargin{3};
    else
        pair = [];
    end
    
else
    help limo_contrast_sessions
    return
end
    
%% check designs are compatible
if size(LIMO1.design.X,2) ~= size(LIMO2.design.X,2) || ...
        LIMO1.design.nb_conditions ~= LIMO2.design.nb_conditions || ...
        LIMO1.design.nb_interactions ~= LIMO2.design.nb_interactions  || ...
        LIMO1.design.nb_continuous ~= LIMO2.design.nb_continuous
    error('Designs don''t match')
else
    LIMO = LIMO1; % need a 'local' copy for group stats 
end

% since those are sessions, they must have a common folder location
LIMO.dir              = fileparts(LIMO1.dir(1:find(LIMO1.dir ~= LIMO2.dir)-1));
LIMO.data             = rmfield(LIMO.data,'data');
LIMO.data.data{1}     = fullfile(LIMO1.dir,'Yr.mat');
LIMO.data.data{2}     = fullfile(LIMO2.dir,'Yr.mat');
LIMO.data             = rmfield(LIMO.data,'data_dir');
LIMO.data.data_dir{1} = LIMO1.dir;
LIMO.data.data_dir{2} = LIMO2.dir;
LIMO.data.Cat         = [];
LIMO.data.Cont        = [];
LIMO.design.X         = [];

%% for each regressor compute the difference

% load the data
try
    LIMO.data.data{1} = fullfile(LIMO1.dir,'Yr.mat');
    data1 = load(fullfile(LIMO1.dir,'Yr.mat')); data1 = data1.Yr;
    if strcmp(LIMO1.Analysis,'Time-Frequency') || strcmp(LIMO1.Analysis,'ITC')
        data1 = limo_tf_4d_reshape(data1);
    end
    
    LIMO.data.data{2} = fullfile(LIMO2.dir,'Yr.mat');
    data2 = load(fullfile(LIMO2.dir,'Yr.mat')); data2 = data2.Yr;
    if strcmp(LIMO2.Analysis,'Time-Frequency') || strcmp(LIMO2.Analysis,'ITC')
        data2 = limo_tf_4d_reshape(data2);
    end
    
    if any(size(data1,[1 2]) ~= size(data2,[1 2]))
        error('Dataset have different dimensions - reprocess sessions with the same parameters')
    end
catch no_load
    error('could not load datasets? \n%s',no_load.message)
end

% ensure match on dimension 1 by getting joint labels in the same order
LIMO.data = rmfield(LIMO.data,'chanlocs');
joint_labels   = intersect(arrayfun(@(x)x.labels,LIMO1.data.chanlocs,'UniformOutput',false),...
    arrayfun(@(x)x.labels,LIMO2.data.chanlocs,'UniformOutput',false),'stable');

% compute for each regressor
if ~isempty(joint_labels)
    for regressors = (size(LIMO1.design.X,2)-1):-1:1
        con      = NaN(size(data1,1),size(data1,2),size(data1,3),5); % dim 3 = mean diff/se/df/t/p
        filesout{regressors} = fullfile(LIMO.dir,['con_' num2str(regressors) 'sess_' pair '.mat']);
        
        for e = 1:length(joint_labels)
            channel = find(arrayfun(@(x) strcmp(x.labels,joint_labels(e)),LIMO1.data.chanlocs));
            tmp     = squeeze(data1(channel,:,:)); % all frames and trials
            index   = intersect(find(~isnan(tmp(1,:))),find(LIMO1.design.X(:,regressors))); % only valid trials for that regressor
            Y1      = tmp(:,index); clear tmp
            if strcmpi(LIMO1.design.method,'WLS')
                W  = squeeze(LIMO1.design.weights(channel,index));
                Y1 = Y1.*repmat(W,[size(Y1,1) 1]);
            elseif strcmpi(LIMO1.design.method,'IRLS')
                W  = squeeze(LIMO1.design.weights(channel,:,index));
                Y1 = Y1.*W;
            end
            
            channel = find(arrayfun(@(x) strcmp(x.labels,joint_labels(e)),LIMO2.data.chanlocs));
            tmp     = squeeze(data2(channel,:,:));
            index   = intersect(find(~isnan(tmp(1,:))),find(LIMO2.design.X(:,regressors)));
            Y2      = tmp(:,index); clear tmp
            if strcmpi(LIMO2.design.method,'WLS')
                W  = squeeze(LIMO2.design.weights(channel,index));
                Y2 = Y2.*repmat(W,[size(Y2,1) 1]);
            elseif strcmpi(LIMO2.design.method,'IRLS')
                W  = squeeze(LIMO2.design.weights(channel,:,index));
                Y2 = Y2.*W;
            end
            
            % we want mean diff/se/df/t/p and we get mean/dfe/ci/std/n/t/p
            LIMO.data.chanlocs(channel) = LIMO2.data.chanlocs(channel); 
            [con(channel,:,1),con(channel,:,3),~,sd,~,con(channel,:,4),...
                con(channel,:,5)]=limo_ttest(2,Y1,Y2,.05);
            sd = sd.^2; a = sd(1,:)./size(Y1,2); b = sd(1,:)./size(Y2,2);
            con(channel,:,2) = sqrt(a + b);
            clear Y1 Y2
        end
        if strcmp(LIMO1.Analysis,'Time-Frequency') ||  strcmp(LIMO1.Analysis,'ITC')
            con = limo_tf_4d_reshape(con);
        end
        save (filesout{regressors},'con', '-v7.3')
    end
    save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
else
    filesout = [];
end




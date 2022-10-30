function [X, nb_conditions, nb_interactions, nb_continuous] = limo_design_matrix_tf(varargin)

% LIMO_DESIGN_MATRIX_TF - called once by LIMO_EEG.m to create the design matrix X
% that is stored in the LIMO.mat file and called for all analyses.
%
% This is a reduced version of limo_design_matrix, with specialised
% handling for 4D time-frequency (ersp) data.
%
% FORMAT
% [X,nb_conditions,nb_interactions,nb_continuous] = limo_design_matrix_tf(Y,LIMO,flag)
% [X,nb_conditions,nb_interactions,nb_continuous] = limo_design_matrix_tf(Y,Cat,Cont,directory,zscoring,full_factorial,flag)
%
% INPUTS:
% Y = EEG data with format electrodes x frames x timewindows x trials/subjects
% LIMO = structure that contains the above information (created by limo_egg and limo_import_tf)
% or input info that are in LIMO.mat
% Cat = vector describing the different conditions (if no conditions Cat = 0)
% Cont = matrix describing the different covariates (if no covariates Cont = 0)
% directory = path of folder where the outputs will be saved (see below)
% zscoring = [0/1] - if 1 (default) continuous regressors are zscored,
% which means that the betas coefficients will have units
% micro-volts per std of the predictor variable.
% full_factorial = [0/1] - if 1 create interaction terms between the
% factors described in Cat
% flag = figure on/off
%
% OUTPUTS:
% X = 2 dimensional matrix that describes the experiments' events
% nb_conditions = vector that returns the number of conditions per factor e.g. [2 2 2]
% nb_interactions = vector that returns the number of conditions perinteraction e.g. [4 4 4 8]
% nb_continuous = scalar that returns the number of continuous variables e.g. [3]
%
% These outputs are written to disk in DIRECTORY (allows dirty check that memery holds) and populated latter:
%
% Yr.mat = NaNs - the EEG data from the .set (reorganized to fit X, that is grouped by conditions if Cat ~=0)
% Yhat.mat = NaNs - the predicted data (same size as Yr)
% Res.mat = NaNs - the residual (non modeled) data (same size as Yr)
% R2.mat = NaNs - the model fit (same size as Yr)
% Beta.mat = NaNs - the beta values (dim: electrode,frame, number of paramters in the model)
%
% See also LIMO_IMPORT, LIMO_EEG, LIMO_GLM
%
% Andrew Stewart v1 Nov 13 % Adapted from Cyril's limo_design_matrix
% ------------------------------
%  Copyright (C) LIMO Team 2019


%% varagin stuffs

type_of_analysis = 'Mass-univariate'; % default to create R2 file, updated using the LIMO structure only

if nargin==3
    Y                = varargin{1};
    Cat              = varargin{2}.data.Cat;
    Cont             = varargin{2}.data.Cont;
    directory        = varargin{2}.dir;
    zscoring         = varargin{2}.design.zscore;
    full_factorial   = varargin{2}.design.fullfactorial;
    chanlocs         = varargin{2}.data.chanlocs;
    type_of_analysis = varargin{2}.design.type_of_analysis;
    size4D           = varargin{2}.data.size4D;
    flag             = varargin{3};
    try
        expected_chanlocs = varargin{2}.data.expected_chanlocs;
    catch ME
        expected_chanlocs = [];
    end
elseif nargin==7 || nargin==9 % allows running without specifying chanlocs and neighbourging matrix
    Y              = varargin{1};
    Cat            = varargin{2};
    Cont           = varargin{3};
    directory      = varargin{4};
    zscoring       = varargin{5};
    full_factorial = varargin{6};
    flag           = varargin{7};
    try
        chanlocs = varargin{8};
        expected_chanlocs = varargin{9};
    catch ME
        chanlocs = [];
        expected_chanlocs = [];
    end
else
    error('wrong number of arguments')
end


%% Dimension check
Y = limo_tf_4d_reshape(Y);

if Cont == 0
    Cont = [];
end

if Cat == 0
    Cat = [];
end

% check Cat
if size(Cat,2) == size(Y,3)
    Cat = Cat';
    disp('Cat has been transposed')
end

% check Cont
if size(Cont,2) == size(Y,3)
    Cont = Cont';
    disp('Cont has been transposed')
end

% overall dimensions check
if isempty(Cat) && ~isempty(Cont)
    if size(Y,3) ~= size(Cont,1)
        error('The number of trials and the covariate(s) length are of different')
    end
elseif isempty(Cont) && ~isempty(Cat)
    if size(Y,3) ~= size(Cat,1)
        error('The number of trials and the number of events are different size')
    end
elseif ~isempty(Cat) && ~isempty(Cont)
    if size(Y,3) ~= size(Cont,1)
        error('The number of trials and the covariate(s) length are of different')
    elseif size(Y,3) ~= size(Cat,1)
        error('The number of trials and the number of events are different size')
    elseif length(Cat) ~= size(Cont,1)
        error('The number of events and the covariate(s) length are different')
    end
end
 
% additional checking for regressions
if isempty(Cat) == 0

    if size(Cont,2)+1 >= size(Y,3)
        error('there are too many regressors for the number of trials, reduce the model size')
        % one cannot compute any statistics if size(Y,3) > size(Cont,2)+1
    end

    if size(Cont,2) >1
        if det(Cont'*Cont) == 0
            errordlg('the regression matrix is singular, 1 or more predictors are function to each other, please consider another model', ...
                'Matrix definition error'); return
        elseif cond(Cont'*Cont,1) == Inf
            errordlg('the regression matrix is close to singular, please consider another model','Matrix definition error'); return
        end

    end
end

% final checking for full factorial
if full_factorial == 1 && ~isempty(Cont)
    error('LIMO does not compute full factorial ANCOVAs with covariate(s), consider a factorial ANCOVA witout interaction')
end

cd(directory);

%% check if NaNs - that is offer the possibility to remove some trials
if ~isempty(Cat)
    check        = find(sum(isnan(Cat),2));
    Cat(check,:) = [];
    if ~isempty(Cont)
        Cont(check,:) = [];
    end
    Y(:,:,check) = [];
end

if ~isempty(Cont)
    check         = find(sum(isnan(Cont),2));
    Cont(check,:) = [];
    if ~isempty(Cat)
        Cat(check,:) = [];
    end
    Y(:,:,check) = [];
end

%% Make the design matrix and create files

if isempty(Cont)
    add_cont = 0; % no continuous regressors
else
    add_cont = 1; % add continuous regressors
end

if isempty(Cat)
    nb_factors = 0;
else
    nb_factors = size(Cat,2); % stands for the nb of factors
end
nb_conditions = 0; % stands for the nb_conditions in each factor
nb_interactions = 0; % stands for the nb of interactions between factors
nb_continuous = 0; % stands for the nb of continuous regressors

disp('making up the design matrix and data files ...')

% deal with covariates
if add_cont == 1

    [l,w]=size(Cont);
    nb_continuous = w;

    if l~=size(Y,3)
        errordlg('the covariate(s) must be the same length as your dependant variable(s)')
        disp('please retry changing the 3rd argument');
    end

    % to make continuous regressors comparable data are zscored
    if zscoring == 1
        for i=1:w
            Cont(:,i) = (Cont(:,i)-mean(Cont(:,i)))./std(Cont(:,i));
        end
    end

    if isempty(Cat)
        X = zeros(l,w+1);
        X(:,1:w) = Cont;
        X(:,w+1) = 1;
        Yr = Y;
    end
        
end % closes if Cont~=0


% Cat is a vector or a matrix and we want X
% if Cont is used; it is reorganized accordingly
% Y is also reorganized accordingly (i.e. grouped per condition)

if ~isempty(Cat)
    
    if size(Cat,1)~=size(Y,3)
        errordlg('the categorical regressor must be the same length as your dependant variables')
        disp('please retry changing the 2nd argument');error('error line 228 in limo_design_matrix');
    end

    % sort Cat column wise - apply to Cont and Y
    [Cat,index] = sortrows(Cat,[1:size(Cat,2)]);
    
    if add_cont == 1
        for i=1:size(Cont,2)
            Cont(:,i) = Cont(index,i);
        end
    end
    
    Yr = NaN(size(Y)); % Y reorganized
    for electrode = 1:size(Y,1)
        for time = 1:size(Y,2)
            Yr(electrode,time,:) = Y(electrode,time,index);
        end
    end
    clear Y
    
    % find out the number of trial per condition per factor
    % create as many columns as conditions
    condition_count = 0;
    index_nb_conditions =1;
    index_indices_conditions =1;
    for c=1:size(Cat,2)
        column = Cat(:,c);
        for i=min(column):max(column)
            tmp = find(column==i);
            if ~isempty(tmp) % use this trick so that any number can be used as event marker
                nb_conditions(index_nb_conditions) = condition_count+1;
                condition_count = condition_count+1;
                indices_conditions{index_indices_conditions}=tmp;
                index_indices_conditions = index_indices_conditions+1;
            end
        end
        index_nb_conditions = index_nb_conditions+1;
        condition_count = 0;
    end

    % create the design matrix for the main conditions
    x = zeros(size(Yr,3),length(indices_conditions));
    for f=1:length(indices_conditions)
        x(indices_conditions{f},f) = 1;
    end
    % basic_design = x;
    
    % add interactions if requested
    if full_factorial == 1
        if nb_factors == 1
            disp('full factorial impossible, only 1 factor specificied')
        else
            [x, nb_interactions] = limo_make_interactions(x, nb_conditions);
        end
    end
            
    % add the continuous regressors and the constant
    if add_cont == 1
        X = [x Cont ones(size(Yr,3),1)];
    else
        X = [x ones(size(Yr,3),1)];
    end
    
    % finally check the balance of the design and adjust X and Yr
    if full_factorial == 1
        start_at = sum(nb_conditions) + sum(nb_interactions(1:end-1)) + 1;
        higher_interaction = X(:,start_at:(end-nb_continuous-1));
        if size(higher_interaction,2) ~= prod(nb_conditions)
            errordlg(sprintf('the design is too unbalanced to be corrected \n can''t run full factorial'))
            nb_interactions = 0;
            X = [basic_design ones(size(Yr,3),1)];
        else
            nb_trials = sum(higher_interaction);
            if length(unique(nb_trials)) > 1
                sample_to_n = min(nb_trials);
                if sample_to_n == 1
                    errordlg(sprintf('the design is over-specified, \n only one observation per condition \n can''t run full factorial'))
                    nb_interactions = 0;
                    X = [basic_design ones(size(Yr,3),1)];
                else % sample
                    s = 1; % now list rows to keep
                    for c=1:size(higher_interaction,2)
                        indices = find(higher_interaction(:,c));
                        if length(indices) > sample_to_n
                            new_sample{s} = randsample(indices,sample_to_n);
                            s = s+1;
                        else
                            new_sample{s} = indices;
                            s = s+1;
                        end
                    end
                    % sample
                    sample_index = [];
                    for s = 1:size(new_sample,2)
                        v = new_sample{s};
                        sample_index = [ sample_index ; v];
                    end
                    X = X(sample_index,:);
                    Yr = Yr(:,:,sample_index);
                    disp('DATA HAVE BEEN RESAMPLED TO PERFORM A FULL FACTORIAL ANALYSIS')
                end
            end
        end
    end
end % closes if Cat ~=0

% if Cat and Cont are empty (ie just comnpute the mean)
if ~exist('Yr','var') 
    Yr = Y;
    X = ones(size(Yr,3),1);
end
clear Y

% --------------------------------------------------------
% create files to be used by LIMO in all cases - dim = expected chanlocs

% Data - update with expected chanlocs
if ~isempty(expected_chanlocs) && size(Yr,1) > 1
    Yr = limo_match_elec(chanlocs,expected_chanlocs,1,size(Yr,2),Yr);
end

% create and hold in memory all matrices to be - allows a crude check a
% memory abilities - (in the next version - catch the error to run slowly
% looping per frequency rathen than holding large matrices)
try
    
    % no matter the analysis we have Beta, Yhat, Res
    Yhat = zeros(size4D);
    Res = zeros(size4D);
    Betas = zeros(size4D(1),size4D(2),size4D(3),size(X,2));
    
    % only for univariate analyses
    if strcmp(type_of_analysis,'Mass-univariate')
        R2 = zeros(size4D(1),size4D(2),size4D(3),3); save R2 R2 -v7.3;
    end
    
    if nb_conditions ~=0
        tmp_Condition_effect = NaN(size(Yr,1),size(Yr,2),length(nb_conditions),2);
    end
    
    if nb_interactions ~=0
        tmp_Interaction_effect = NaN(size(Yr,1),size(Yr,2),length(nb_interactions),2);
    end
    
    if nb_continuous ~=0
        tmp_Covariate_effect = NaN(size(Yr,1),size(Yr,2),nb_continuous,2);
    end
    
    % note - overwritten in limo_eeg_tf
    disp('saving 4D files to disk, be patient ...')
    Yr=limo_tf_4d_reshape(Yr);
    save Yr Yr -v7.3; clear Yr
    save Yhat Yhat -v7.3; clear Yhat
    save Res Res -v7.3; clear Res
    save Betas Betas -v7.3; clear Betas
    if strcmp(type_of_analysis,'Mass-univariate'); clear R2; end
    if nb_conditions ~=0; clear tmp_Condition_effect; end
    if nb_interactions ~=0; clear tmp_Interaction_effect; end
    if nb_continuous ~=0; clear tmp_Covariate_effect; end
    
catch FileError
    sprintf(FileError.message)
    error('error while memory mapping futur results - it is likely you don''t have enough RAM')
end
    
% ------
% figure
if flag == 1
    figure('Name','LIMO design','Color','w','NumberTitle','off')
    Xdisplay = X;
    if add_cont == 1
        cat_column = sum(nb_conditions)+sum(nb_interactions);
        REGdisplay = Cont + abs(min(Cont(:)));
        Xdisplay(:,cat_column+1:size(X,2)-1) = REGdisplay ./ max(REGdisplay(:));
    end
    imagesc(Xdisplay); colormap('gray'); drawnow;
    title('Design matrix'); xlabel('regressors');ylabel('trials');
    set(gca,'XTick',1:size(X,2))
end
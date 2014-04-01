function new_cont = limo_split_continuous(cat,cont)

% simple routine to create normalized continuous regressors (1st level
% analysis) split per condition -- such regressors are useful if one is
% interested in testing if the continuous variable has different effects
% according to conditions at hand (ie it demonstrates an interaction)
%
% FORMAT limo_split_continuous
%        limo_split_continuous(cat,cont)
%
% INPUT if no input the user is prompted to select a categorical variable
% file (*.mat or *.txt) and a continous variable file (*.mat or *.txt)
%       cat and cont if given as input are 2 vectors
%
% OUPUT new_cont is the continuous regressor normalized and split per
% condition - if no output specified, user is prompted to save
%
% Cyril Pernet 01-April-2014
% -------------------------------
% Copyright (C) LIMO Team 2014

%% INPUT CHECK
if nargin == 0
    [cat,cpath,filt]=uigetfile({'*.mat','MAT-file (*.mat)'; '*.txt','Text (*.txt)'},...
        'Select your categorical variable'); CAT = load([cpath cat]);
    if strcmp(cat(end-3:end),'.mat')
        CAT = getfield(CAT,cell2mat(fieldnames(CAT)));
    end
    
    [cont,cpath,filt]=uigetfile({'*.mat','MAT-file (*.mat)'; '*.txt','Text (*.txt)'},...
        'Select your categorical variable'); CONT = load([cpath cont]);
    if strcmp(cont(end-3:end),'.mat')
        CONT = getfield(CONT,cell2mat(fieldnames(CONT)));
    end
    
elseif nargin ==2
    CAT = varargin{1};
    CONT = varargin{2};
else
    error('2 arguments in expected')
end

if ~isvector(CAT) || ~isvector(CONT)
    error('single vectors are expected for the categorical and continuous variables')
end

if size(CAT,2) ~=1; CAT = CAT'; end
if size(CONT,2) ~=1; CONT = CONT'; end
if length(unique(CAT)) > length(unique(CONT))
    disp('Categorical and continuous regressors likely inverted - reversing them')
    tmp = CAT; CAT=CONT; CONT=tmp; clear tmp
end

%% DO THE STUFF
cat_index = unique(CAT);
% new_cat = zeros(length(CAT),length(unique(CAT)));
new_cont = zeros(length(CAT),length(unique(CAT)));
for n=1:length(unique(CAT))
    tmp = zscore(CONT(CAT==cat_index(n)));
    new_cont(CAT==cat_index(n),n) = tmp;
    % new_cat(CAT==cat_index(n),n) = 1;
end

% figure('Name','Design')
% imagesc([new_cat new_cont]); 
% colormap('gray')

%% possibly save
if nargout == 0
    name = {'new_cont'};
    uisave(name,'split_continuous_regressor')
end


    
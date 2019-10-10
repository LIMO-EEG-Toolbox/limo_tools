function go = limo_contrast_checking(varargin)

% Routine to check contrasts are valid
%
% FORMAT
% go = limo_contrast_checking(C,X); % <-- check the statistical validity of C
% go = limo_contrast_checking(LIMO.dir); % <-- add zeros to the last coded contrast
% go = limo_contrast_checking(LIMO.dir, LIMO.design.X, C);
%
% INPUT
% C a vector or matrix of contrasts
% X the design matrix (also LIMO.design.X)
% LIMO.dir the directory where LIMO.mat is
%
% OUTPUT
% if input C and X; go=1 if valid contrast or 0 if invalid
% if input dir, X and C, go = C corrected with extra 0s
%
% Cyril Pernet, v5. 2019
% Valid constrasts are a sum (ones), a weighted sum (ones/N) 
% or projection invariant (orthogonal usualy).
% -----------------------------------------------------------
%  Copyright (C) LIMO Team 2019


%% update the contrast with 0s
% ----------------------------------------
if nargin == 1 || nargin == 3

    limo_path = varargin{:,1};
    if ~exist(limo_path,'dir') 
        error('%s doesn''t exist',limo_path)
    end
    
    if nargin == 1
        LIMO = load(fullfile(limo_path,'LIMO.mat'));
        X = LIMO.LIMO.design.X; 
        C = LIMO.LIMO.contrast{end}.C;
    else
        X = varargin{:,2};
        C = varargin{:,3};
    end
        
    [l,w]=size(C);
    if w ~= size(X,2)
        % could be that need transpose
        if l == size(X,2)
            C = C';
            disp('C has been transposed')
            [~,w]=size(C);
        end

        % could be that cst term not coded
        if w < size(X,2)
            tmp = zeros(size(C,1),size(X,2));
            tmp(:,1:(size(C,2))) = C;
            C = tmp; clear tmp;
            disp('zeros added for the constant term')
            [~,w]=size(C);
        end

        if w ~= size(X,2)
            disp('c must have the same numner of columns as X')
            error('dimensions must agree')
        end
    end
    
    go = C;
    if nargin == 1
        LIMO = LIMO.LIMO;
        LIMO.LIMO.contrast{end}.C = C;
        save(fullfile(limo_path,'LIMO.mat'),'LIMO');
    end
   
    
%% check if contrast is valid
% --------------------------------------    
elseif nargin == 2
    
    X = varargin{2};
    for s=1:size(varargin{1},1)
        C = varargin{1}(s,:);
        N = sum(X(:,C~=0)); % N per condition 

        if sum(C) == length(find(C)) 
            go = 1; % if contrast with ones only to add parameters
        elseif ~all(int16(C(C~=0) - N/sum(N))) 
            go =1; % also a sum but equal 0 weighted by total number of observations
        else
%             check = int16(C*single((pinv(X'*X))*(X'*X)));
%             check = (check == C);
%             if sum(check) ~= length(C)
%                 go = 0; % if contrast not ones nor invariant
%             else
%                 go = 1;
%             end
            
            P = X*pinv(X); % projection
            lambda = X*C'; % contrast
            check = int16(P*lambda); % contrast onto X
            check = (check == int16(lambda)); % it is invariant
            if sum(check) ~= size(X,1)
                go = 0; % if contrast invariant
            else
                go = 1;
            end
        end
    end
else
    error('the number of arguments must be 2 or 3 ; error in limo_contrast_checking')
end

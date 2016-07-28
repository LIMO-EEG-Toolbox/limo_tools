function go = limo_contrast_checking(varargin)

% quick routine to check contrasts are ok
% used to check the validity of C
% also used to make sure contrasts have the right size
%
% FORMAT
% go = limo_contrast_checking(C,X)
% go = limo_contrast_checking(LIMO.dir, LIMO.design.X, C);
%
% INPUT
% C a vector or matrix of contrasts
% X the design matrix
%
% OUTPUT
% if input C and X; go=1 if valid contrast or 0 if invalid
% if input dir, X and C, go = C corrected with extra 0s
%
% Cyril Pernet, v4. 25-04-2010
% -----------------------------
%  Copyright (C) LIMO Team 2010


% 3 inputs = update the contrast 
% ------------------------------
if nargin == 3

    try
        cd (varargin{:,1});
    catch
            disp('can''t CD, using current directory')
    end
    load Yr; Y = Yr; clear Yr;
    
    X = varargin{:,2};
    C = varargin{:,3};
    [l,w]=size(C);

    if w ~= size(X,2)
        % could be that need transpose
        if l == size(X,2)
            C = C';
            disp('C has been transposed')
            [l,w]=size(C);
        end

        % could be that cst term not coded
        if w < size(X,2)
            tmp = zeros(size(C,1),size(X,2));
            tmp(:,1:(size(C,2))) = C;
            C = tmp; clear tmp;
            disp('zeros added for the constant term')
            [l,w]=size(C);
        end

        if w ~= size(X,2)
            disp('c must have the same numner of columns as X')
            error('dimensions must agree')
        end
    end
    
    go = C;
   
    
% 2 inputs = check if contrast is valid
% --------------------------------------    
elseif nargin == 2
    
    X = varargin{2};
    for s=1:size(varargin{1},1)
        C = varargin{1}(s,:);
        
        if sum(C) == length(find(C))
            go = 1; % if contrast with ones only to add parameters
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
                go = 0; % if contrast not ones nor invariant
            else
                go = 1;
            end
        end
    end
else
    error('the number of arguments must be 2 or 3 ; error in limo_contrast_checking')
end

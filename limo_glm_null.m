function null_y = limo_glm_null(y,X,nb_conditions,nb_interactions)

% routine for limo_glm and limo_contrast
% nullify Y based in the design
%
% FORMAT limo_glm_null(Y,LIMO)
%
% INPUT Y is the data matrix
%       X is the design matrix
%       nb_X are the design parameters 
% OUTPUT null_y is your H0 Y data
%
% Cyril Pernet - recycled from limo_glm_boot
% -------------------------------------------
% Copyright (C) LIMO Team 2021

% ---------------
%% Make null data
% ---------------

% if categorical variables, center data 1st overwise nothing to do
% -------------------------------------------------------------
if nb_conditions ==0
    null_y = y(randperm(size(y,1),size(y,1)),:); % just shuffled then resample
else
    null_y = NaN(size(y,1),size(y,2));
    if nb_interactions ~= 0
        % look up the last interaction to get unique groups
        if length(nb_interactions) == 1
            start_at = sum(nb_conditions);
        else
            start_at = sum(nb_conditions)+sum(nb_interactions(1:end-1));
        end
        
        for cel=(start_at+1):(start_at+nb_interactions(end))
            index = find(X(:,cel));
            null_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),length(index),1);
        end
        
    elseif size(nb_conditions,2) == 1 
        % no interactions because just 1 factor
        for cel=1:nb_conditions
            index = find(X(:,cel));
            null_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),length(index),1);
        end
        
    elseif size(nb_conditions,2) > 1
        % create fake interaction to get groups
        [tmpX, interactions] = limo_make_interactions(X(:,1:sum(nb_conditions)), nb_conditions);
        if length(interactions) == 1
            start_at = sum(nb_conditions);
        else
            start_at = sum(nb_conditions)+sum(interactions(1:end-1));
        end
        
        for cel=(start_at+1):(start_at+interactions(end))
            index = find(tmpX(:,cel));
            null_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),size(y(index,:),1),1);
        end
    end
end


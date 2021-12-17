function limo_review(varargin)

% Routine to display the design X with correlations
% 
% FORMAT limo_review
%        limo_review(LIMO)
%
% Cyril Pernet 
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% Varargin

if isempty(varargin)
    [file,dir] = uigetfile('LIMO.mat','select a LIMO.mat file');
    if file == 0
        return
    else
        try
            cd (dir); load LIMO.mat;
        catch
            error('file not supported')
        end
    end
else
   LIMO = cell2mat(varargin); 
end

%% Display
figure('Name','Review Design')
set(gcf,'Color','w');

if size(LIMO.design.X,2) == 1
    imagesc(LIMO.design.X/2); colormap(gca, gray);
    title('Design with the constant only (=mean)')
else
    % display a scaled version of X
    add_subplots = 1;
    if ~isempty(LIMO.design.X)
        if isfield(LIMO.design, 'weights')
            W = LIMO.design.weights;
        else
            W = ones(size(LIMO.design.X,1),1);
        end
        
        if ndims(W) == 2
            W = mean(W,1)';
        elseif ndims(W) == 3
            W = squeeze(mean(mean(W,1),2));
        end
        X = LIMO.design.X.*repmat(W,1,size(LIMO.design.X,2));
        Xdisplay = LIMO.design.X;
        N = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
        if isfield(LIMO.design, 'nb_continuous')
            if  prod(LIMO.design.nb_continuous) ~= 0
                REGdisplay = Xdisplay(:,N+1:size(X,2)-1);
                REGdisplay = REGdisplay ./ max(abs(min(REGdisplay)));
                Xdisplay(:,N+1:size(X,2)-1) = REGdisplay ./ max(max(REGdisplay));
            end
        end
        subplot(3,3,[1 2 4 5]); imagesc(Xdisplay); colormap(gca, gray);
        title('Design matrix','FontSize',14); ylabel('trials / subjects');caxis([0 1+eps])
    else
        if strncmp(LIMO.design.name,'one sample',10)
            add_subplots = 0; image(ones(size(LIMO.data.data,2),1)); colormap(gca, gray); caxis([0 1])
            title('Design matrix','FontSize',14); ylabel('trials / subjects'); set(gca,'XTicks','1')
        end
    end
    
    if add_subplots == 1
        % add the covariance matrix
        subplot(3,3,[3 6]); C = cov(X); imagesc(C); r = min(C(:))-max(C(:)); colormap(gca, parula);
        title('Covariance matrix','FontSize',14); xlabel('regressors');caxis([r+min(C(:)) -r])
        
        % add the orthogonality matrix
        orth_matrix = eye(size(X,2));
        combinations = nchoosek([1:size(X,2)],2);
        for i=1:size(combinations,1)
            orth_matrix(combinations(i,1),combinations(i,2)) = abs(X(:,combinations(i,1))'*X(:,combinations(i,2))) / (norm(X(:,combinations(i,1)))*norm(X(:,combinations(i,2))));
            orth_matrix(combinations(i,2),combinations(i,1)) = orth_matrix(combinations(i,1),combinations(i,2));
        end
        subplot(3,3,[7 8]); imagesc(orth_matrix); colormap(gca, gray);
        title('Orthogonality matrix','FontSize',14); xlabel('regressors'); caxis([0 1+eps])
    end
end

% if the matrix has been estimated we can also compute useful info
% http://en.wikipedia.org/wiki/Variance_inflation_factor


function [Selection, Tested_Models] = LIMO_model_selection (LIMO,output)

% Select the best possible design matrix or matrices subset.
% Following Hilton, 1983 (in Rencher, 2002) one uses C to select a
% subset of regressors (Monte Carlo simulation in univariate case)
% C as defined by Myers (1990) is the expected squarred error, i.e.
% the var(y)+bias due to the fact the X isn't perfect and C =
% SSE/Var(y) - (n-2p) // here we use a multivariate extension and
% the tr(C) as criteria - this is a preliminary function and not used yet
%
% Cyril Pernet, 23 Feb. 2009
% -----------------------------
%  Copyright (C) LIMO Team 2010

cd(LIMO.dir)
X = LIMO.design.X;
nb_conditions = LIMO.design.nb_conditions;
nb_continuous = LIMO.design.nb_continuous;

% selection of one electrode (max R2)
load R2; V = squeeze(R2(:,:,1)); clear R2
[v1,p1]=max(V'); [v2,p2]=max(v1);
LIMO.multivariate.model_selection.ref_electrode = p2;
LIMO.multivariate.model_selection.R2_value = v2;
LIMO.multivariate.model_selection.ref_time = p1(p2);
load Yr; Y=squeeze(Yr(p2,:,:))'; clear Yr

   
% get the usual sum of squares
T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
R        = eye(size(Y,1)) - (X*pinv(X));                                      % Residual matrix
E        = (Y'*R*Y);                                                          % SS Error
Betas    = pinv(X)*Y;                                                         % Parameters of the model

% check dof
p = size(X,2) - 1;
n = size(Y,1);
m = size(Y,2);

% criteria for the full initial matrix
Sk = E/(size(Y,1)-size(X,2));
I = eye(size(inv(Sk)*E,1),size(inv(Sk)*E,2));
Ck = (inv(Sk)*E)+(2*p-n)*I;
CriterionTrace = trace(Ck);
CriterionValue = p*m;
Tested_Models{1} = X;
if output == 1
        color_coding{1} = [nb_conditions,nb_continuous];
end


% -----------------------------
% find all possible subset of X
% -----------------------------

% 1 --> define the number of 'columns' to combine
if nb_conditions == 0 && nb_continuous ~= 0
    nb_columns = nb_continuous;
elseif nb_conditions ~=0 && nb_continuous ~= 0
    nb_columns = (nb_continuous+1);  % because we do not 'permute' columns of cat but all columns together
elseif nb_conditions ~=0 && nb_continuous == 0
    error ('there are no regressors to permute, the model cannot be reduced, error line 44 in LIMO_model_selection');
end

nb_permutation = factorial(nb_columns); % this is the real number of permutation to do

% 2 --> get the different possible permutations
for s=1:nb_columns
    comb{s}=nchoosek((nb_conditions):(nb_conditions+nb_columns-1),s);
end

index = 1;
for i=(nb_columns-1):-1:1   % start at nb_columns-1 as the 1st one is the full model (Ck)
    for j=1:length(comb{i})
        subset{index} = sort(comb{i}(j,:)); % same as above but easier to read
        index = index +1;
    end
end

% 3 --> compute the criteria C using E and R
Cst = X(:,size(X,2));

for s=1:nb_permutation

    % loop for having all conditions (gp) together
    keep_gp = subset{s}(1);
    if keep_gp == nb_conditions
        keep_gp = 1:nb_conditions;
    else
        keep_gp = 0;
    end

    % get other columns
    if keep_gp == 0
        keep_cov = subset{s}(1:size(subset{s},2));
        Xp = [X(:,keep_cov) Cst];
    else
        keep_cov = subset{s}(2:size(subset{s},2));
        Xp = [X(:,keep_gp) X(:,keep_cov) Cst];
    end

    Tested_Models{s+1} = Xp;
    Rp  = eye(size(Y,1)) - (Xp*pinv(Xp));
    Ep  = (Y'*Rp*Y);
    I = eye(size(inv(Sk)*Ep,1),size(inv(Sk)*Ep,2));
    Cp = (inv(Sk)*Ep)+(2*p-n)*I;
    CriterionTrace(s+1) = trace(Cp);
    CriterionValue(s+1) = (size(Xp,2) - 1)*m;

    if output == 1
        nb_gp = sum(keep_gp);
        if nb_gp ~= 0
            nb_gp = nb_conditions;
        end
        nb_cov =length(keep_cov);
        color_coding{s+1} = [nb_gp,nb_cov];
    end

end

Selection = CriterionTrace < CriterionValue;


% -----------------------------
% graphical output
% -----------------------------
if output == 1

    figure
    even_nb = ceil(size(Selection,2)/ 2);
    nb_rows = (even_nb /4)*2; % because 4 design matrix per row
    diff = CriterionTrace - CriterionValue; [value,position]=min(diff);
    
    model_nb = 1; make_plot = 1; i = 1;
    while make_plot == 1
        subplot(nb_rows, 4, i);
        Xdisplay = Tested_Models{model_nb};
        model_nb = model_nb +1;

        if  color_coding{i}(2) ~= 0;
            REGdisplay = Xdisplay(:,color_coding{i}(1)+1:size(Xdisplay,2)-1);
            REGdisplay = REGdisplay + max(abs(min(REGdisplay)));
            Xdisplay(:,color_coding{i}(1)+1:size(Xdisplay,2)-1) = REGdisplay ./ max(max(REGdisplay));
        end

        imagesc(Xdisplay); colormap('gray'); drawnow
        if Selection(i) == 1
            if position == i
                title('Best Model');
            else
                title('Valid Model')
            end
        else
            title('Invalid Model')
        end
        
        if i == length(Selection)
            make_plot = 0;
        else
            i = i+1;
        end
        
    end
end







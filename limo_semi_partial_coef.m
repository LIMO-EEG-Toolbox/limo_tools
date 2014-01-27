function limo_semi_partial_coef(LIMO)

% This function returns the semi-partial correlation coefficient
% A coef is simply the difference between the full and a reduced model
% See limo_glm for details on computation
% Cyril Pernet 12-06-2009
% Cyril Pernet 12-11-2013 update for bootstrap and factorial designs
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2010

%% simple checking

if LIMO.design.nb_continuous == 0 && length(LIMO.design.nb_conditions) == 1
    errordlg('there is only 1 categorical factor, nothing to compute')
    return
elseif sum(LIMO.design.nb_conditions)==0 && length(LIMO.design.nb_continuous) == 1
    errordlg('there is only 1 continuous factor, nothing to compute')
    return
end

%% move to the right directoty and load data
cd (LIMO.dir);
load Yr; load R2;

% do the computation
% -----------------------
tmp = compute(LIMO,Yr,R2,1);
clear R2

for i=1:size(tmp,3)
    name = sprintf('semi_partial_coef_%g',i);
    semi_partial_coef = squeeze(tmp(:,:,i,:));
    save(name,'semi_partial_coef','-v7.3')
end

% same under H0 and TFCE
% -------------------------
if LIMO.design.bootstrap == 1
    
    % Compute under H0
    % ---------------
    
    % get the boot table used, the R2 under H0
    cd('H0'); load('boot_table'); load('H0_R2');
    H0 = NaN(size(Yr,1),size(Yr,2),size(tmp,3),3,size(boot_table,2)); clear tmp
    
    % compute nboot times
    for b=1:size(boot_table,2)
        fprintf('computing semi parital coef under H0 - bootstrap %g\n',b);
        H0(:,:,:,:,b) = compute(LIMO,Yr(:,:,boot_table(:,b)),squeeze(H0_R2(:,:,1,b)),0);
    end
        
    for i=1:size(H0,3)
        name = sprintf('H0_semi_partial_coef_%g',i);
        H0_semi_partial_coef = squeeze(H0(:,:,i,[2 3],:)); % only keep F and p
        save(name,'H0_semi_partial_coef','-v7.3')
    end
    clear H0 H0_semi_partial_coef
    
    % run TFCE
    % ----------
    
    if LIMO.design.tfce == 1
        
         % check / get TFCE scores for all covariates
        cd (LIMO.dir); names = dir('semi_partial_coef*.mat');
       
        for i=1:size(names,1)
            name = sprintf('semi_partial_coef_%g',i); load(name)
            fprintf('Thresholding observed semi partial coef %g using TFCE \n',i);
            tfce_score = limo_tfce(squeeze(semi_partial_coef(:,:,2)), LIMO.data.neighbouring_matrix);
            cd TFCE; newname = sprintf('tfce_%s',name); save(newname,'tfce_score');
            clear tfce_score; cd ..
            
            % check / get H0 for TFCE
            cd H0; disp('Thresholding H0 semi partial coef using TFCE');
            newname = sprintf('H0_%s',name); load(newname)
            tfce_H0_score = limo_tfce(squeeze(H0_semi_partial_coef(:,:,i,:)), LIMO.data.neighbouring_matrix);
            newname = sprintf('tfce_%s',newname); save(newname,'tfce_H0_score');
            clear tfce_H0_score; cd ..
        end
    end
end

clear Yr
disp('semi-partial coefficent analysis done ...')


end

%% Compute
function Partial_coef = compute(LIMO,Yr,R2,flag)
% LIMO is the LIMO.mat structure
% Yr ther data electrode * time frame * trials
% R2 is the R2 data
% flag is to indicate to print what's going on or not

N = 1:size(LIMO.design.X,2);
df = LIMO.model.model_df(1);

% set the columns of the reduced designs
% --------------------------------------
index = 1;

% categorical variables
for i=1:length(LIMO.design.nb_conditions)
    if i == 1
        effect = index:LIMO.design.nb_conditions(i);
    else
        effect = index:index+LIMO.design.nb_conditions(i)-1;
    end
    regressors{i} = setdiff(N,effect);
    index = index+LIMO.design.nb_conditions(i);
end

% interactions 
if length(LIMO.design.nb_interactions)>1 && LIMO.design.nb_interactions(1) ~= 0
    for j=1:length(LIMO.design.nb_interactions)
        if j == 1
            effect = index:LIMO.design.nb_interactions(i);
        else
            effect = index:index+LIMO.design.nb_interactions(i)-1;
        end
        regressors{i+j} = setdiff(N,effect);
        index = index+LIMO.design.nb_interactions(i);
    end
end
 
% continuous variables
if LIMO.design.nb_continuous~=0
    s = length(regressors);
    for k=1:length(LIMO.design.nb_continuous)
        regressors{s+k} = setdiff(N,index);
        index = index+1;
    end
end

% compute reduced R2 and do the stat
% ----------------------------------

Partial_coef = NaN(size(Yr,1),size(Yr,2),length(regressors),3);


for r=1:length(regressors)
    
    if flag ~=0
        if r<length(LIMO.design.nb_conditions)+length(LIMO.design.nb_interactions)
            fprintf('computing coef. for the categorical variable %g ,',r);
        else
            fprintf('computing coef. for the continuous variable %g ,',r-length(LIMO.design.nb_conditions)-length(LIMO.design.nb_interactions));
        end
    end
    
    X = LIMO.design.X(:,regressors{r});
    df_reduced = rank(X)-1;
    dfe = (df-df_reduced);
    
    % fit the model
    for i = 1:size(Yr,1)
        if flag ~= 0
            fprintf('...electrode %g \n',i);
        end
        Y = squeeze(Yr(i,:,:))';
        T = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        R = eye(size(Y,1)) - (X*pinv(X));
        Betas = pinv(X)*Y;
        C = eye(size(X,2));
        C(:,size(X,2)) = 0;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0;
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;
        H    = (Betas'*X'*M*X*Betas);
        R2_reduced = diag(H)./diag(T);      % dim time frames 
        Partial_coef(i,:,r,1)= squeeze(R2(i,:,1)) -  R2_reduced';  % dim electrode i * frame : * regressor r -->R2
        Partial_coef(i,:,r,2) = ((size(Yr,3)-df-1).*squeeze(Partial_coef(i,:,r,1))) ./ ((df-df_reduced).*(1-squeeze(R2(i,:,1)))); % --> F
        Partial_coef(i,:,r,3) = 1 - fcdf(squeeze(Partial_coef(i,:,r,2)), df, dfe); % --> p
    end
end
end

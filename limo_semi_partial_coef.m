function limo_semi_partial_coef(varargin)

% This function returns the semi-partial correlation coefficient
% A coef is simply the difference between the full and a reduced model
%
% FORMAT limo_semi_partial_coef(LIMO)
%        limo_semi_partial_coef(data_path,analysis type,observations,X,W,nb_conditions,nb_interactions,nb_continuous,method,model df)
%
% INPUTS LIMO is the LIMO.mat structure after running limo_glm
%        data_path       = where to find the data and R2
%        analysis type   = 'Time', 'Frequency' or 'Time-Frequency'
%        observations    - 'Channels' or 'IC'
%        X               = the design matrix
%        W               = weights to apply to X and Y
%        nb_conditions   = a vector indicating the number of conditions per factor
%        nb_interactions = a vector indicating number of columns per interactions
%        nb_continuous   = number of covariates
%        method          = 'OLS', 'WLS', 'IRLS'
%        model df        = [df dfe]
%
% See limo_glm for details on computation
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2019

%% chech data input
if nargin == 1
    LIMO                        = varargin{1};
else
    LIMO.dir                    = varargin{1};
    LIMO.Analysis               = varargin{2};
    LIMO.Type                   = varargin{3};
    LIMO.design.X               = varargin{4};
    LIMO.design.weights         = varargin{5};
    LIMO.design.nb_conditions   = varargin{6};
    LIMO.design.nb_interactions = varargin{7};
    LIMO.design.nb_continuous   = varargin{8};
    LIMO.design.method          = varargin{9};
    LIMO.model.model_df         = varargin{10};
end
clear varargin

%% simple checking

if LIMO.design.nb_continuous == 0 && length(LIMO.design.nb_conditions) == 1
    errordlg('there is only 1 categorical factor, nothing to compute')
    return
elseif sum(LIMO.design.nb_conditions)==0 && LIMO.design.nb_continuous == 0
    errordlg('there is only 1 continuous factor, nothing to compute')
    return
end

%% move to the right directoty and load data
disp('------------------------------------------')
disp('Starting semi_partial_coefficient analyses')
disp('------------------------------------------')

cd (LIMO.dir);

% do the computation
% -----------------------
if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'Frequency')
    load Yr; load R2;
    tmp = compute(LIMO,Yr,R2,1);
    clear R2
    
    for i=1:size(tmp,3)
        name = sprintf('semi_partial_coef_%g.mat',i);
        semi_partial_coef = squeeze(tmp(:,:,i,:));
        save(name,'semi_partial_coef')
    end
    
elseif strcmp(LIMO.Analysis,'Time-Frequency')
    
    disp('loading 4D data - be patient')
    load Yr; load R2;
    Yr = limo_tf_4d_reshape(Yr);
    R2 = limo_tf_4d_reshape(R2);
    tmp = compute(LIMO,Yr,R2,1);
    clear R2
    
    for i=1:size(tmp,3)
        name = sprintf('semi_partial_coef_%g.mat',i);
        tmp2 = squeeze(tmp(:,:,i,:));
        semi_partial_coef = limo_tf_4d_reshape(tmp2);
        save(name,'semi_partial_coef','-v7.3')
    end
    
end

% same under H0 and TFCE
% -------------------------
if isfield(LIMO.design,'boostrap')
    if LIMO.design.bootstrap == 1
        disp('loading data under H0')
        
        % Compute under H0
        % ---------------
        
        % get the boot table used, the R2 under H0
        cd('H0'); load('boot_table'); load('H0_R2');
        H0 = NaN(size(Yr,1),size(Yr,2),size(tmp,3),3,size(boot_table,2)); clear tmp
        
        if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'Frequency')
            % compute nboot times
            for b=1:size(boot_table,2)
                fprintf('computing semi parital coef under H0 - bootstrap %g\n',b);
                H0(:,:,:,:,b) = compute(LIMO,Yr(:,:,boot_table(:,b)),squeeze(H0_R2(:,:,1,b)),0);
            end
            
            for i=1:size(H0,3)
                name = sprintf('H0_semi_partial_coef_%g.mat',i);
                H0_semi_partial_coef = squeeze(H0(:,:,i,[2 3],:)); % only keep F and p
                save(name,'H0_semi_partial_coef','-v7.3')
            end
            clear H0 H0_semi_partial_coef
            
        elseif strcmp(LIMO.Analysis,'Time-Frequency')
            % compute nboot times
            for b=1:size(boot_table,2)
                fprintf('computing semi parital coef under H0 - bootstrap %g\n',b);
                staked_R2 = limo_tf_4d_reshape(squeeze(H0_R2(:,:,:,:,b)));
                H0(:,:,:,:,b) = compute(LIMO,Yr(:,:,boot_table(:,b)),squeeze(staked_R2(:,:,1)),0);
            end
            
            for i=1:size(H0,3)
                name = sprintf('H0_semi_partial_coef_%g.mat',i);
                tmp = squeeze(H0(:,:,i,[2 3],:)); % only keep F and p
                H0_semi_partial_coef = limo_tf_5d_reshape(tmp);
                save(name,'H0_semi_partial_coef','-v7.3')
            end
            clear H0 H0_semi_partial_coef
            
        end
        
        % run TFCE
        % ----------
        
        if LIMO.design.tfce == 1
            
            % check / get TFCE scores for all - transform F values
            cd (LIMO.dir); names = dir('semi_partial_coef*.mat');
            
            for i=1:size(names,1)
                name = sprintf('semi_partial_coef_%g',i); load(name)
                fprintf('Thresholding observed semi partial coef %g using TFCE \n',i);
                if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'Frequency')
                    tfce_score = limo_tfce(2,squeeze(semi_partial_coef(:,:,2)), LIMO.data.neighbouring_matrix);
                else
                    tfce_score = limo_tfce(3,squeeze(semi_partial_coef(:,:,:,2)), LIMO.data.neighbouring_matrix);
                end
                cd TFCE; newname = sprintf('tfce_%s',name); save(newname,'tfce_score');
                clear tfce_score; cd ..
                
                % check / get H0 for TFCE
                cd H0; disp('Thresholding H0 semi partial coef using TFCE');
                newname = sprintf('H0_%s',name); load(newname)
                if strcmp(LIMO.Analysis,'Time') || strcmp(LIMO.Analysis,'Frequency')
                    tfce_H0_score = limo_tfce(2,squeeze(H0_semi_partial_coef(:,:,1,:)), LIMO.data.neighbouring_matrix);
                else
                    tfce_H0_score = limo_tfce(3,squeeze(H0_semi_partial_coef(:,:,:,1,:)), LIMO.data.neighbouring_matrix);
                end
                newname = sprintf('tfce_%s',newname); save(newname,'tfce_H0_score');
                clear tfce_H0_score; cd ..
            end
        end
    end
end

clear Yr
disp('semi-partial coefficent analysis done ...')
cd(LIMO.dir)

end

%% Compute
function Partial_coef = compute(LIMO,Yr,R2,flag)
% LIMO is the LIMO.mat structure
% Yr the data electrode * [freq/time] * trials
% R2 is the R2 data
% flag is to indicate to print what's going on or not

N  = 1:size(LIMO.design.X,2);

% set the columns of the reduced designs
% --------------------------------------

index      = 1;
% categorical variables
for i=1:length(LIMO.design.nb_conditions)
    if i == 1
        effect = index:LIMO.design.nb_conditions(i);
    else
        effect = index:index+LIMO.design.nb_conditions(i)-1;
    end
    
    if ~isempty(effect)
        regressors{i} = setdiff(N,effect); % all but the effect of interest
        index = index+LIMO.design.nb_conditions(i);
    end
end

% interactions
if LIMO.design.nb_interactions ~= 0
    for j=1:length(LIMO.design.nb_interactions)
        if j == 1
            effect = index:LIMO.design.nb_interactions(i);
        else
            effect = index:index+LIMO.design.nb_interactions(i)-1;
        end
        
        if ~isempty(effect)
            regressors{i+j} = setdiff(N,effect);
            index = index+LIMO.design.nb_interactions(i);
        end
    end
end

% continuous variables
if LIMO.design.nb_continuous~=0
    if exist('regressors','var')
        s = length(regressors);
    else
        s = 0;
    end
    
    for k=1:LIMO.design.nb_continuous
        regressors{s+k} = setdiff(N,index);
        index = index+1;
    end
end


% compute reduced R2 and do the stat
% ----------------------------------
Partial_coef = NaN(size(Yr,1),size(Yr,2),length(regressors),3);
for r=1:length(regressors)
    
    if flag ~=0
        fprintf('computing coef. for variable %g ,',r);
    end
    
    X  = LIMO.design.X(:,regressors{r});
    if strcmp(LIMO.design.method,'OLS') || strcmp(LIMO.design.method,'WLS')
        
        for i = 1:size(Yr,1)
            if flag ~= 0
                fprintf('...%s %g \n',LIMO.Type,i);
            end
            Y                     = squeeze(Yr(i,:,:))';
            T                     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
            if strcmp(LIMO.design.method,'OLS')
                Betas             = pinv(X)*Y;
                W                 = ones(size(Y,1),1);
            else
                [Betas,W]         = limo_WLS(X,Y);
            end
            WX                    = [X(:,1:end-1).*repmat(W,1,size(X,2)-1) X(:,end)];
            R                     = eye(size(Y,1)) - (WX*pinv(WX));
            C                     = eye(size(X,2));
            C(:,size(X,2))        = 0;
            C0                    = eye(size(X,2)) - C*pinv(C);
            X0                    = WX*C0;
            R0                    = eye(size(Y,1)) - (X0*pinv(X0));
            M                     = R0 - R;
            H                     = (Betas'*X'*M*X*Betas);
            R2_reduced            = diag(H)./diag(T);  % dim [freq/time]
            Partial_coef(i,:,r,1) = abs(squeeze(R2(i,:,1))- R2_reduced');  % dim electrode i * [freq/time] * regressor r -->R2
            df                    = rank([LIMO.design.X(:,1:end-1).*repmat(LIMO.design.weights(i,:)',1,size(LIMO.design.X,2)-1) LIMO.design.X(:,end)])-1;
            df_reduced            = rank(WX)-1;
            dfe                   = (df-df_reduced);
            Partial_coef(i,:,r,2) = ((size(Yr,3)-df-1).*squeeze(Partial_coef(i,:,r,1))) ./ ((df-df_reduced).*(1-squeeze(R2(i,:,1)))); % --> F
            if strcmp(LIMO.design.method,'OLS') 
                Partial_coef(i,:,r,3) = 1 - fcdf(squeeze(Partial_coef(i,:,r,2)), df, dfe); % --> p
            end
        end
        
    else % IRLS

        for i = 1:size(Yr,1)
            if flag ~= 0
                fprintf('...%s %g \n',LIMO.Type,i);
            end
            % got to iterate for each cell
            for j=1:size(Yr,2)
                Y                     = squeeze(Yr(i,j,:))';
                T                     = (Y-mean(Y))'*(Y-mean(Y));
                WX                    = [X(:,1:end-1).*repmat(LIMO.design.weights(i,j,:),1,size(X,2)-1) X(:,end)];
                R                     = eye(size(Y,1)) - (WX*pinv(WX));
                C                     = eye(size(X,2));
                C(:,size(X,2))        = 0;
                C0                    = eye(size(X,2)) - C*pinv(C);
                X0                    = WX*C0;
                R0                    = eye(size(Y,1)) - (X0*pinv(X0));
                M                     = R0 - R;
                Betas                 = pinv(WX)*(Y.*LIMO.design.weights(i,j,:)');
                H                     = (Betas'*X'*M*X*Betas);
                R2_reduced            = H./T;  % dim [freq/time]
                Partial_coef(i,j,r,1) = squeeze(R2(i,j,1)) -  R2_reduced';  % dim electrode i * [freq/time] * regressor r -->R2
                
                df_reduced            = rank(WX)-1;
                dfe                   = (df-df_reduced);
                Partial_coef(i,j,r,2) = ((size(Yr,3)-df-1).*squeeze(Partial_coef(i,j,r,1))) ./ ((df-df_reduced).*(1-squeeze(R2(i,j,1)))); % --> F
                if strcmp(LIMO.design.method,'OLS')
                    Partial_coef(i,j,r,3) = 1 - fcdf(squeeze(Partial_coef(i,j,r,2)), df, dfe); % --> p
                end
            end
        end
    end
end
end
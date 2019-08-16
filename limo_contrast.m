function result = limo_contrast(varargin)

% LIMO_CONTRAST - computes contrasts using outputs from LIMO_GLM
% This function uses the parameters computed with LIMO_glm.m
% simply read the matrix format and perform a series of contrasts
% this can be performed automatically by entering the contrast to perform
% in C or if C is not specified via a GUI (limo_contrast_manager)
% the function used here is similar as the one used in limo_glm1.m
%
% FORMAT:
% result = limo_contrast(Y, Betas, LIMO, contrast type, analysis type)
%
% INPUT:
% Y              = 3D data
% Betas          = betas computed in limo_glm1 
% LIMO           = the LIMO.mat with the design matrix and contrast
% contrast type  = 0 for T test, 1 for F test
% analysis type  = 1 for 1st level and 2nd level ANOVA/ANCOVA 
%                  2 for 1/2nd level bootrapped ANOVA/ANCOVA
%
% FORMAT:
% result = limo_contrast(Yr,LIMO,3);
%
% INPUT:
% Y              = 3D data
% LIMO           = the LIMO.mat with the design matrix and contrast
% analysis type  = 3 for 2nd level repeated measures ANOVA 
%                  4 for 2nd level bootrapped repeated measures ANOVA
%
% OUTPUT
% con/ess maps saved on disk
% result is the multivariate output (not in use for now)
%
% *****************************************************
% See also limo_glm1, limo_results, limo_contrast_manager
%
% Cyril Pernet v4 26/09/2010
% updated 21-16-2013
% -----------------------------
%  Copyright (C) LIMO Team 2010


%% nargin stuff

type = varargin{end};



%% Analyses

switch type
    
    case{1}
        % -----------------------------------------------------------------
        % Contrast for 1st level analyses and 2nd level regression/ANOVA/ANCOVA
        % -----------------------------------------------------------------
        
        Y           = varargin{1}; % 3D data
        Betas       = varargin{2}; % 3D betas
        LIMO        = varargin{3};
        X           = LIMO.design.X;
        nb_beta     = size(LIMO.design.X,2);
        contrast_nb = size(varargin{3}.contrast,2);
        C           = varargin{3}.contrast{size(varargin{3}.contrast,2)}.C;
        Method      = varargin{3}.design.type_of_analysis;
        try
            dfe     = varargin{3}.model.model_df(2);
        catch ME
            dfe     = size(Y,1)-rank(X); %% hapens for 2nd level N-way ANOVA or ANCOVA
        end
        Test        = varargin{4};
        
 
        % compute Projection onto the error
        load Res; % rather than projecting Y onto error use Res because Res depends on how the GLM was done (OLS,WLS, IRLS)

        if strcmp(LIMO.Analysis ,'Time-Frequency')
            Y = limo_tf_4d_reshape(Y);
            Betas = limo_tf_4d_reshape(Betas);
            Res = limo_tf_4d_reshape(Res);
        end
                
        if strcmp(Method,'Mass-univariate')
            
            % create con or ess file
            if Test == 0
                con = NaN(size(Y,1),size(Y,2),5); % dim 3 = C*Beta/se/df/t/p
                filename = sprintf('con_%g.mat',size(varargin{3}.contrast,2));
            else
                ess = NaN(size(Y,1),size(Y,2),size(C,1)+4); % dim 3 = C*Beta/se/df/F/p
                filename = sprintf('ess_%g.mat',size(varargin{3}.contrast,2));
            end
            
            % update con/ess file
            array = find(~isnan(Y(:,1,1))); % skip empty electrodes
            for e = 1:length(array)
                electrode = array(e); warning off;
                if strcmp(LIMO.Type,'Channels')
                    fprintf('analyzing channel %g/%g \n',e,size(array,1));
                else
                    fprintf('analyzing component %g/%g \n',e,size(array,1));
                end
                                
                % T contrast
                % -----------
                if Test == 0
                    
                    var   = (squeeze(Res(electrode,:,:))*squeeze(Res(electrode,:,:))') / dfe;
                    % Update con file [mean value, se, df, t, p]
                    con(electrode,:,1) = C*squeeze(Betas(electrode,:,:))' ;  
                    con(electrode,:,2) = sqrt(diag(var)'.*(C*pinv(X'*X)*C')); 
                    con(electrode,:,3) = dfe; 
                    con(electrode,:,4) = (C*squeeze(Betas(electrode,:,:))') ./ sqrt(diag(var)'.*(C*pinv(X'*X)*C')); 
                    con(electrode,:,5) = 1-tcdf(squeeze(con(electrode,:,4)), dfe); 
                    
                    % F contrast
                    % ----------
                else
                    % Update ess file [mean value, se, df, F, p]
                    ess(electrode,:,1:size(C,1)) = (C*squeeze(Betas(electrode,:,:))')' ; % contrast
                    
                    R = eye(size(Y,3)) - (X*pinv(X));
                    E = (squeeze(Res(electrode,:,:))*squeeze(Res(electrode,:,:))');
                    c = zeros(length(C));
                    for n=1:length(C)
                        c(n,n) = C(n);
                    end
                    
                    try
                        C0 = eye(rank(X)+1) - c*pinv(c);
                    catch ME
                        C0 = eye(rank(X)) - c*pinv(c);
                    end
                    X0 = X*C0;
                    R0 = eye(size(Y,3)) - (X0*pinv(X0));
                    M = R0 - R;
                    H = (squeeze(Betas(electrode,:,:))*X'*M*X*squeeze(Betas(electrode,:,:))');
                    if rank(c) == 1
                        df = 1;
                    else
                        df = rank(c) - 1;
                    end
                    ess(electrode,:,end-3) = diag(E)/dfe;  
                    ess(electrode,:,end-2) = df;   
                    ess(electrode,:,end-1) = (diag(H)/df)./(diag(E)/dfe);  % F value
                    ess(electrode,:,end)   = 1 - fcdf(ess(electrode,:,end-1), rank(c)-1, dfe);   % p value
                end
                
            end
            result = [];
            
            if Test == 0
                if strcmp(LIMO.Analysis ,'Time-Frequency')
                    con = limo_tf_4d_reshape(con);
                end
                save ([filename], 'con'); clear con
            else
                if strcmp(LIMO.Analysis ,'Time-Frequency')
                    ess = limo_tf_4d_reshape(ess);
                end
                save ([filename], 'ess'); clear ess
            end
            
            
        elseif strcmp(Method,'Multivariate')
            % ------------------------------
            
            con = NaN(size(Y,2),2); %  F /p values (always the same no matter Roy or Pillai)
            for time = 1:size(Y,2)
                fprintf('time frame %g \n',time);
                
                E = (squeeze(Y(:,time,:))*R*squeeze(Y(:,time,:))');
                c = zeros(length(C));
                for n=1:length(C)
                    c(n,n) = C(n);
                end
                
                try
                    C0 = eye(rank(X)+1) - c*pinv(c);
                catch ME
                    C0 = eye(rank(X)) - c*pinv(c);
                end
                X0 = X*C0;
                R0 = eye(size(Y,2)) - (X0*pinv(X0));
                M = R0 - R;
                H = (squeeze(Betas(:,time,:))*X'*M*X*squeeze(Betas(:,time,:))');
                
                multivariate.EV    = limo_decomp(E,H);
                multivariate.theta = max(multivariate.EV) / (1+max(multivariate.EV));
                multivariate.V     = sum(multivariate.EV ./ (1+multivariate.EV));
                multivariate.df    = size(Y,2);
                multivariate.dfe   = abs(size(Y,1) - (nb_beta-1) - (multivariate.df-1));
                multivariate.T_contrast    = sqrt((dfe*max(multivariate.EV))/multivariate.df);
                multivariate.pval_contrast = 1-fcdf(multivariate.T_contrast, multivariate.df, abs(dfe));
                result = multivariate;
            end
        end
                
    case{2} 
        % -----------------------------------------------------------------
        % ANOVA/ANCOVA/Regression bootstraps
        % -----------------------------------------------------------------
        
        % INPUTS
        % ------
        y           = varargin{1}; % 3D original data
        Betas       = varargin{2}; % 4D Betas under H0
        nboot       = size(Betas,4);
        LIMO        = varargin{3}; % LIMO
        X           = LIMO.design.X; % design
        nb_beta     = size(LIMO.design.X,2);
        contrast_nb = size(varargin{3}.contrast,2);
        C           = varargin{3}.contrast{size(varargin{3}.contrast,2)}.C;
        Method      = varargin{3}.design.type_of_analysis;
        try
            dfe     = varargin{3}.model.model_df(2);
        catch ME
            dfe     = size(Y,1)-rank(X); %% hapens for 2nd level N-way ANOVA or ANCOVA
        end
        Test        = varargin{4};
        
      
        % make data files
        % ----------------
        if Test == 0
            H0_con = NaN(size(y,1),size(y,2),3,nboot); % dim 3 = C*Beta/t/p
            filename = sprintf('H0_con_%g.mat',size(varargin{3}.contrast,2));
        else
            H0_ess = NaN(size(y,1),size(y,2),size(C,1)+2,nboot); % dim 3 = C*Beta/F/p
            filename = sprintf('H0_ess_%g.mat',size(varargin{3}.contrast,2));
        end

        
        % prepare data for bootstrap
        % --------------------------
        % if categorical design, center data 1st
        % ---------------------------------------
        if LIMO.design.nb_continuous == 0
            for e=1:size(y,1)
                centered_y = NaN(size(y,1),size(y,2),size(y,3));
                if LIMO.design.nb_interactions ~=0
                    % look up the last interaction to get unique groups
                    if length(LIMO.design.nb_interactions) == 1
                        start_at = sum(LIMO.design.nb_conditions);
                    else
                        start_at = sum(LIMO.design.nb_conditions)+sum(LIMO.design.nb_interactions(1:end-1));
                    end
                    
                    for cel=(start_at+1):(start_at+LIMO.design.nb_interactions(end))
                        index = find(X(:,cel));
                        centered_y(e,:,index) = squeeze(y(e,:,index)) - repmat(mean(squeeze(y(e,:,index)),2),1,length(index));
                    end
                    
                elseif size(LIMO.design.nb_conditions,2) == 1
                    % no interactions because just 1 factor
                    for cel=1:LIMO.design.nb_conditions
                        index = find(X(:,cel));
                        centered_y(e,:,index) = squeeze(y(e,:,index)) - repmat(nanmean(squeeze(y(e,:,index)),2),1,length(index));
                    end
                    
                else
                    % create fake interaction to get groups
                    [tmpX interactions] = make_interactions(X, LIMO.design.nb_conditions);
                    if length(interactions) == 1
                        start_at = sum(LIMO.design.nb_conditions);
                    else
                        start_at = sum(LIMO.design.nb_conditions)+sum(LIMO.design.interactions(1:end-1));
                    end
                    
                    for cel=(start_at+1):(start_at+interactions(end))
                        index = find(X(:,cel));
                        centered_y(e,:,index) = squeeze(y(e,:,index)) - repmat(mean(squeeze(y(e,:,index)),2),1,[size(y(index,:),1)]);
                    end
                end
            end
        end
        
        
        % start the analysis
        % -------------------
        load boot_table
        array = find(~isnan(y(:,1,1))); % skip empty electrodes
        design = X;
        
        if strcmp(Method,'Mass-univariate')
            % ---------------------------------
            for e = 1:length(array)
                electrode = array(e); warning off;
                fprintf('compute bootstrap electrode %g ... \n',electrode)
                for B = 1:nboot
                    if ~iscell(boot_table)
                        resampling_index = boot_table(:,B); % 1st level boot_table all the same ever
                    else
                        resampling_index = boot_table{electrode}(:,B);
                    end
                    
                    % create data under H0
                    if LIMO.design.nb_continuous == 0
                        % sample from the centered data in categorical designs
                        Y = squeeze(centered_y(electrode,:,resampling_index))';
                        X = design(resampling_index,:); % resample X as Y
                    else
                        % sample and break the link between Y and (regression and AnCOVA designs)
                        Y = squeeze(y(electrode,:,resampling_index))';
                        X = design(find(~isnan(y(electrode,1,:))),:);
                        if LIMO.design.zscore == 1 % rezscore the covariates
                            N = LIMO.design.nb_conditions + LIMO.design.nb_interactions;
                            if N==0
                                if sum(mean(X(:,1:end-1),1)) > 10e-15
                                    X(:,1:end-1) = zscore(X(:,1:end-1));
                                end
                            else
                                if sum(mean(X(:,N+1:end-1),1)) > 10e-15
                                    X(:,N+1:end-1) = zscore(X(:,N+1:end-1));
                                end
                            end
                        end
                    end
                    
                    % compute Projection onto the error
                    R = eye(size(Y,1)) - (X*pinv(X));
                    
                    % T contrast
                    % -----------
                    if Test == 0
                        
                        var   = ((R*Y)'*(R*Y)) / dfe; % error of H0 data
                        H0_con(electrode,:,1,B) = C*squeeze(Betas(electrode,:,:,B))' ;  % contrast using betas H0
                        H0_con(electrode,:,2,B) = (C*squeeze(Betas(electrode,:,:,B))') ./ sqrt(diag(var)'.*(C*pinv(X'*X)*C')); % T value
                        H0_con(electrode,:,3,B) = 1-tcdf(squeeze(H0_con(electrode,:,2,B)), dfe); % p value
                        
                        % F contrast
                        % ----------
                    else
                        H0_ess(electrode,:,1:size(C,1),B) = (C*squeeze(Betas(electrode,:,:,B))')' ; % contrast
                        
                        E = (Y'*R*Y);
                        c = zeros(length(C));
                        for n=1:length(C)
                            c(n,n) = C(n);
                        end
                        
                        try
                            C0 = eye(rank(X)+1) - c*pinv(c);
                        catch ME
                            C0 = eye(rank(X)) - c*pinv(c);
                        end
                        X0 = X*C0;
                        R0 = eye(size(Y,1)) - (X0*pinv(X0));
                        M = R0 - R;
                        H = (squeeze(Betas(electrode,:,:,B))*X'*M*X*squeeze(Betas(electrode,:,:,B))');
                        if rank(c) == 1
                            df = 1;
                        else
                            df = rank(c) - 1;
                        end
                        H0_ess(electrode,:,end-1,B)    = (diag(H)/df)./(diag(E)/dfe);  % F value
                        H0_ess(electrode,:,end,B) = 1 - fcdf(H0_ess(electrode,:,end-1,B), rank(c)-1, dfe);   % p value
                    end
                end
            end
            
            result = [];
            if Test == 0;
                save ([filename], 'H0_con'); clear H0_con
            else
                save ([filename], 'H0_ess'); clear H0_ess
            end
        end
        
        % ----------------------------------------
        if strcmp(Method,'Multivariate')
            % ----------------------------------------
            
            for e = 1:size(y,1)
                electrode = array(e); warning off;
                fprintf('compute bootstrap electrode %g ... \n',electrode)
                for B = 1:nboot
                    % create data under H0
                    if LIMO.design.nb_continuous == 0
                        % sample from the centered data in categorical designs
                        Y = centered_y(boot_table(:,B));
                        X = design(boot_table(:,B)); % resample X as Y
                    else
                        % sample and break the link between Y and (regression and AnCOVA designs)
                        Y = y(boot_table(:,B));
                        if LIMO.design.zscore == 1 % rezscore the covariates
                            N = LIMO.design.nb_conditions + LIMO.design.nb_interactions;
                            if sum(mean(X(:,N+1:end-1),1)) ~= 0
                                X(:,N+1:end-1) = zscore(X(:,N+1:end-1));
                            end
                        end
                    end
                    
                    % compute Projection onto the error
                    R = eye(size(Y,1)) - (X*pinv(X));
                    
                    E = (Y'*R*Y);
                    c = zeros(length(C));
                    for n=1:length(C)
                        c(n,n) = C(n);
                    end
                    
                    try
                        C0 = eye(rank(X)+1) - c*pinv(c);
                    catch ME
                        C0 = eye(rank(X)) - c*pinv(c);
                    end
                    X0 = X*C0;
                    R0 = eye(size(Y,1)) - (X0*pinv(X0));
                    M = R0 - R;
                    H = (Betas'*X'*M*X*Betas);
                    
                    multivariate.EV    = limo_decomp(E,H);
                    multivariate.theta = max(multivariate.EV) / (1+max(multivariate.EV));
                    multivariate.V     = sum(multivariate.EV ./ (1+multivariate.EV));
                    multivariate.df    = size(Y,2);
                    multivariate.dfe   = abs(size(Y,1) - (nb_beta-1) - (multivariate.df-1));
                    multivariate.T_contrast    = sqrt((dfe*max(multivariate.EV))/multivariate.df);
                    multivariate.pval_contrast = 1-fcdf(multivariate.T_contrast, multivariate.df, abs(dfe));
                    result = multivariate;
                end
            end
        end
        
    case(3)
        % --------------------------------------------
        %              Repeated Measure ANOVA
        % ---------------------------------------------
        
        Yr = varargin{1};
        LIMO = varargin{2};
        if LIMO.Level == 1
            error('1st level Analysis detected - limo_contrast line 434 wrong case');
        end
        gp_values = LIMO.design.nb_conditions;
        index = size(LIMO.contrast,2);
        C = LIMO.contrast{index}.C;
        
        % [mean value, se, df, F, p])
        if gp_values == 1
            ess = zeros(size(Yr,1),size(Yr,2),5); 
            for electrode = 1:size(Yr,1)
                fprintf('electrode %g \n',electrode); 
                % Inputs
                tmp = squeeze(Yr(electrode,:,:,:));
                Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                gp = LIMO.data.Cat(find(~isnan(tmp(1,:,1))),:);
                % mean, se, df
                n = size(Y,2);
                g=floor((20/100)*n); 
                for time=1:size(Y,1)
                    ess(electrode,time,1) = nanmean(C(1:size(Y,3))*squeeze(Y(time,:,:))',2);
                    ess(electrode,time,2) = sqrt(C(1:size(Y,3))*cov(squeeze(Y(time,:,:)))*C(1:size(Y,3))');
                end
                df  = rank(C); dfe = n-df;
                ess(electrode,:,3) = dfe;
                % F and p
                result = limo_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(Y,3)));
                ess(electrode,:,4) = result.F;
                ess(electrode,:,5) = result.p;
            end
        else
            ess = zeros(size(Yr,1),size(Yr,2),5); % dim rep measures, F,p
            ess2 = zeros(size(Yr,1),size(Yr,2),5); % dim gp*interaction F,p
            % design matrix for gp effects
            k = LIMO.design.nb_conditions;
            gp_vector = LIMO.data.Cat;
            gp_values = unique(gp_vector); k = length(gp_values); X = NaN(size(gp_vector,1),k+1);
            for g =1:k; X(:,g) = gp_vector == gp_values(g); end; X(:,end) = 1; % design matrix for gp effects

            % call rep anova
            for electrode = 1:size(Yr,1)
                fprintf('electrode %g \n',electrode); 
                % Inputs
                tmp = squeeze(Yr(electrode,:,:,:));
                Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                gp = LIMO.data.Cat(find(~isnan(tmp(1,:,1))),:);
                XB = X(find(~isnan(tmp(1,:,1))),:);
                % mean, se, df
                n = size(Y,2);
                g=floor((20/100)*n); 
                for time=1:size(Y,1)
                    [v,indices] = sort(squeeze(Y(time,:,:))); % sorted data
                    TD(time,:,:) = v((g+1):(n-g),:); % trimmed data 
                    ess(electrode,time,1) = nanmean(C(1:size(TD,3))*squeeze(TD(time,:,:))',2);
                    I = zeros(1,1,n); I(1,1,:) = (C(1:size(TD,3))*squeeze(Y(time,:,:))')'; % interaction
                    ess2(electrode,time,1) = limo_trimmed_mean(I);
                    v(1:g+1,:)=repmat(v(g+1,:),g+1,1);
                    v(n-g:end,:)=repmat(v(n-g,:),g+1,1); % winsorized data
                    [~,reorder] = sort(indices);
                    for j = 1:size(Y,3), SD(:,j) = v(reorder(:,j),j); end % restore the order of original data
                    S(time,:,:) = cov(SD); % winsorized covariance
                    ess(electrode,time,2) = sqrt(C(1:size(TD,3))*squeeze(S(time,:,:))*C(1:size(TD,3))');
                    ess2(electrode,time,2) = NaN;
                end
                df  = rank(C); dfe = n-df;
                ess(electrode,:,3) = dfe;
                % F and p values
                result = limo_rep_anova(Y, gp, LIMO.design.repeated_measure, C(1:size(TD,3)),XB);
                ess(electrode,:,1,4) = result.repeated_measure.F;
                ess(electrode,:,1,5) = result.repeated_measure.p;
                ess2(electrode,:,2,4) = result.interaction.F;
                ess2(electrode,:,2,5) = result.interaction.p;
            end
        end
        
        filename = sprintf('ess_repeated_measure_%g.mat',index);
        save ([filename], 'ess');
        if exist('ess2','var')
            ess = ess2; filename = sprintf('ess_interaction_gp_repeated_measure_%g.mat',index);
            save ([filename], 'ess2');
        end
end




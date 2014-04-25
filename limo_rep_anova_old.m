function results = limo_rep_anova_old(Y,gp,factors,flag)

% This function allows to compute a repeated measure ANOVA (up to 3
% factors) with several groups of subjects (same sample sizes)
% Note that this function is not used by LIMO EEG as repeated measures are
% analyzed using a Hotteling t-test which accounts for sphericity.
%
% FORMAT
% results = limo_rep_anova(Y,[ ],[2 3],1)
%
% INPUT
% Y:       the data as a 3D matrix (frame,subjects,conditions)
%          This input is like a standard repeated measures table i.e. in
%          rows are the subjects and in columns the conditions organized
%          following the factor structure - the only difference is that dim.
%          1 marks the frames on which to repeat the analysis
% factor:  a vector with the levels of the different factors - for instance
%          [6] means Y has 6 columns for 6 levels and 1 factor
%          [2 3] means Y has 6 columns, the 3 first and 3 last are the 2 levels of factor 1 and within these levels we have 3 levels for factor 2 
%          WARNING Max 3 factors because the way interactions are coded - using projections should solve this
% gp:      a vector of values indicating group belonging (e.g. 1 1 2 2) or [ ] if there is only one group
%          WARNING for now this function only handles groups of same size
%          (same size because the sphericity tool doesn't allows different N)
% flag:    0/1 to deactivate/activate the generation of figures
%
% OUTPUT
% results is a mat file with the design matrix X, the names of the effects
% and corresponding F and p values corrected for sphericity
%
% Different codes are used for ANOVAs with 1 repeated factor vs. multiple
% factors - when only one factor exist one can compute all effects and
% errors on all frames in one pass vs. need to pool some data to
% compute the error when there are many factors, slowing down the
% computation. Separate codes (see type below) therefore reflect this
% special need. Ultimatetly, this code will be rewriten using projection
% matrices which should make easier to compute the various effects and make
% the code faster.
%
% The sphericity correction used is the Huynh-Feldt one - this done using
% the GenCalcHFEps.m function writen by Matthew Nelson
%
% Cyril Pernet v1 April 2010
%              v2 July  2010 - changed the Y structure to accommodate
%              sphericity computation
% -----------------------------------------------------------------
%  Copyright (C) LIMO Team 2010


% ------------------------------
%% Get some information from the input
% -----------------------------------
results = [];
if numel(size(Y)) == 2; 
    tmp = NaN(1,size(Y,1),size(Y,2)); tmp(1,:,:) = Y;
    clear Y; Y = tmp; clear tmp
end

nb_factors    = size(factors,2);
nb_effects    = (2^nb_factors - 1);
nb_conditions = prod(factors);

if length(unique(gp)) == 1
    gp =[];
end

if isempty(gp)
    nb_interactions = 0;
    if nb_factors == 1
        type = 1;
    else
        type = 2;
    end
    nb_subjects   = size(Y,2);
    S = repmat([1:nb_subjects],1,nb_conditions)'; % subject vector
else
    nb_interactions = nb_effects;
    if nb_factors == 1
        type = 3;
    else
        type = 4;
    end
    % need to accomodate different sample sizes
    index = 1;
    for i=min(gp):max(gp)
        nb_subjects(index) = sum(gp==i); index = index+1;
    end
    
    v = [];
    for i=1:length(nb_subjects)
        if i == 1; v = [1:nb_subjects(i)]; m = max(v);s = [v];
        else v = m+[1:nb_subjects(i)]; m=max(v); s = [s,v]; end
    end
    S = repmat(s',nb_conditions,1);
    nb_subjects = sum(nb_subjects);
end

% reshape Y 
for dim=1:size(Y,1)
    tmp = squeeze(Y(dim,:,:));
    data(:,dim) = tmp(:);
end
clear tmp Y; Y = data; clear data;

%% Compute the Huynh-Feldt and Geisser-Greenhouse corrections
% ----------------------------------------------------------
% disp('computing sphericity ...')

% create f matrices per repeated measures (re-used below if several Fs)
dim = size(Y,1);
for f=1:nb_factors
    if f ==1;
        F{f} = kron(eye(factors(f)),ones((nb_subjects*(nb_conditions/factors(f))),1));
        Name{f} = sprintf('Factor %g',f); N(f) = nb_subjects*(nb_conditions/factors(f));
        df_factors{f}   = rank(F{f}) -1;
    else
        dim = dim / factors(f-1);
        vec = kron(eye(factors(f)),ones((dim / factors(f)),1));
        F{f} = repmat(vec,size(Y,1)/size(vec,1),1);
        Name{f} = sprintf('Factor %g',f); N(f) = dim / factors(f);
        df_factors{f}   = rank(F{f}) -1;
    end
end

% use F above to create a matrix of effects (repeated measure)
WInFacs = zeros(size(Y,1),nb_factors);
for i = 1:size(F,2)
    tmp = F{i};
    for  j = 1:size(tmp,2)
        tmp(tmp(:,j)==1)=j;
    end
    WInFacs(:,i) = repmat(tmp(:,1),size(Y,1)/length(tmp(:,1)),1);
end
clear tmp

% create a vector for the groups
if gp == 0
    BTFacs = [];
else
    BTFacs = repmat(gp,nb_conditions,1);
end

% get correction values
for n=1:size(Y,2)
    if factors == 2
        EpsHF{n} = 1;
    else
        [EpsHF{n} EpsList{n} EpsGG{n}]=GenCalcHFEps(Y(:,n),BTFacs,WInFacs,S);
    end
end


%% Conpute ANOVA results
% ----------------------

switch type
    case {1}
        
        % ------------------------------------------------------------------------------------------
        %                   One Sample Repeated measure ANOVA
        %                   Y = XB + E with X = [Factor / Subject]
        % ------------------------------------------------------------------------------------------
        
            
            df  = nb_conditions -1;
            dfe = size(Y,1)  - nb_subjects - df;
            
            % create the design matrix for the different levels
            % ------------------------------------------------
            
            x = kron(eye(nb_conditions),ones(nb_subjects,1));  % effect
            x0 = repmat(eye(nb_subjects),nb_conditions,1); % error
            X = [x x0];
            if flag == 1
                figure('Name','Design matrix'); set(gcf,'Color','w'); imagesc(X);
                colormap('gray');  title('ANOVA model','FontSize',16);xlabel('regressors'); 
                ylabel('scaled model values'); drawnow; 
                go = questdlg('start the analysis?');
                if strcmp(go,'No')
                    return
                end
            end
            
            % Compute the repeated measure effect
            % ------------------------------------
            M     = X*pinv(X'*X)*X';
            R     = eye(size(Y,1)) - M;
            SSe   = diag(Y'*R*Y); % in fact the subject effect
            Betas = pinv(x)*Y;  % compute without cst/subjects
            yhat  = x*Betas;
            SS    = diag((yhat-repmat(mean(yhat),size(yhat,1),1))'*(yhat-repmat(mean(yhat),size(yhat,1),1)));
            
            % get the stats
            % -------------
            F_values = (SS'./(cell2mat(EpsHF)*df)) ./ (SSe'/(cell2mat(EpsHF)*dfe));
            p_values = 1 - fcdf(F_values, cell2mat(EpsHF)*df, cell2mat(EpsHF)*dfe);
            
            % output
            % ------
            results.Design = X;
            results.name = 'Factor 1';
            results.Betas{1} = Betas;
            results.F = F_values;
            results.p = p_values;
            

        
        
    case {2}
        % ------------------------------------------------------------------------------------------
        %                   Repeated measure ANOVA with several factors
        %                   Y = XB + E with X = [Factors Interactions / Subject]
        % ------------------------------------------------------------------------------------------
        
        % create the design matrix for the different factors
        % --------------------------------------------------
        
        Subjects = repmat(eye(nb_subjects),nb_conditions,1);
        
        % make all the interactions between factors
        % compute all possible interactions
        for nb_comb = 2:nb_factors
            combinations{nb_comb-1} = nchoosek([1:nb_factors],nb_comb);
        end
        
        % reorganize as a list
        index = 1; v = f+1;
        for i=1:size(combinations,2)
            selection = single(combinations{i});
            for j=1:size(selection,1)
                interaction{index}= selection(j,:);
                Name{v} = sprintf('Interaction %s',num2str(interaction{index}));
                v = v+1;index= index +1;
            end
        end
        
        % multiply factors to make the interaction submatrices
        for i=1:size(interaction,2)
            f = f+1;
            
            % get matrices to multiply into M
            for j=1:size(interaction{i},2)
                M{j} = F{interaction{i}(j)};
                nb_columns(j) = size(M{j},2);
                tmpN(j) = N(interaction{i}(j));
                tmp_df(j) = df_factors(interaction{i}(j));
            end
            df_factors{f} = prod(cell2mat(tmp_df));
            
            % multiply column-wise
            if size(M,2) == 2 % if only 2 matrices
                
                index = 1;
                [value,matrix1]=min(tmpN);
                [value,matrix2]=max(tmpN);
                
                for k = 1:size(M{matrix2},2) % rotates on the matrix with largest N
                    tmp{index} = M{matrix1} .* repmat(M{matrix2}(:,k),1,size(M{matrix1},2));
                    index = index + 1;
                end
                F{f} = cell2mat(tmp); clear tmp;
                
            else
                
                if i == size(interaction,2) %% all factors
                    F{f} = kron(eye(nb_conditions),ones(nb_subjects,1));
                    
                else
                    index = 1;
                    [value,matrix1]=min(tmpN);
                    [value,matrix2]=max(tmpN);
                    
                    for k = 1:size(M{matrix2},2) % rotates on the matrix with largest N
                        tmp{index} = M{matrix1} .* repmat(M{matrix2}(:,k),1,size(M{matrix1},2));
                        index = index + 1;
                    end
                    tmp2 = cell2mat(tmp); clear tmp;
                    
                    [r,c,v] = find((((interaction{i}-matrix1 == 0) + (interaction{i}-matrix2 == 0)) == 0) .* interaction{i});
                    
                    for l=v
                        index = 1;
                        for k=1:size(M{l},2)
                            tmp{index} = tmp2 .* repmat(M{l}(:,k),1,size(tmp2,2));
                            index = index + 1;
                        end
                        tmp3 = cell2mat(tmp); clear tmp;
                        clear tmp2; tmp2 = tmp3; clear tmp3
                    end
                    F{f} = tmp2; clear tmp2;
                end
            end
            clear M nb_columns tmpN
        end
        
        % assemble the parts
        X = [];
        for c = 1:size(F,2)
            X = [X F{c}];
        end
        
        if flag == 1
            figure('Name','Design matrix'); set(gcf,'Color','w'); imagesc([X Subjects]);
            colormap('gray');  title('ANOVA model','FontSize',16);xlabel('regressors');
            ylabel('scaled model values'); drawnow;
            go = questdlg('start the analysis?');
            if strcmp(go,'No')
                return
            end
        end
                    
        % analyse data
        % -----------------------------
        i = 0;
        
        % SS subject
        Betas = pinv(Subjects)*Y;
        yS    = Subjects*Betas;
        SSs  = diag((yS-repmat(mean(yS),size(yS,1),1))'*(yS-repmat(mean(yS),size(yS,1),1)));
        
        % SS factors and error
        for e=1:nb_effects
            
            % SS effect - since this is the same model everywhere, we do
            % not have to run it per frame
            fprintf('compute effect and error %g out of %g',e, nb_effects); disp(' ');
            Betas = pinv(F{e})*Y;
            Rep_betas{e} = Betas;
            yhat  = F{e}*Betas;
            df = df_factors{e};
            dfe = (nb_subjects-1)*df;
            SS{e} = diag((yhat-repmat(mean(yhat),size(yhat,1),1))'*(yhat-repmat(mean(yhat),size(yhat,1),1)));
            if e > nb_factors
                i = i+1;
                if length(interaction{i}) <= 2
                    main_effects = cell2mat(SS(interaction{i}));
                    SS{e} = SS{e} - sum(main_effects,2); % removes sum of SS main effects
                else
                    main_effects_interactions = cell2mat(SS(1:e-1));
                    SS{e} = SS{e} - sum(main_effects_interactions,2); % removes sum of SS main effects and lower interactions
                end
            end
            
            % SS error - here it depends on Y, ie we cannot run it per frame
            for frame = 1:size(Y,2)
                ye = repmat(Y(:,frame),1,size(F{e},2)).*F{e}; % data into the F{e} format
                avg_number = [1:nb_subjects];
                tmp = ye(avg_number,:);
                for avg = 2:nb_conditions
                    avg_number = avg_number + nb_subjects;
                    tmp = tmp+ye(avg_number,:);
                end
                tmp = tmp./(nb_conditions/size(F{e},2)); % produces the data averaged per level
                ye = tmp(:);
                SSt = (nb_conditions/size(F{e},2)) .* (diag((ye-repmat(mean(ye),size(ye,1),1))'*(ye-repmat(mean(ye),size(ye,1),1))));
                
                if e <= nb_factors
                    tmpSSe(frame) = SSt - SSs(frame) - SS{e}(frame);
                else
                    if length(interaction{i}) <= 2
                        error_effects = cell2mat(SSe(interaction{i})); error_effects = error_effects(frame,:);
                        tmpSSe(frame) = SSt - SSs(frame) - sum(main_effects(frame,:),2) - SS{e}(frame) - sum(error_effects,2); % error of interactions computed this way contains main effects/errors to be removed
                        
                    else
                        effects = cell2mat(SS); effects=effects(frame,:);
                        errors = cell2mat(SSe); errors=errors(frame,:);
                        tmpSSe(frame) = SSt - SSs(frame) - sum(effects) - sum(errors);
                    end
                end
            end
            SSe{e} = tmpSSe'; clear tmpSSe;
            
            % Stat
            disp('compute F and p values ...');
            for frame = 1:size(Y,2)
                F_values{e}(frame) = (SS{e}(frame)./(EpsHF{frame}(e)*df)) ./ (SSe{e}(frame)./(EpsHF{frame}(e).*dfe));
                p_values{e}(frame) = 1 - fcdf(F_values{e}(frame), (EpsHF{frame}(e)*df), (EpsHF{frame}(e).*dfe));
            end
        end
        
        
        % output
        % ------
        results.Design = [X Subjects];
        results.name = Name;
        results.Betas = Rep_betas;
        results.F = F_values;
        results.p = p_values;
        
        
    case{3}
        %% ------------------------------------------------------------------------------------------
        %                   Gp * Repeated measure ANOVA
        %                   Y = XB + E with X = [Factor  Gp*Factor / Gp Subject]
        % ------------------------------------------------------------------------------------------
        
        df  = nb_conditions -1;
        dfe = size(Y,1)  - nb_subjects - df;
        
        % create the design matrix for the different levels
        % ------------------------------------------------
        x = kron(eye(nb_conditions),ones(nb_subjects,1));
        S = repmat(eye(nb_subjects),nb_conditions,1);
        
        % Gp allowing different number of subjects per group (store this info in nb_items)
        Gp = zeros(size(gp,1),max(gp));
        index1 = 1; index2 = 1;
        for i=1:max(gp)
            tmp = gp(gp==i);
            nb_items(index1) = length(tmp);
            Gp(index2:index2+nb_items(index1)-1,i) = 1;
            index2 = nb_items(index1) + 1;
            index1 = index1 +1;
        end
        Gp = repmat(Gp,nb_conditions,1);
        
        % Interaction
        XGp = [];
        for i=1:size(Gp,2)
            XGp = [XGp repmat(Gp(:,i),1,nb_conditions).*x];
        end
        
        X = [Gp x XGp S];
        if flag == 1
            figure('Name','Design matrix'); set(gcf,'Color','w'); imagesc(X);
            colormap('gray');  title('ANOVA model','FontSize',16);xlabel('regressors');
            ylabel('scaled model values'); drawnow;
            go = questdlg('start the analysis?');
            if strcmp(go,'No')
                return
            end
        end
        
        
        % the results are
        % --------------
        % error / subjects
        Betas = pinv(S)*Y;
        yS    = S*Betas;
        SS_S  = diag((yS-repmat(mean(yS),size(yS,1),1))'*(yS-repmat(mean(yS),size(yS,1),1)));
        
        % Gp
        Betas = pinv(Gp)*Y;
        Rep_betas{1} = Betas;
        yGp   = Gp*Betas;
        SS_Gp  = diag((yGp-repmat(mean(yGp),size(yGp,1),1))'*(yGp-repmat(mean(yGp),size(yGp,1),1)));
        SS_S   = SS_S - SS_Gp;
        
        df_gp  = max(gp) -1;
        dfe_gp = nb_subjects - 1 - df_gp; % -1 for subjects per se -gp;
        F_values(:,1)   = (SS_Gp./df_gp) ./ (SS_S./dfe_gp);
        p_values(:,1)   = 1 - fcdf(F_values(:,1), df_gp, dfe_gp);
        
        % effect / factor
        Betas = pinv(x)*Y;
        Rep_betas{2} = Betas;
        yF    = x*Betas;
        SS_F  = diag((yF-repmat(mean(yF),size(yF,1),1))'*(yF-repmat(mean(yF),size(yF,1),1)));
        
        % interaction
        Betas = pinv(XGp)*Y;
        Rep_betas{3} = Betas;
        yXGp   = XGp*Betas;
        SS_XGp = diag((yXGp-repmat(mean(yXGp),size(yXGp,1),1))'*(yXGp-repmat(mean(yXGp),size(yXGp,1),1)));
        SS_XGp = SS_XGp - SS_Gp - SS_F;
        df_interaction = df*df_gp;
        dfe = dfe - df_interaction;
        
        % Residual error
        SS_T = diag((Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1)));
        SS_E = SS_T -SS_F - SS_Gp - SS_S - SS_XGp;
        
        % stats
        % factor
        F_values(:,2) = (SS_F./(EpsHF{1}*df)) ./ (SS_E./(EpsHF{1}*dfe));
        p_values(:,2) = 1 - fcdf(F_values(:,2),EpsHF{1}*df,EpsHF{1}*dfe);
        % interaction
        F_values(:,3) = (SS_XGp./(EpsHF{1}*df_interaction)) ./ (SS_E./(EpsHF{1}*dfe));  % there is only one correction factor
        p_values(:,3) = 1 - fcdf(F_values(:,3),EpsHF{1}*df,EpsHF{1}*dfe);
        
        
        % output
        % ------
        results.Design = X;
        results.name = ['Gp ','Factor 1 ','Interaction'];
        results.Betas = Rep_betas;
        results.F = F_values;
        results.p = p_values;
        
        
        
    case{4}
        
        %% ------------------------------------------------------------------------------------------
        %                   Gp * Repeated measure ANOVA with several factors
        %                   Y = XB + E with X = [Factors  Gp*Factors / Gp Subject]
        % ------------------------------------------------------------------------------------------
        
        % create the design matrix for the different factors
        % --------------------------------------------------
        
        % Start as case 2 for all repeated factors
        % ----------------------------------------
        
        Subjects = repmat(eye(nb_subjects),nb_conditions,1);
        
        % make all the interactions between factors
        % compute all possible interactions
        for nb_comb = 2:nb_factors
            combinations{nb_comb-1} = nchoosek([1:nb_factors],nb_comb);
        end
        
        % reorganize as a list
        index = 1; v = f+1;
        for i=1:size(combinations,2)
            selection = single(combinations{i});
            for j=1:size(selection,1)
                interaction{index}= selection(j,:);
                Name{v} = sprintf('Interaction %s',num2str(interaction{index}));
                v = v+1;index= index +1;
            end
        end
        
        % multiply factors
        for i=1:size(interaction,2)
            f = f+1;
            
            % get matrices to multiply into M
            for j=1:size(interaction{i},2)
                M{j} = F{interaction{i}(j)};
                nb_columns(j) = size(M{j},2);
                tmpN(j) = N(interaction{i}(j));
                tmp_df(j) = df_factors(interaction{i}(j));
            end
            df_factors{f} = prod(cell2mat(tmp_df));
            
            % multiply column-wise
            if size(M,2) == 2 % if only 2 matrices
                
                index = 1;
                [value,matrix1]=min(tmpN);
                [value,matrix2]=max(tmpN);
                
                for k = 1:size(M{matrix2},2) % rotates on the matrix with largest N
                    tmp{index} = M{matrix1} .* repmat(M{matrix2}(:,k),1,size(M{matrix1},2));
                    index = index + 1;
                end
                F{f} = cell2mat(tmp); clear tmp;
                
            else
                
                if i == size(interaction,2) %% all factors
                    F{f} = kron(eye(nb_conditions),ones(nb_subjects,1));
                    
                else
                    index = 1;
                    [value,matrix1]=min(tmpN);
                    [value,matrix2]=max(tmpN);
                    
                    for k = 1:size(M{matrix2},2) % rotates on the matrix with largest N
                        tmp{index} = M{matrix1} .* repmat(M{matrix2}(:,k),1,size(M{matrix1},2));
                        index = index + 1;
                    end
                    tmp2 = cell2mat(tmp); clear tmp;
                    
                    [r,c,v] = find((((interaction{i}-matrix1 == 0) + (interaction{i}-matrix2 == 0)) == 0) .* interaction{i});
                    
                    for l=v
                        index = 1;
                        for k=1:size(M{l},2)
                            tmp{index} = tmp2 .* repmat(M{l}(:,k),1,size(tmp2,2));
                            index = index + 1;
                        end
                        tmp3 = cell2mat(tmp); clear tmp;
                        clear tmp2; tmp2 = tmp3; clear tmp3
                    end
                    F{f} = tmp2; clear tmp2;
                end
            end
            clear M nb_columns tmpN
        end
        
        % Now get the gp effect and compute all interactions
        % ---------------------------------------------------
        
        % Gp allowing different nuber of subjects per group (store this info in nb_items)
        Gp = zeros(size(gp,1),max(gp));
        index1 = 1; index2 = 1;
        for i=1:max(gp)
            tmp = gp(gp==i);
            nb_items(index1) = length(tmp);
            Gp(index2:index2+nb_items(index1)-1,i) = 1;
            index2 = nb_items(index1) + 1;
            index1 = index1 +1;
        end
        Gp = repmat(Gp,nb_conditions,1);
        
        for i=1:nb_interactions
            tmp = [];
            for j=1:size(Gp,2)
                tmp = [tmp repmat(Gp(:,j),1,size(F{i},2)).*F{i}];
            end
            XGp{i} = tmp;
        end
        
        % Assemble the parts
        % -------------------
        X = [];
        for c = 1:size(F,2)
            X = [X F{c}];
        end
        
        X = [Gp X];
        for c=1:size(XGp,2)
            X = [X XGp{c}];
        end
        
        if flag == 1
            figure('Name','Design matrix'); set(gcf,'Color','w'); imagesc([X Subjects]);
            colormap('gray');  title('ANOVA model','FontSize',16);xlabel('regressors');
            ylabel('scaled model values'); drawnow;
            go = questdlg('start the analysis?');
            if strcmp(go,'No')
                return
            end
        end
        
        
        % analyse data
        % -----------------------------
        i = 0;
        
        % SS subject
        Betas = pinv(Subjects)*Y;
        yS    = Subjects*Betas;
        SSs  = diag((yS-repmat(mean(yS),size(yS,1),1))'*(yS-repmat(mean(yS),size(yS,1),1)));
        
        % Gp
        Betas = pinv(Gp)*Y;
        Rep_betas{1} = Betas;
        yGp   = Gp*Betas;
        SS_Gp  = diag((yGp-repmat(mean(yGp),size(yGp,1),1))'*(yGp-repmat(mean(yGp),size(yGp,1),1)));
        SSs   = SSs - SS_Gp;
        
        df_gp  = max(gp) -1;
        dfe_gp = nb_subjects - 1 - df_gp; % -1 for subjects per se - gp
        F_values{1}   = (SS_Gp./df_gp) ./ (SSs./dfe_gp);
        p_values{1}   = 1 - fcdf(F_values{1}, df_gp, dfe_gp);
        
        
        % SS factors and error
        index = size(F,2) +1; % use to store interaction betas; starts at+1 because gp is in column 1
        for e=1:nb_effects
            
            % SS effect - since this is the same model everywhere, we do
            % not have to run it per frame
            fprintf('compute effect, interaction and error %g out of %g',e, nb_effects); disp(' ');
            
            % Factor e or interaction
            Betas = pinv(F{e})*Y;
            Rep_betas{e} = Betas;
            yhat  = F{e}*Betas;
            SS{e} = diag((yhat-repmat(mean(yhat),size(yhat,1),1))'*(yhat-repmat(mean(yhat),size(yhat,1),1)));
            
            % interaction with Gp
            Betas = pinv(XGp{e})*Y;
            Rep_betas{index} = Betas;
            yhat  = XGp{e}*Betas;
            SS{index} = diag((yhat-repmat(mean(yhat),size(yhat,1),1))'*(yhat-repmat(mean(yhat),size(yhat,1),1)));
            SS{index} = SS{index} - SS_Gp - SS{e};
            
            df = df_factors{e};
            dfe = (nb_subjects-1)*df-(df*df_gp); % remove the df interaction from df error
            
            if e > nb_factors
                i = i+1;
                if length(interaction{i}) <= 2
                    main_effects = cell2mat(SS(interaction{i}));
                    SS{e} = SS{e} - sum(main_effects,2); % removes sum of SS main effects
                    interactions = cell2mat(SS(interaction{i}+size(F,2)));
                    SS{index} = SS{index}  - sum(interactions,2); % - sum(main_effects,2)
                else
                    main_effects_interactions = cell2mat(SS(1:e-1));
                    SS{e} = SS{e} - sum(main_effects_interactions,2); % removes sum of SS main effects and lower interactions
                    interactions = cell2mat(SS(nb_effects+1:nb_effects+e-1));
                    SS{index} = SS{index} - sum(interactions,2);
                end
            end
            index = index + 1;
            
            
            % SS error - here it depends on Y, ie we do have to run it per frame
            for frame = 1:size(Y,2)
                ye = repmat(Y(:,frame),1,size(F{e},2)).*F{e}; % data into the F{e} format
                avg_number = [1:nb_subjects];
                tmp = ye(avg_number,:);
                for avg = 2:nb_conditions
                    avg_number = avg_number + nb_subjects;
                    tmp = tmp+ye(avg_number,:);
                end
                tmp = tmp./(nb_conditions/size(F{e},2)); % produces the data averaged per level
                ye = tmp(:);
                SSt = (nb_conditions/size(F{e},2)) .* (diag((ye-repmat(mean(ye),size(ye,1),1))'*(ye-repmat(mean(ye),size(ye,1),1))));
                
                if e <= nb_factors
                    tmpSSe(frame) = SSt - SSs(frame) - SS_Gp(frame) - SS{e}(frame) - SS{index-1}(frame);
                else
                    if length(interaction{i}) <= 2
                        error_effects = cell2mat(SS(interaction{i}))+ cell2mat(SS(interaction{i}+size(F,2))) + cell2mat(SSe(interaction{i})); error_effects = error_effects(frame,:);
                        tmpSSe(frame) = SSt - SSs(frame) - SS_Gp(frame) - sum(error_effects,2) - SS{e}(frame) - SS{index-1}(frame); % error of interactions computed this way contains main effects/errors to be removed
                        
                    else
                        error_effects = cell2mat(SS(1:e-1))+ cell2mat(SS(nb_interactions+1:nb_interactions+e-1))+ cell2mat(SSe(1:e-1)); error_effects = error_effects(frame,:); % this is a bit crazy as it will work only up to 3 factors ; more need to be able to call the right cells in SS{e}
                        interaction_effects = cell2mat(SS(size(F,2)+1:size(F,2)+e)); interaction_effects = interaction_effects(frame,:);
                        tmpSSe(frame) = SSt - SSs(frame) - SS_Gp(frame) - sum(error_effects,2) - SS{e}(frame) - SS{index-1}(frame);
                    end
                end
            end
            SSe{e} = tmpSSe'; clear tmpSSe;
            
            % Stat
            disp('compute F and p values ...');
            for frame = 1:size(Y,2)
                F_values{e+1}(frame) = (SS{e}(frame)./(EpsHF{frame}(e)*df)) ./ (SSe{e}(frame)./(EpsHF{frame}(e).*dfe));
                p_values{e+1}(frame) = 1 - fcdf(F_values{e+1}(frame), (EpsHF{frame}(e)*df), (EpsHF{frame}(e).*dfe));
                F_values{index}(frame) = (SS{index-1}(frame)./(EpsHF{frame}(e)*df*df_gp)) ./ (SSe{e}(frame)./(EpsHF{frame}(e).*dfe));
                p_values{index}(frame) = 1 - fcdf(F_values{index}(frame), (EpsHF{frame}(e)*df*df_gp), (EpsHF{frame}(e).*dfe));
            end
        end
        
        
        % output
        % ------
        results.Design = [X Subjects];
        clear tmp
        for i=1:size(Name,2)
            tmp{i} = [Name{i} ' * Gp'];
        end
        results.name = [Name tmp];
        results.Betas = Rep_betas;
        results.F = F_values;
        results.p = p_values;
        
end


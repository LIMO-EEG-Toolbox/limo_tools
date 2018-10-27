function limotest_glm(option)

% this is the unit test of the limo_glm function
% no input is required, instead a series of known values are passed on for
% the various options available and compared to the known output - know
% input/output comes from alternative software or solutions (eg robustfit)

nmc = 10000;

if strcmpi(option,'OLS')
    %% Simple Regression
    
    Y = randn(100,1);
    X = [randn(100,1) ones(100,1)];
    
    model = limo_glm(Y,X,0,0,1,'OLS','Time',0,1);
    [B,~,~,~,STATS] = regress(Y,X);
    if any((single(model.betas) == single(B)))
        disp('Simple regression OLS: same parameter estimates as Matlab regress')
    else
        disp('Simple regression OLS: different parameter estimates as Matlab regress')
    end
    
    if ~any([single(model.R2_univariate) == single(STATS(1)), ...
            single(model.F) == single(STATS(2)), ...
            single(model.p) == single(STATS(3))])
        disp('Simple regression OLS: different model stats values as Matlab regress')
    else
        disp('Simple regression OLS: same model stats values as Matlab regress')
    end
    
    
    %% Multiple Regression
    
    Y = randn(100,1);
    X = [randn(100,3) ones(100,1)];
    
    model = limo_glm(Y,X,0,0,3,'OLS','Time',0,5);
    [B,~,~,~,STATS] = regress(Y,X);
    if any((single(model.betas) == single(B)))
        disp('Multiple regression OLS: same parameter estimates as Matlab regress')
    else
        disp('Multiple regression OLS: different parameter estimates as Matlab regress')
    end
    
    if ~any([single(model.R2_univariate) == single(STATS(1)), ...
            single(model.F) == single(STATS(2)), ...
            single(model.p) == single(STATS(3))])
        disp('Multiple regression OLS: different stats values as Matlab regress')
    else
        disp('Multiple regression OLS: same model stats values as Matlab regress')
    end
    
    %% One-Way ANOVA
    
    Y = randn(90,1);
    X = [kron(eye(3),ones(30,1)) ones(90,1)];
    
    model = limo_glm(Y,X,1,0,0,'OLS','Time',0,1);
    group = X(:,1)+X(:,2).*2+X(:,3).*3;
    [~,T]=anovan(Y,group,'display','off');
    if ~any([single(T{2,6}) == single(model.F), ...
            single(T{2,7}) == single(model.p), ...
            single(T{2,3}) == single(model.df(1)), ...
            single(T{3,3}) == single(model.df(2)), ...
            single(T{2,2}/T{4,2}) == single(model.R2_univariate)])
        disp('One Way ANOVA OLS: different model stats values as Matlab anovan')
    else
        disp('One Way ANOVA OLS: same model stats values as Matlab anovan')
    end
    
    
    %% One-Way ANCOVA
    
    Y = randn(90,1);
    X = [kron(eye(3),ones(30,1)) randn(90,2) ones(90,1)];
    
    model = limo_glm(Y,X,3,0,2,'OLS','Time',0,1);
    group = X(:,1)+X(:,2).*2+X(:,3).*3;
    cov1  = X(:,4); cov2  = X(:,5);
    [~,T]=anovan(Y,{group,cov1,cov2},'continuous',[2 3],'display','off');
    
    if ~any([single(T{2,6}) == single(model.conditions.F), ...
            single(T{2,7}) == single(model.conditions.p), ...
            single(T{3,6}) == single(model.continuous.F(1)), ...
            single(T{3,7}) == single(model.continuous.p(1)), ...
            single(T{4,6}) == single(model.continuous.F(2)), ...
            single(T{4,7}) == single(model.continuous.p(2)), ...
            single(T{5,2}/T{6,2}) == single(model.R2_univariate)])
        disp('One Way ANCOVA OLS: different model stats values as Matlab anovan')
    else
        disp('One Way ANCOVA OLS: same model stats values as Matlab anovan')
    end
    
    
    %% N-ways ANOVA without interactions
    Y = randn(240,1);
    X = [kron(eye(3),ones(80,1)) repmat(kron(eye(2),ones(40,1)),3,1) ...
        repmat(kron(eye(2),ones(20,1)),6,1) ones(240,1)];
    
    model = limo_glm(Y,X,[3 2 2],0,0,'OLS','Time',0,1);
    group1 = X(:,1)+X(:,2).*2+X(:,3).*3;
    group2 = X(:,4)+X(:,5).*2;
    group3 = X(:,6)+X(:,7).*2;
    [~,T]=anovan(Y,{group1,group2,group3},'display','off');
    
    if ~any([single(T{2,6}) == single(model.conditions.F(1)), ...
            single(T{2,7}) == single(model.conditions.p(1)), ...
            single(T{3,6}) == single(model.conditions.F(2)), ...
            single(T{3,7}) == single(model.conditions.p(2)), ...
            single(T{4,6}) == single(model.conditions.F(3)), ...
            single(T{4,7}) == single(model.conditions.p(3)), ...
            single((T{2,2}+T{3,2}+T{4,2})/T{6,2}) == single(model.R2_univariate)])
        disp('Three Way ANOVA OLS: different model stats values as Matlab anovan')
    else
        disp('Three Way ANOVA OLS: same model stats values as Matlab anovan')
    end
    
    
    %% N-ways ANCOVA without interactions
    Y = randn(100,1);
    X = [kron(eye(2),ones(50,1)) repmat(kron(eye(2),ones(25,1)),2,1) ...
        randn(100,2) ones(100,1)];
    
    model = limo_glm(Y,X,[2 2],0,2,'OLS','Time',0,1);
    group1 = X(:,1)+X(:,2).*2;
    group2 = X(:,3)+X(:,4).*2;
    cov1   = X(:,5); cov2 = X(:,6);
    [~,T]=anovan(Y,{group1,group2,cov1,cov2},'continuous',[3 4],'display','off');
    
    if ~any([single(T{2,6}) == single(model.conditions.F(1)), ...
            single(T{2,7}) == single(model.conditions.p(1)), ...
            single(T{3,6}) == single(model.conditions.F(2)), ...
            single(T{3,7}) == single(model.conditions.p(2)), ...
            single(T{4,6}) == single(model.continuous.F(1)), ...
            single(T{4,7}) == single(model.continuous.p(1)), ...
            single(T{5,6}) == single(model.continuous.F(2)), ...
            single(T{5,7}) == single(model.continuous.p(2)), ...
            single((T{2,2}+T{3,2}+T{4,2}+T{5,2})/T{7,2}) == single(model.R2_univariate)])
        disp('Two Way ANCOVA OLS: different model stats values as Matlab anovan')
    else
        disp('Two Way ANCOVA OLS: same model stats values as Matlab anovan')
    end
    
    
    %% N-ways ANOVA with interactions (full factorial)
    
    Y = randn(360,1);
    tmpX = [kron(eye(3),ones(120,1)) repmat(kron(eye(2),ones(60,1)),3,1) ...
        repmat(kron(eye(2),ones(30,1)),6,1)];
    [X,interactions] = limo_make_interactions(tmpX, [3 2 2]);
    X = [X ones(360,1)];
    
    model = limo_glm(Y,X,[3 2 2],interactions,0,'OLS','Time',0,1);
    group1 = X(:,1)+X(:,2).*2+X(:,3).*3;
    group2 = X(:,4)+X(:,5).*2;
    group3 = X(:,6)+X(:,7).*2;
    [~,T]=anovan(Y,{group1,group2,group3},'model','full','display','off');
    
    if ~any([single(T{2,6}) == single(model.conditions.F(1)), ...
            single(T{2,7}) == single(model.conditions.p(1)), ...
            single(T{3,6}) == single(model.conditions.F(2)), ...
            single(T{3,7}) == single(model.conditions.p(2)), ...
            single(T{4,6}) == single(model.conditions.F(3)), ...
            single(T{4,7}) == single(model.conditions.p(3)), ...
            single((T{2,2}+T{3,2}+T{4,2})/T{6,2}) == single(model.R2_univariate)])
        disp('Full factorial three way ANOVA OLS: different model stats values as Matlab anovan')
    else
        disp('Full factorial three way ANOVA OLS: same model stats values as Matlab anovan')
    end
    
else strcmpi(option,'IRLS')
    
    %% simply check the data strcuture is the same
    % Simple Regression
    Y = randn(100,5);
    X = [randn(100,1) ones(100,1)];
    modelo = limo_glm(Y,X,0,0,1,'OLS','Time',0,1);
    modeli = limo_glm(Y,X,0,0,1,'IRLS','Time',0,1);
    fields = fieldnames(modelo);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'W')
            if any(size(getfield(modelo,fields{f})) ~= size(getfield(modeli,fields{f})))
                error('different field size model.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    % Multiple Regression
    Y = randn(100,5);
    X = [randn(100,3) ones(100,1)];
    modelo = limo_glm(Y,X,0,0,3,'OLS','Time',0,5);
    modeli = limo_glm(Y,X,0,0,3,'IRLS','Time',0,1);

    fields = fieldnames(modelo);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'W')
            if any(size(getfield(modelo,fields{f})) ~= size(getfield(modeli,fields{f})))
                error('different field size model.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.continuous);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.continuous,fields{f})) ~= size(getfield(modeli.continuous,fields{f})))
                error('different field size model.continuous.%s between models OLS vs IRLS',fields{f})
            end
        end
    end

    % One-Way ANOVA
    Y = randn(90,5);
    X = [kron(eye(3),ones(30,1)) ones(90,1)];
    modelo = limo_glm(Y,X,1,0,0,'OLS','Time',0,1);
    modeli = limo_glm(Y,X,1,0,0,'IRLS','Time',0,1);
   
    fields = fieldnames(modelo);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'W')
            if any(size(getfield(modelo,fields{f})) ~= size(getfield(modeli,fields{f})))
                error('different field size model.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.conditions);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.conditions,fields{f})) ~= size(getfield(modeli.conditions,fields{f})))
                error('different field size model.conditions.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    % One-Way ANCOVA
    Y = randn(90,5);
    X = [kron(eye(3),ones(30,1)) randn(90,2) ones(90,1)];
    modelo = limo_glm(Y,X,3,0,2,'OLS','Time',0,1);
    modeli = limo_glm(Y,X,3,0,2,'IRLS','Time',0,1);
    
    fields = fieldnames(modelo);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'W')
            if any(size(getfield(modelo,fields{f})) ~= size(getfield(modeli,fields{f})))
                error('different field size model.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.conditions);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.conditions,fields{f})) ~= size(getfield(modeli.conditions,fields{f})))
                error('different field size model.conditions.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.continuous);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.continuous,fields{f})) ~= size(getfield(modeli.continuous,fields{f})))
                error('different field size model.continuous.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    % N-ways ANOVA without interactions
    Y = randn(240,5);
    X = [kron(eye(3),ones(80,1)) repmat(kron(eye(2),ones(40,1)),3,1) ...
        repmat(kron(eye(2),ones(20,1)),6,1) ones(240,1)];
    modelo = limo_glm(Y,X,[3 2 2],0,0,'OLS','Time',0,1);
    modeli = limo_glm(Y,X,[3 2 2],0,0,'IRLS','Time',0,1);

    fields = fieldnames(modelo);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'W')
            if any(size(getfield(modelo,fields{f})) ~= size(getfield(modeli,fields{f})))
                error('different field size model.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.conditions);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.conditions,fields{f})) ~= size(getfield(modeli.conditions,fields{f})))
                error('different field size model.conditions.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    % N-ways ANCOVA without interactions
    Y = randn(100,5);
    X = [kron(eye(2),ones(50,1)) repmat(kron(eye(2),ones(25,1)),2,1) ...
        randn(100,2) ones(100,1)];
    modelo = limo_glm(Y,X,[2 2],0,2,'OLS','Time',0,1);
    modeli = limo_glm(Y,X,[2 2],0,2,'IRLS','Time',0,1);
    
    fields = fieldnames(modelo);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'W')
            if any(size(getfield(modelo,fields{f})) ~= size(getfield(modeli,fields{f})))
                error('different field size model.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.conditions);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.conditions,fields{f})) ~= size(getfield(modeli.conditions,fields{f})))
                error('different field size model.conditions.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.continuous);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.continuous,fields{f})) ~= size(getfield(modeli.continuous,fields{f})))
                error('different field size model.continuous.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    % N-ways ANOVA with interactions (full factorial)
    Y = randn(360,5);
    tmpX = [kron(eye(3),ones(120,1)) repmat(kron(eye(2),ones(60,1)),3,1) ...
        repmat(kron(eye(2),ones(30,1)),6,1)];
    [X,interactions] = limo_make_interactions(tmpX, [3 2 2]);
    X = [X ones(360,1)];
    modelo = limo_glm(Y,X,[3 2 2],interactions,0,'OLS','Time',0,1);
    modeli = limo_glm(Y,X,[3 2 2],interactions,0,'IRLS','Time',0,1);
 
    fields = fieldnames(modelo);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'W')
            if any(size(getfield(modelo,fields{f})) ~= size(getfield(modeli,fields{f})))
                error('different field size model.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    fields = fieldnames(modelo.conditions);
    for f=1:size(fields,1)
        if ~strcmp(fields{f},'df')
            if any(size(getfield(modelo.conditions,fields{f})) ~= size(getfield(modeli.conditions,fields{f})))
                error('different field size model.conditions.%s between models OLS vs IRLS',fields{f})
            end
        end
    end
    
    disp('Data structure OLS and IRLS are equivalent (df and weights differ - not tested)')

    %% run simulations
    % since IRLS implementation depends on the tuning function, scaling, etc we
    % cannot expect the same results - but we can make sure we return the right
    % type 1 error rate (i.e. if not the same there is a mistake)
    q = questdlg('run simlations? (takes few days)','Monte Carlo validation','yes','no','no');
    
    if strcmp(q,'yes')
        
        %% Simple Regression
        ols  = NaN(1,100);
        irls = NaN(1,100);
        for avg = 1:100
            Y2 = randn(100,nmc);
            X2 = [randn(100,1) ones(100,1)];
            pvaluesO = NaN(1,nmc);
            pvaluesI = NaN(1,nmc);
            
            parfor MC=1:nmc
                model = limo_glm(Y2(:,MC),X2,0,0,1,'OLS','Time',0,1);
                pvaluesO(MC) = model.p;
                model = limo_glm(Y2(:,MC),X2,0,0,1,'IRLS','Time',0,1);
                pvaluesI(MC) = model.p;
            end
            pvaluesO(isnan(pvaluesO)) = []; ols(avg) = sum(pvaluesO<0.05)/length(pvaluesO);
            pvaluesI(isnan(pvaluesI)) = []; irls(avg) = sum(pvaluesI<0.05)/length(pvaluesI);
        end
        ols = ols.*100;
        irls = irls.*100;
        
        figure;
        subplot(1,2,1);
        scatter(ones(100,1),ols,30,'k'); hold on;
        scatter(ones(100,1).*2,irls,30,'k');
        for pair = 1:100
            plot([1,2],[ols(pair) irls(pair)]);
        end
        axis([0.8 2.2 4 6]); box on; grid on;
        [h,CI] = rst_1ttest([ols' irls'],'mean',0);
        title(sprintf('Type 1 error rate \n OLS %g [%g %g] \n IRLS %g [%g %g]', ...
            mean(ols),CI(1,1),CI(2,1),mean(irls),CI(1,2),CI(2,2)),'FontSize',12)
        ylabel('type 1 error rate')
        
        subplot(1,2,2);
        ols_irls_diff = (ols-irls);
        [meandiff,cidiff]=rst_data_plot(ols_irls_diff,'estimator','mean','newfig','no');
        title(sprintf('Simple regression \n OLS - IRLS = %g \n CI=[%g %g]', ...
            meandiff,cidiff(1),cidiff(2)),'FontSize',12)
        ylabel('difference in type 1 error rate')
        
        % --------------------------------
        if all(CI(:,1)<5) || all(CI(:,1)>5)
            fprintf('Simple regression OLS is not at type 1 error 5%% %g [%g %g]\n',mean(ols),CI(1,1),CI(2,1))
        else
            fprintf('Simple regression OLS is at type 1 error 5%% %g [%g %g]\n',mean(ols),CI(1,1),CI(2,1))
        end
        
        if all(CI(:,2)<5) || all(CI(:,2)>5)
            fprintf('Simple regression IRLS is not at type 1 error 5%% %g [%g %g]\n',mean(ols),CI(1,1),CI(2,1))
        else
            fprintf('Simple regression IRLS is at type 1 error 5%% %g [%g %g]\n',mean(irls),CI(1,2),CI(2,2))
        end
        
        if all(cidiff<0) || all(cidiff>0)
            fprintf('Simple regression: there is a significant difference in type 1 error \n between OLS %g and IRLS %g diff=%g%%\n',mean(ols),mean(irls),meandiff)
        else
            fprintf('Simple regression: there is no significant difference in type 1 error \n between OLS %g and IRLS %g diff=%g%%\n',mean(ols),mean(irls),meandiff)
        end
        
        %% Multiple Regression
        
        ols  = NaN(4,100);
        irls = NaN(4,100);
        for avg = 1:100
            Y2 = randn(100,nmc);
            X2 = [randn(100,3) ones(100,1)];
            pvaluesO = NaN(nmc,4);
            pvaluesI = NaN(nmc,4);
            
            parfor MC=1:nmc
                model = limo_glm(Y2(:,MC),X2,0,0,3,'OLS','Time',0,1);
                pvaluesO(MC,:) = [model.p model.continuous.p];
                model = limo_glm(Y2(:,MC),X2,0,0,3,'IRLS','Time',0,1);
                pvaluesI(MC,:) = [model.p model.continuous.p];
            end
            
            for test = 1:4
                tmp = pvaluesO(:,test); tmp(isnan(tmp)) = [];
                ols(test,avg) = sum(tmp<0.05)/length(tmp);
                tmp = pvaluesI(:,test); tmp(isnan(tmp)) = [];
                irls(test,avg) = sum(tmp<0.05)/length(tmp);
            end
        end
        ols = ols.*100;
        irls = irls.*100;
        
        figure;
        subplot(1,2,1);
        scatter(ones(100,1),ols(1,:),30,'k'); hold on;
        scatter(ones(100,1).*2,irls(1,:),30,'k');
        for pair = 1:100
            plot([1,2],[ols(1,pair) irls(1,pair)]);
        end
        axis([0.8 2.2 4 6]); box on; grid on;
        [~,CI] = rst_1ttest([ols' irls'],'mean',0);
        title(sprintf('Type 1 error rate \n OLS %g [%g %g] \n IRLS %g [%g %g]', ...
            mean(ols(1,:)),CI(1,1),CI(2,1),mean(irls(1,:)),CI(1,5),CI(2,5)),'FontSize',12)
        ylabel('type 1 error rate')
        
        subplot(1,2,2);
        ols_irls_diff = (ols'-irls');
        [meandiff,cidiff]=rst_data_plot(ols_irls_diff(:,1),'estimator','mean','newfig','no');
        title(sprintf('Multiple regression \n OLS - IRLS = %g \n CI=[%g %g]', ...
            meandiff,cidiff(1),cidiff(2)),'FontSize',12)
        ylabel('difference in type 1 error rate')
        
        % --------------------------------
        if sum(all(CI(:,1:4)<5))~=0 || sum(all(CI(:,1:4)>5))~=0
            fprintf('Multiple regression OLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        else
            fprintf('Multiple regression OLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        end
        
        if sum(all(CI(:,5:8)<5))~=0 || sum(all(CI(:,5:8)>5))~=0
            fprintf('Multiple regression IRLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        else
            fprintf('Multiple regression IRLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        end
        
        [~,cidiff] = rst_1ttest(ols_irls_diff,'mean',0);
        meandiff = mean(ols_irls_diff);
        
        nb_diff = max([sum(all(cidiff<0)) sum(all(cidiff>0))]);
        if nb_diff~=0
            if nb_diff == 1
                disp('Multiple regression: there is 1 significant difference in type 1 error between OLS and IRLS ')
            else
                fprintf('Multiple regression: there are %g significant differences in type 1 error between OLS and IRLS \n',nb_diff)
            end
        else
            disp('Multiple regression: there is no significant difference in type 1 error between OLS and IRLS')
        end
        
        T = array2table([mean(ols')' CI(:,1:4)' mean(irls')' CI(:,5:8)' meandiff' cidiff(1,:)' cidiff(2,:)'],...
            'VariableNames',{'MeanOLS','CIOlow','CIOhigh','MeanIRLS','CIIlow','CIIhigh' 'MeanDiff','CIDlow','CIDhigh'});
        disp(T)
        save('multiple_reg_table','T')
        
        %% One-Way ANOVA
        
        ols  = NaN(1,100);
        irls = NaN(1,100);
        X = [kron(eye(3),ones(30,1)) ones(90,1)];
        for avg = 1:100
            Y2 = randn(90,nmc);
            pvaluesO = NaN(1,nmc);
            pvaluesI = NaN(1,nmc);
            
            parfor MC=1:nmc
                model = limo_glm(Y2(:,MC),X,1,0,0,'OLS','Time',0,1);
                pvaluesO(MC) = model.p;
                model = limo_glm(Y2(:,MC),X,1,0,0,'IRLS','Time',0,1);
                pvaluesI(MC) = model.p;
            end
            pvaluesO(isnan(pvaluesO)) = []; ols(avg) = sum(pvaluesO<0.05)/length(pvaluesO);
            pvaluesI(isnan(pvaluesI)) = []; irls(avg) = sum(pvaluesI<0.05)/length(pvaluesI);
        end
        ols = ols.*100;
        irls = irls.*100;
        
        figure;
        subplot(1,2,1);
        scatter(ones(100,1),ols,30,'k'); hold on;
        scatter(ones(100,1).*2,irls,30,'k');
        for pair = 1:100
            plot([1,2],[ols(pair) irls(pair)]);
        end
        axis([0.8 2.2 4 6]); box on; grid on;
        [~,CI] = rst_1ttest([ols' irls'],'mean',0);
        title(sprintf('Type 1 error rate \n OLS %g [%g %g] \n IRLS %g [%g %g]', ...
            mean(ols),CI(1,1),CI(2,1),mean(irls),CI(1,2),CI(2,2)),'FontSize',12)
        ylabel('type 1 error rate')
        
        subplot(1,2,2);
        ols_irls_diff = (ols-irls);
        [meandiff,cidiff]=rst_data_plot(ols_irls_diff,'estimator','mean','newfig','no');
        title(sprintf('One Way-ANOVA \n OLS - IRLS = %g \n CI=[%g %g]', ...
            meandiff,cidiff(1),cidiff(2)),'FontSize',12)
        ylabel('difference in type 1 error rate')
        
        % --------------------------------
        if all(CI(:,1)<5) || all(CI(:,1)>5)
            fprintf('One Way-ANOVA OLS is not at type 1 error 5%% %g [%g %g]\n',mean(ols),CI(1,1),CI(2,1))
        else
            fprintf('One Way-ANOVA OLS is at type 1 error 5%% %g [%g %g]\n',mean(ols),CI(1,1),CI(2,1))
        end
        
        if all(CI(:,2)<5) || all(CI(:,2)>5)
            fprintf('One Way-ANOVA IRLS is not at type 1 error 5%% %g [%g %g]\n',mean(ols),CI(1,1),CI(2,1))
        else
            fprintf('One Way-ANOVA IRLS is at type 1 error 5%% %g [%g %g]\n',mean(irls),CI(1,2),CI(2,2))
        end
        
        if all(cidiff<0) || all(cidiff>0)
            fprintf('One Way-ANOVA: there is a significant difference in type 1 error \n between OLS %g and IRLS %g diff=%g%%\n',mean(ols),mean(irls),meandiff)
        else
            fprintf('One Way-ANOVA: there is no significant difference in type 1 error \n between OLS %g and IRLS %g diff=%g%%\n',mean(ols),mean(irls),meandiff)
        end
        
        %% One-Way ANCOVA
        
        ols  = NaN(4,100);
        irls = NaN(4,100);
        for avg = 1:100
            Y2 = randn(90,nmc);
            X2 = [kron(eye(3),ones(30,1)) randn(90,2) ones(90,1)];
            pvaluesO = NaN(nmc,4);
            pvaluesI = NaN(nmc,4);
            
            parfor MC=1:nmc
                model = limo_glm(Y2(:,MC),X2,3,0,2,'OLS','Time',0,1);
                pvaluesO(MC,:) = [model.p model.conditions.p model.continuous.p];
                model = limo_glm(Y2(:,MC),X2,3,0,2,'IRLS','Time',0,1);
                pvaluesI(MC,:) = [model.p model.conditions.p model.continuous.p];
            end
            
            for test = 1:4
                tmp = pvaluesO(:,test); tmp(isnan(tmp)) = [];
                ols(test,avg) = sum(tmp<0.05)/length(tmp);
                tmp = pvaluesI(:,test); tmp(isnan(tmp)) = [];
                irls(test,avg) = sum(tmp<0.05)/length(tmp);
            end
        end
        ols = ols.*100;
        irls = irls.*100;
        
        figure;
        subplot(1,2,1);
        scatter(ones(100,1),ols(1,:),30,'k'); hold on;
        scatter(ones(100,1).*2,irls(1,:),30,'k');
        for pair = 1:100
            plot([1,2],[ols(1,pair) irls(1,pair)]);
        end
        axis([0.8 2.2 4 6]); box on; grid on;
        [~,CI] = rst_1ttest([ols' irls'],'mean',0);
        title(sprintf('Type 1 error rate \n OLS %g [%g %g] \n IRLS %g [%g %g]', ...
            mean(ols(1,:)),CI(1,1),CI(2,1),mean(irls(1,:)),CI(1,5),CI(2,5)),'FontSize',12)
        ylabel('type 1 error rate')
        
        subplot(1,2,2);
        ols_irls_diff = (ols'-irls');
        [meandiff,cidiff]=rst_data_plot(ols_irls_diff(:,1),'estimator','mean','newfig','no');
        title(sprintf('One-Way ANCOVA \n OLS - IRLS = %g \n CI=[%g %g]', ...
            meandiff,cidiff(1),cidiff(2)),'FontSize',12)
        ylabel('difference in type 1 error rate')
        
        % --------------------------------
        if sum(all(CI(:,1:4)<5))~=0 || sum(all(CI(:,1:4)>5))~=0
            fprintf('One-Way ANCOVA OLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        else
            fprintf('One-Way ANCOVA OLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        end
        
        if sum(all(CI(:,5:8)<5))~=0 || sum(all(CI(:,5:8)>5))~=0
            fprintf('One-Way ANCOVA IRLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        else
            fprintf('One-Way ANCOVA IRLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        end
        
        [~,cidiff] = rst_1ttest(ols_irls_diff,'mean',0);
        meandiff = mean(ols_irls_diff);
        
        nb_diff = max([sum(all(cidiff<0)) sum(all(cidiff>0))]);
        if nb_diff~=0
            if nb_diff == 1
                disp('One-Way ANCOVA: there is 1 significant difference in type 1 error between OLS and IRLS ')
            else
                fprintf('One-Way ANCOVA: there are %g significant differences in type 1 error between OLS and IRLS \n',nb_diff)
            end
        else
            disp('One-Way ANCOVA: there is no significant difference in type 1 error between OLS and IRLS')
        end
        
        T = array2table([mean(ols')' CI(:,1:4)' mean(irls')' CI(:,5:8)' meandiff' cidiff(1,:)' cidiff(2,:)'],...
            'VariableNames',{'MeanOLS','CIOlow','CIOhigh','MeanIRLS','CIIlow','CIIhigh' 'MeanDiff','CIDlow','CIDhigh'});
        disp(T)
        save('one-way-ANCOVA_table','T')
        
        
        %% N-ways ANOVA without interactions
        ols  = NaN(4,100);
        irls = NaN(4,100);
        X = [kron(eye(3),ones(80,1)) repmat(kron(eye(2),ones(40,1)),3,1) ...
            repmat(kron(eye(2),ones(20,1)),6,1) ones(240,1)];
        for avg = 1:100
            Y2 = randn(240,nmc);
            pvaluesO = NaN(nmc,4);
            pvaluesI = NaN(nmc,4);
            
            parfor MC=1:nmc
                model = limo_glm(Y2(:,MC),X,[3 2 2],0,0,'OLS','Time',0,1);
                pvaluesO(MC,:) = [model.p model.conditions.p'];
                model = limo_glm(Y2(:,MC),X,[3 2 2],0,0,'IRLS','Time',0,1);
                pvaluesI(MC,:) = [model.p model.conditions.p'];
            end
            
            for test = 1:4
                tmp = pvaluesO(:,test); tmp(isnan(tmp)) = [];
                ols(test,avg) = sum(tmp<0.05)/length(tmp);
                tmp = pvaluesI(:,test); tmp(isnan(tmp)) = [];
                irls(test,avg) = sum(tmp<0.05)/length(tmp);
            end
        end
        ols = ols.*100;
        irls = irls.*100;
        
        figure;
        subplot(1,2,1);
        scatter(ones(100,1),ols(1,:),30,'k'); hold on;
        scatter(ones(100,1).*2,irls(1,:),30,'k');
        for pair = 1:100
            plot([1,2],[ols(1,pair) irls(1,pair)]);
        end
        axis([0.8 2.2 4 6]); box on; grid on;
        [~,CI] = rst_1ttest([ols' irls'],'mean',0);
        title(sprintf('Type 1 error rate \n OLS %g [%g %g] \n IRLS %g [%g %g]', ...
            mean(ols(1,:)),CI(1,1),CI(2,1),mean(irls(1,:)),CI(1,5),CI(2,5)),'FontSize',12)
        ylabel('type 1 error rate')
        
        subplot(1,2,2);
        ols_irls_diff = (ols'-irls');
        [meandiff,cidiff]=rst_data_plot(ols_irls_diff(:,1),'estimator','mean','newfig','no');
        title(sprintf('Three-Way ANOVA \n OLS - IRLS = %g \n CI=[%g %g]', ...
            meandiff,cidiff(1),cidiff(2)),'FontSize',12)
        ylabel('difference in type 1 error rate')
        
        % --------------------------------
        if sum(all(CI(:,1:4)<5))~=0 || sum(all(CI(:,1:4)>5))~=0
            fprintf('Three-Way ANOVA OLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        else
            fprintf('Three-Way ANOVA OLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        end
        
        if sum(all(CI(:,5:8)<5))~=0 || sum(all(CI(:,5:8)>5))~=0
            fprintf('Three-Way ANCOVA IRLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        else
            fprintf('Three-Way ANCOVA IRLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        end
        
        [~,cidiff] = rst_1ttest(ols_irls_diff,'mean',0);
        meandiff = mean(ols_irls_diff);
        
        nb_diff = max([sum(all(cidiff<0)) sum(all(cidiff>0))]);
        if nb_diff~=0
            if nb_diff == 1
                disp('Three-Way ANCOVA: there is 1 significant difference in type 1 error between OLS and IRLS ')
            else
                fprintf('Three-Way ANCOVA: there are %g significant differences in type 1 error between OLS and IRLS \n',nb_diff)
            end
        else
            disp('Three-Way ANCOVA: there is no significant difference in type 1 error between OLS and IRLS')
        end
        
        T = array2table([mean(ols')' CI(:,1:4)' mean(irls')' CI(:,5:8)' meandiff' cidiff(1,:)' cidiff(2,:)'],...
            'VariableNames',{'MeanOLS','CIOlow','CIOhigh','MeanIRLS','CIIlow','CIIhigh' 'MeanDiff','CIDlow','CIDhigh'});
        disp(T)
        save('Three-Way-ANOVA_table','T')
        
        
        %% N-ways ANCOVA without interactions
        ols  = NaN(5,100);
        irls = NaN(5,100);
        for avg = 1:100
            fprintf('running simulation %g/100 \n',avg)
            Y2 = randn(100,nmc);
            X2 = [kron(eye(2),ones(50,1)) repmat(kron(eye(2),ones(25,1)),2,1) ...
                randn(100,2) ones(100,1)];
            pvaluesO = NaN(nmc,5);
            pvaluesI = NaN(nmc,5);
            
            parfor MC=1:nmc
                model = limo_glm(Y2(:,MC),X2,[2 2],0,2,'OLS','Time',0,1);
                pvaluesO(MC,:) = [model.p model.conditions.p' model.continuous.p];
                model = limo_glm(Y2(:,MC),X2,[2 2],0,2,'IRLS','Time',0,1);
                pvaluesI(MC,:) = [model.p model.conditions.p' model.continuous.p];
            end
            
            for test = 1:5
                tmp = pvaluesO(:,test); tmp(isnan(tmp)) = [];
                ols(test,avg) = sum(tmp<0.05)/length(tmp);
                tmp = pvaluesI(:,test); tmp(isnan(tmp)) = [];
                irls(test,avg) = sum(tmp<0.05)/length(tmp);
            end
        end
        ols = ols.*100;
        irls = irls.*100;
        
        figure;
        subplot(1,2,1);
        scatter(ones(100,1),ols(1,:),30,'k'); hold on;
        scatter(ones(100,1).*2,irls(1,:),30,'k');
        for pair = 1:100
            plot([1,2],[ols(1,pair) irls(1,pair)]);
        end
        axis([0.8 2.2 4 6]); box on; grid on;
        [~,CI] = rst_1ttest([ols' irls'],'mean',0);
        title(sprintf('Type 1 error rate \n OLS %g [%g %g] \n IRLS %g [%g %g]', ...
            mean(ols(1,:)),CI(1,1),CI(2,1),mean(irls(1,:)),CI(1,5),CI(2,5)),'FontSize',12)
        ylabel('type 1 error rate')
        
        subplot(1,2,2);
        ols_irls_diff = (ols'-irls');
        [meandiff,cidiff]=rst_data_plot(ols_irls_diff(:,1),'estimator','mean','newfig','no');
        title(sprintf('Two-Way ANCOVA \n OLS - IRLS = %g \n CI=[%g %g]', ...
            meandiff,cidiff(1),cidiff(2)),'FontSize',12)
        ylabel('difference in type 1 error rate')
        
        % --------------------------------
        if sum(all(CI(:,1:5)<5))~=0 || sum(all(CI(:,1:5)>5))~=0
            fprintf('Two-Way ANCOVA OLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        else
            fprintf('Two-Way ANCOVA OLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4))
        end
        
        if sum(all(CI(:,6:10)<5))~=0 || sum(all(CI(:,6:10)>5))~=0
            fprintf('Two-Way ANCOVA IRLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        else
            fprintf('Two-Way ANCOVA IRLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        end
        
        [~,cidiff] = rst_1ttest(ols_irls_diff,'mean',0);
        meandiff = mean(ols_irls_diff);
        
        nb_diff = max([sum(all(cidiff<0)) sum(all(cidiff>0))]);
        if nb_diff~=0
            if nb_diff == 1
                disp('Two-Way ANCOVA: there is 1 significant difference in type 1 error between OLS and IRLS ')
            else
                fprintf('Two-Way ANCOVA: there are %g significant differences in type 1 error between OLS and IRLS \n',nb_diff)
            end
        else
            disp('Two-Way ANCOVA: there is no significant difference in type 1 error between OLS and IRLS')
        end
        
        T = array2table([mean(ols')' CI(:,1:5)' mean(irls')' CI(:,6:10)' meandiff' cidiff(1,:)' cidiff(2,:)'],...
            'VariableNames',{'MeanOLS','CIOlow','CIOhigh','MeanIRLS','CIIlow','CIIhigh' 'MeanDiff','CIDlow','CIDhigh'});
        disp(T)
        save('Two-Way ANCOVA_table','T')
        
        
        %% N-ways ANOVA with interactions (full factorial)
        
        ols  = NaN(8,100);
        irls = NaN(8,100);
        tmpX = [kron(eye(3),ones(240,1)) repmat(kron(eye(2),ones(120,1)),3,1) ...
            repmat(kron(eye(2),ones(60,1)),6,1)];
        [X,interactions] = limo_make_interactions(tmpX, [3 2 2]);
        X = [X ones(720,1)];
        
        for avg = 1:100
            fprintf('running simulation %g/100 \n',avg)
            Y2 = randn(720,nmc);
            pvaluesO = NaN(nmc,8);
            pvaluesI = NaN(nmc,8);
            
            parfor MC=1:nmc
                model = limo_glm(Y2(:,MC),X,[3 2 2],interactions,0,'OLS','Time',0,1);
                pvaluesO(MC,:) = [model.p model.conditions.p' model.interactions.p'];
                model = limo_glm(Y2(:,MC),X,[3 2 2],interactions,0,'IRLS','Time',0,1);
                pvaluesI(MC,:) = [model.p model.conditions.p' model.interactions.p'];
            end
            
            for test = 1:8
                tmp = pvaluesO(:,test); tmp(isnan(tmp)) = [];
                ols(test,avg) = sum(tmp<0.05)/length(tmp);
                tmp = pvaluesI(:,test); tmp(isnan(tmp)) = [];
                irls(test,avg) = sum(tmp<0.05)/length(tmp);
            end
        end
        
        ols = ols.*100;
        irls = irls.*100;
        
        figure;
        subplot(1,2,1);
        scatter(ones(100,1),ols(1,:),30,'k'); hold on;
        scatter(ones(100,1).*2,irls(1,:),30,'k');
        for pair = 1:100
            plot([1,2],[ols(1,pair) irls(1,pair)]);
        end
        axis([0.8 2.2 4 6]); box on; grid on;
        [~,CI] = rst_1ttest([ols' irls'],'mean',0);
        title(sprintf('Type 1 error rate \n OLS %g [%g %g] \n IRLS %g [%g %g]', ...
            mean(ols(1,:)),CI(1,1),CI(2,1),mean(irls(1,:)),CI(1,5),CI(2,5)),'FontSize',12)
        ylabel('type 1 error rate')
        
        subplot(1,2,2);
        ols_irls_diff = (ols'-irls');
        [meandiff,cidiff]=rst_data_plot(ols_irls_diff(:,1),'estimator','mean','newfig','no');
        title(sprintf('Full factorial three-way ANOVA \n OLS - IRLS = %g \n CI=[%g %g]', ...
            meandiff,cidiff(1),cidiff(2)),'FontSize',12)
        ylabel('difference in type 1 error rate')
        
        % --------------------------------
        if sum(all(CI(:,1:8)<5))~=0 || sum(all(CI(:,1:8)>5))~=0
            fprintf('Full factorial three-way ANOVA OLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4),CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        else
            fprintf('Full factorial three-way ANOVA OLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,1),CI(2,1),CI(1,2),CI(2,2),CI(1,3),CI(2,3),CI(1,4),CI(2,4),CI(1,5),CI(2,5),CI(1,6),CI(2,6),CI(1,7),CI(2,7),CI(1,8),CI(2,8))
        end
        
        if sum(all(CI(:,9:16)<5))~=0 || sum(all(CI(:,9:16)>5))~=0
            fprintf('Full factorial three-way ANCOVA IRLS is not at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,9),CI(2,9),CI(1,10),CI(2,10),CI(1,11),CI(2,11),CI(1,12),CI(2,12), CI(1,13),CI(2,13),CI(1,14),CI(2,14),CI(1,15),CI(2,15),CI(1,16),CI(2,16))
            
        else
            fprintf('Full factorial three-way ANCOVA IRLS is at type 1 error 5%% \n [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] [%g %g] \n',CI(1,9),CI(2,9),CI(1,10),CI(2,10),CI(1,11),CI(2,11),CI(1,12),CI(2,12), CI(1,13),CI(2,13),CI(1,14),CI(2,14),CI(1,15),CI(2,15),CI(1,16),CI(2,16))
        end
        
        [~,cidiff] = rst_1ttest(ols_irls_diff,'mean',0);
        meandiff = mean(ols_irls_diff);
        
        nb_diff = max([sum(all(cidiff<0)) sum(all(cidiff>0))]);
        if nb_diff~=0
            if nb_diff == 1
                disp('Full factorial three-way ANCOVA: there is 1 significant difference in type 1 error between OLS and IRLS ')
            else
                fprintf('Full factorial three-way ANCOVA: there are %g significant differences in type 1 error between OLS and IRLS \n',nb_diff)
            end
        else
            disp('Full factorial three-way ANCOVA: there is no significant difference in type 1 error between OLS and IRLS')
        end
        
        T = array2table([mean(ols')' CI(:,1:8)' mean(irls')' CI(:,9:16)' meandiff' cidiff(1,:)' cidiff(2,:)'],...
            'VariableNames',{'MeanOLS','CIOlow','CIOhigh','MeanIRLS','CIIlow','CIIhigh' 'MeanDiff','CIDlow','CIDhigh'});
        disp(T)
        save('Full-factorial-three-way-ANOVA_table','T')
    end
end






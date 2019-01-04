function limotest_wls

% this is the unit test of the limo_glm function when calling limo_glm with
% the WLS option. Given null data in, only the nominal level of error should
% be obtained - to test that, we must create time vectors since this is a
% multivariate technique
% see also limotest_glm

nmc = 1000;

%% Simple Regression

Y = randn(100,20);
X = [randn(100,1) ones(100,1)];

wls  = NaN(1,100);
for avg = 1:100
    X2 = [randn(100,1) ones(100,1)];
    pvalues = NaN(20,nmc);
    
    parfor MC=1:nmc
        Y2 = randn(100,20);
        model = limo_glm(Y2,X2,0,0,1,'WLS','Time',0,1);
        pvalues(:,MC) = model.p;
    end
    pvalues = pvalues(:);
    pvalues(isnan(pvalues)) = [];
    wls(avg) = sum(pvalues<0.05)/length(pvalues);
end
wls = wls.*100;
[meanwls,ciwls]=rst_data_plot(wls,'estimator','mean','newfig','yes');

% one way ANOVA
wls = NaN(1,100);
X = [kron(eye(3),ones(30,1)) ones(90,1)];
for avg = 1:100
    pvalues = NaN(20,nmc);
    
    parfor MC=1:nmc
        Y2 = randn(90,20);
        model = limo_glm(Y2,X,1,0,0,'WLS','Time',0,20);
        pvalues(:,MC) = model.p;
    end
    pvalues = pvalues(:);
    pvalues(isnan(pvalues)) = [];
    wls(avg) = sum(pvalues<0.05)/length(pvalues);
end
wls = wls.*100;

[meanwls,ciwls]=rst_data_plot(wls,'estimator','mean','newfig','yes');



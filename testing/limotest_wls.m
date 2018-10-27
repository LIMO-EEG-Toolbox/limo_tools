function limotest_wls

% this is the unit test of the limo_glm function when calling limo_wls
% given null data in, only the nominal level of error should be obtained
% to test limo_wls, we must create time vectors since this is a
% multivariate technique, i.e. the outlier detection will be wrong
% see also limotest_glm

nmc = 1000;

%% Simple Regression

Y = randn(100,nmc);
X = [randn(100,1) ones(100,1)]; 
pvalues = NaN(1,nmc);
parfor MC=1:nmc
model = limo_glm(Y(:,MC),X,[],[],1,'OLS','Time',0,1);
pvalues(MC) = model.p;
end
pvalues(isnan(pvalues)) = [];
FP_rate = sum(pvalues<0.05)/length(pvalues);

% WLS
% ---
pvalues = NaN(1,nmc);
parfor MC=1:nmc
model = limo_glm(Y(:,MC),X,[],[],1,'WLS','Time',0,1);
pvalues(MC) = model.p;
end
pvalues(isnan(pvalues)) = [];
FP_rate = sum(pvalues<0.05)/length(pvalues);

% IRLS
% ----
pvalues = NaN(1,nmc);
parfor MC=1:nmc
model = limo_glm(Y(:,MC),X,[],[],1,'IRLS','Time',0,1);
pvalues(MC) = model.p;
end
pvalues(isnan(pvalues)) = [];


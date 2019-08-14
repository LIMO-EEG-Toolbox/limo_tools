% Testing frame for univariate statitistics in time with 1 weight per time series
% The weight is the multivariate distance from a time series to the others defined 
% by it's projection on the principal compoments that explain 90% of the data.

%% 1st let's check the distributions of stat values when no adjustment is performed

Rsquare_ols   = NaN(20,1000);
F_Rsquare_ols = NaN(20,1000);
p_Rsquare_ols = NaN(20,1000);
Rsquare_wls   = NaN(20,1000);
F_Rsquare_wls = NaN(20,1000);
p_Rsquare_wls = NaN(20,1000);

% draw 1000 random Y to regress X onto
parfor MC=1:10000
    Y = randn(100,20);
    X = [randn(100,1) ones(100,1)];
    Betas_ols = X\Y;
    [Betas_wls,W] = limo_WLS(X,Y);
    WX = [X(:,1:end-1).*repmat(W,1,size(X,2)-1) X(:,end)];
    
    T       = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total (the data)
    HM_ols  = X*pinv(X);                                                         % Hat matrix, projection onto X
    R_ols   = eye(size(Y,1)) - HM_ols;                                           % Projection onto E
    E_ols   = Y'*R_ols*Y;                                                        % SS Error
    HM_wls  = WX*pinv(WX);
    R_wls   = eye(size(Y,1)) - HM_wls;
    E_wls   = Y'*R_wls*Y;
    
    % degrees of freedom
    % -------------------
    % makes no difference ols or wls the difference in diagonals^2 is close to 0
    df(MC)  = trace(HM_ols'*HM_ols)^2/trace(HM_ols'*HM_ols*HM_ols'*HM_ols)-1; % Satterthwaite approximation
    dfe(MC) = trace((eye(size(HM_ols))-HM_ols)'*(eye(size(HM_ols))-HM_ols));
        
    % model R^2
    % -----------
    C                   = eye(size(X,2));
    C(:,size(X,2))      = 0;                              % all columns but the constant
    C0                  = eye(size(X,2)) - C*pinv(C);     % only the constant
    % OLS
    X0                  = X*C0;                           % Reduced model design matrix
    R0                  = eye(size(Y,1)) - (X0*pinv(X0)); % Projection onto error
    M                   = R0 - R_ols;                     % Projection matrix onto Xc
    H                   = (Betas_ols'*X'*M*X*Betas_ols);  % SS Effects
    Rsquare_ols(:,MC)   = diag(H)./diag(T);               % Variances explained
    F_Rsquare_ols(:,MC) = (diag(H)./df(MC)) ./ (diag(E_ols)/dfe(MC));
    p_Rsquare_ols(:,MC) = 1 - fcdf(F_Rsquare_ols(:,MC), df(MC), dfe(MC));
    % WLS
    X0                  = WX*C0;
    R0                  = eye(size(Y,1)) - (X0*pinv(X0));
    M                   = R0 - R_wls;
    H                   = (Betas_wls'*X'*M*X*Betas_wls);
    Rsquare_wls(:,MC)   = diag(H)./diag(T);
    F_Rsquare_wls(:,MC) = (diag(H)./df(MC)) ./ (diag(E_wls)/dfe(MC));
    p_Rsquare_wls(:,MC) = 1 - fcdf(F_Rsquare_wls(:,MC), df(MC), dfe(MC));
end
        
% vizualize and compare results for expectations under the null
[~,values]=rst_hist(F_Rsquare_ols(:));
Y = pdf('f',values,df(1),dfe(1));
figure; 
subplot(2,3,1); rst_hist(F_Rsquare_ols(:)); 
hold on; plot(values,Y,'r','LineWidth',2); 
axis([0 20 0 0.5]); title('F OLS')
subplot(2,3,2); rst_hist(F_Rsquare_wls(:)); 
hold on; plot(values,Y,'r','LineWidth',2); 
axis([0 20 0 0.5]); title('F WLS')
subplot(2,3,3); rst_density_hist(F_Rsquare_wls(:)-F_Rsquare_ols(:)); 
axis([-15 15 0 0.5]); title('WLS-OLS');

subplot(2,3,4); rst_hist(Rsquare_ols(:)); 
axis([0 0.12 0 60]); title('R OLS')
subplot(2,3,5); rst_hist(Rsquare_wls(:)); 
axis([0 0.12 0 60]); title('R WLS')
subplot(2,3,6); rst_density_hist(Rsquare_wls(:)-Rsquare_ols(:)); 
axis([-0.1 0.1 0 60]); title('WLS-OLS');

%% repeat adjusting for dfe

% draw 1000 random Y to regress X onto
parfor MC=1:10000
    Y = randn(100,20);
    X = [randn(100,8) ones(100,1)];
    Betas_ols = X\Y;
    [Betas_wls,W] = limo_WLS(X,Y);
    WX = [X(:,1:end-1).*repmat(W,1,size(X,2)-1) X(:,end)];
    
    T       = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total (the data)
    HM_ols  = X*pinv(X);                                                         % Hat matrix, projection onto X
    R_ols   = eye(size(Y,1)) - HM_ols;                                           % Projection onto E
    E_ols   = Y'*R_ols*Y;                                                        % SS Error
    HM_wls  = WX*pinv(WX);
    R_wls   = eye(size(Y,1)) - HM_wls;
    E_wls   = Y'*R_wls*Y;
    
    % degrees of freedom
    % -------------------
    % makes no difference ols or wls the difference in diagonals^2 is close to 0
    df(MC)  = trace(HM_ols'*HM_ols)^2/trace(HM_ols'*HM_ols*HM_ols'*HM_ols)-1; % Satterthwaite approximation
    dfe(MC) = trace((eye(size(HM_ols))-HM_ols)'*(eye(size(HM_ols))-HM_ols));
    % dfe_wls(MC) = size(Y,1)-size(Y,2)+df(MC);  
    
    p = size(Y,2);      % = number of variables (dimension)
    q = rank(WX)-1;      % -1 because of intercept column
    s = min(p,q);       % df
    n = size(Y,1);      % nb of observations (dfe)
    m = (abs(q-p)-1)/2;
    N = (n-q-p-2)/2;
    ve = n - rank(WX);
    df_wls(MC) = max(p,q); % (2*m+s+1); % 
    dfe_wls(MC) = ve-max(p,q)+q; % (2*N+s+1); % size(Y,1)-size(Y,2)-1
    
    
    % model R^2
    % -----------
    C                   = eye(size(X,2));
    C(:,size(X,2))      = 0;                              % all columns but the constant
    C0                  = eye(size(X,2)) - C*pinv(C);     % only the constant
    % OLS
    X0                  = X*C0;                           % Reduced model design matrix
    R0                  = eye(size(Y,1)) - (X0*pinv(X0)); % Projection onto error
    M                   = R0 - R_ols;                     % Projection matrix onto Xc
    H                   = (Betas_ols'*X'*M*X*Betas_ols);  % SS Effects
    Rsquare_ols(:,MC)   = diag(H)./diag(T);               % Variances explained
    F_Rsquare_ols(:,MC) = (diag(H)./df(MC)) ./ (diag(E_ols)/dfe(MC));
    p_Rsquare_ols(:,MC) = 1 - fcdf(F_Rsquare_ols(:,MC), df(MC), dfe(MC));
    % WLS
    X0                  = WX*C0;
    R0                  = eye(size(Y,1)) - (X0*pinv(X0));
    M                   = R0 - R_wls;
    H                   = (Betas_wls'*X'*M*X*Betas_wls);
    Rsquare_wls(:,MC)   = diag(H)./diag(T);
    F_Rsquare_wls(:,MC) = (diag(H)./df(MC)) ./ (diag(E_wls)/dfe_wls(MC));
    p_Rsquare_wls(:,MC) = 1 - fcdf(F_Rsquare_wls(:,MC), df(MC), dfe_wls(MC));
end
        
% vizualize and compare results for expectations under the null
[error_rates, CI]=rst_data_plot([mean(p_Rsquare_ols < .05,2) mean(p_Rsquare_wls < .05,2)].*100,'estimator','mean','bubble','on');

[~,values]=rst_hist(F_Rsquare_ols(:));
Y = pdf('f',values,df(1),dfe(1));
figure; 
subplot(2,3,1); rst_hist(F_Rsquare_ols(:)); 
hold on; plot(values,Y,'r','LineWidth',2); 
axis([0 20 0 0.5]); title('F OLS')
subplot(2,3,2); rst_hist(F_Rsquare_wls(:)); 
hold on; plot(values,Y,'r','LineWidth',2); 
axis([0 20 0 0.5]); title('F WLS')
subplot(2,3,3); rst_density_hist(F_Rsquare_wls(:)-F_Rsquare_ols(:)); 
axis([-15 15 0 0.5]); title('WLS-OLS');

subplot(2,3,4); rst_hist(Rsquare_ols(:)); 
axis([0 0.12 0 60]); title('R OLS')
subplot(2,3,5); rst_hist(Rsquare_wls(:)); 
axis([0 0.12 0 60]); title('R WLS')
subplot(2,3,6); rst_density_hist(Rsquare_wls(:)-Rsquare_ols(:)); 
axis([-0.1 0.1 0 60]); title('WLS-OLS');



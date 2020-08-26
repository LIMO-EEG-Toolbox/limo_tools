function [Difference,CI,SSanova,SSregression] = limo_non_linear(varargin)

% General Linear model for of EEG data
% The model consider trials as independent observations
% Two analyses are performed - a regression and an ANOVA
% the difference between sum of squares indicates non linear effects
%
% FORMAT:
% [Difference,CI,SSanova,SSregression] = limo_non_linear(Y,X1,X2,nb_conditions,nb_interactions,method)
%
% INPUTS: Y  = 2D matrix of EEG data with format trials x frames
%         X1 = 2 dimensional design matrix for the ANOVA model
%         X2 = 2 dimensional design matrix for the regression model
%         nb_conditions = a vector indicating the number of conditions per factor
%         nb_interactions = a vector indicating number of columns per interactions
%         method = 'OLS', 'WLS', 'IRLS' (bisquare)
%
% See also
% LIMO_DESIGN_MATRIX, LIMO_PCOUT, LIMO_IRLS
%
% Cyril Pernet v1 20-08-2012
% -----------------------------
%  Copyright (C) LIMO Team 2012

%% varagin
nboot = 600; % 600 resamples

if nargin == 6 
    y               = varargin{1};
    X1              = varargin{2};
    X2              = varargin{3};
    nb_conditions   = varargin{4};
    nb_interactions = varargin{5};
    method          = varargin{6};
    boot_table      = randi(size(y,1),size(y,1),nboot); 
else
    error('varargin error')
end

nb_factors = numel(nb_conditions);
if nb_factors == 1 && nb_conditions == 0
    nb_factors = 0;
end

% ----------- 
%% Data check
% -----------

if size(y,1)~=size(X1,1) || size(y,1)~=size(X2,1)
    error('The number of events in Y and the design matrix are different')
end

if nb_interactions == 0
    nb_interactions = [];
end

% ----------
%% Bootstrap
% -----------

h = waitbar(0,'bootstraping trials','name','% done');

for B = 1:nboot
    waitbar(B/nboot)
    
    % create data under H0 by breaking the link between Y and X
    Y = y(boot_table(:,B),:); % sample randomly from all Y
    
    % ------------------------------
    % Compute ANOVA model
    % ------------------------------
    
    % -------------------------
    if nb_factors == 1   %  1-way ANOVA
        % -------------------------
        
        % total sum of squares, projection matrix for errors, residuals and betas
        % -----------------------------------------------------------------------
        T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
        R     = eye(size(Y,1)) - (X1*pinv(X1));                                      % Projection on E
        E     = (Y'*R*Y);                                                          % SS Error
        
        % compute Beta parameters and weights
        if strcmp(method,'OLS')
            W = ones(size(Y,1),1);
            Betas = pinv(X1)*Y;
        elseif strcmp(method,'WLS')
            W = ones(size(Y,1),1);
            % use princ. comp. method for each condition
            for i = 1:nb_conditions
                index = find(X1(:,i));
                W(index,1)  = limo_pcout(Y(index,:));
            end
            WY = Y .* repmat(W,1,size(Y,2));
            WX = X1 .* repmat(W,1,size(X1,2));
            Betas = pinv(WX)*WY;
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X1,Y);
        end
        
        % compute model SS
        % -----------------
        C = eye(size(X1,2));
        C(:,size(X1,2)) = 0;
        C0 = eye(size(X1,2)) - C*pinv(C);
        X0 = X1*C0;  % Reduced model
        R0 = eye(size(Y,1)) - (X0*pinv(X0));
        M  = R0 - R;  % Projection matrix onto Xc
        SSanova(B,:)  = diag(Betas'*X1'*M*X1*Betas);  % SS Hypothesis (Effect)
        
        % ------------------------------------------------
    elseif nb_factors > 1  && isempty(nb_interactions) % N-ways ANOVA without interactions
        % ------------------------------------------------
        
        % compute basic SS total, projection matrices and parameters
        T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        R        = eye(size(Y,1)) - (X1*pinv(X1));
        E        = (Y'*R*Y);
        % compute Beta parameters and weights
        if strcmp(method,'OLS')
            W = ones(size(Y,1),1);
            Betas = pinv(X1)*Y;
        elseif strcmp(method,'WLS')
            W = ones(size(Y,1),1);
            % use princ. comp. method for each condition
            for i = 1:sum(nb_conditions)
                index = find(X1(:,i));
                W(index,1)  = limo_pcout(Y(index,:));
            end
            WY = Y .* repmat(W,1,size(Y,2));
            WX = X .* repmat(W,1,size(X1,2));
            Betas = pinv(WX)*WY;
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X1,Y);
        end
        
        % compute model SS
        % --------------------
        C = eye(size(X1,2));
        C(:,size(X1,2)) = 0;
        C0   = eye(size(X1,2)) - C*pinv(C);
        X0   = X1*C0; % Reduced model (i.e. only intercept)
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;      % M is the projection matrix onto Xc
        SSanova(B,:)    = diag(Betas'*X1'*M*X1*Betas);   % SS Hypothesis (Effect)
        
        % ------------------------------------------------
    elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways ANOVA with interactions
        % ------------------------------------------------
        
        % compute basic SS total, projection matrices and parameters
        T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        R        = eye(size(Y,1)) - (X1*pinv(X1));
        E        = (Y'*R*Y);
        if strcmp(method,'OLS')
            W = ones(size(Y,1),1);
            Betas = pinv(X1)*Y;
        elseif strcmp(method,'WLS')
            W = ones(size(Y,1),1);
            % use princ. comp. method for each condition,
            % ie higher order interaction
            for i = 1:nb_interaction(end)
                j = sum(nb_conditions)+i;
                index = find(X1(:,j));
                W(index,1)  = limo_pcout(Y(index,:));
            end
            WY = Y .* repmat(W,1,size(Y,2));
            WX = X1 .* repmat(W,1,size(X1,2));
            Betas = pinv(WX)*WY;
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X1,Y);
        end
        
        % compute model SS
        % --------------------
        C = eye(size(X1,2));
        C(:,size(X1,2)) = 0;
        C0   = eye(size(X1,2)) - C*pinv(C);
        X0   = X1*C0; % Reduced model (i.e. only intercept)
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;      % M is the projection matrix onto Xc
        SSanova(B,:)    = diag(Betas'*X1'*M*X1*Betas);   % SS Hypothesis (Effect)
        
    end
    
    
    % -----------------------------------
    %% Compute regression
    % -----------------------------------
    
            T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
            R     = eye(size(Y,1)) - (X2*pinv(X2));
            E     = (Y'*R*Y);
            if strcmp(method,'OLS')
                W = ones(size(Y,1),1);
                Betas = pinv(X2)*Y;
            elseif strcmp(method,'WLS')
                W = limo_pcout(Y);
                WY = Y .* repmat(W,1,size(Y,2));
                WX = X .* repmat(W,1,size(X2,2));
                Betas = pinv(WX)*WY;
            elseif strcmp(method,'IRLS')
                [Betas,W] = limo_IRLS(X2,Y);
            end
            
            % compute model R^2
            % -----------------
            C = eye(size(X2,2));
            C(:,size(X2,2)) = 0;
            C0 = eye(size(X2,2)) - C*pinv(C);
            X0 = X2*C0;
            R0 = eye(size(Y,1)) - (X0*pinv(X0));
            M  = R0 - R;
            SSregression(B,:)  = diag(Betas'*X2'*M*X2*Betas);
        
end
close(h)

%% now compare ANOVA / regression
Difference = mean(SSanova-SSregression);
Diff = sort(SSanova-SSregression);
low = round(nboot*5/100/2);
high = nboot - low; 
CI = [Diff(low,:); Diff(high,:)];
if sum(CI(1,:)>0)>0 || sum(CI(2,:)<0)>0
    disp('non linear effect detected')
end

% hist(SSanova(:,frame)-SSregression(:,frame)); grid on % is the histogram we want

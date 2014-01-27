function [F,p] = limo_robust_1way_anova(Y,X,percent)

% A heteroscedastic one-way ANOVA for trimmed means
% using a generalization of Welch's method.
% See Wilcox 2012: 
% Introduction to Robust Estimation and Hypothesis testing p 293
%
% FORMAT: [F,p] = limo_robust_1way_anova(Y,X)
%
% INPUT: Y is a 2D matrix times frames x trials/subjects
%        X is a design matrix of 1 and 0, one column per group
%           example with 3 groups of 3 subjects:
%           X=[1 1 1 0 0 0 0 0 0;0 0 0 1 1 1 0 0 0;0 0 0 0 0 0 1 1 1];
%        percent [0 100] is the amount of trimming - default 20
%        Warning: do not use this function to compare medians (percent=50)
%
% OUTPUT: F and p values for trimmed mean differences 
%
% EXAMPLE: 200 time frames x 5 groups of 20 subjects
%          Y = randn(200,100);
%          X = kron(eye(5),ones(20,1));
%          that is
%          X1 = [ones(20,1);zeros(80,1)];
%          X2 = [zeros(20,1);ones(20,1);zeros(60,1)];
%          X3 = [zeros(40,1);ones(20,1);zeros(40,1)];
%          X4 = [zeros(60,1);ones(20,1);zeros(20,1)];
%          X5 = [zeros(80,1);ones(20,1)];
%          X = [X1 X2 X3 X4 X5];
%
% Cyril Pernet v1 11/06/2013
% GAR 12/06/2013: edited help, fixed code and added percent input
% -----------------------------------------------------------------
%  Copyright (C) LIMO Team 2010

if nargin<3 || isempty(percent)
percent = 20;
end

J = size(X,2);
Nf = size(Y,1);
h = zeros(1,J); % same for all frames
w = zeros(Nf,J); % inverse of the adjusted variance 
xbar = cell(J); % trimmed means

% for each group, trim the data, get the sample size and winsorized variance
% compute w the inverse of the adjusted variance d
for gp = 1:J
    data = Y(:,X(:,gp)==1);
    na = size(data,2); % how many subjects
    ga=floor((percent/100)*na);% number of items to trim / winsorize
    asort=sort(data,2);
    trimdata=asort(:,(ga+1):(na-ga),:); % trimmed data
    h(gp) = size(trimdata,2); % effective sample size
    xbar{gp} = mean(trimdata,2);
    wa=asort; 
    wa(:,1:ga+1)=repmat(asort(:,ga+1),1, ga+1);
    wa(:,na-ga:end)=repmat(asort(:,na-ga),1, ga+1);
    wva=var(wa,0,2); % windsorized variances
    
    d = ((na-1).*wva) ./ (h(gp)*(h(gp)-1));
    w(:,gp) = 1./d;
end
    
U = sum(w,2);

% weighted marginal mean
% ----------------------
M=zeros(Nf,J);
for gp = 1:J
    M(:,gp) = sum(repmat(w(:,gp),1,size(xbar{gp},2)).*xbar{gp},2);
end
MM = sum(M,2) ./U;

% weighted variance from the marginal mean
% ----------------------------------------
AM=zeros(Nf,J);
for gp = 1:J
    AM(:,gp) = sum(repmat(w(:,gp),1,size(xbar{gp},2)).*(xbar{gp}-repmat(MM,1,size(xbar{gp},2))).^2,2);
end
A = sum(AM,2) ./ (J-1);

% error
% -------
f = 2*(J-2) / (J^2-1);  % reduces to 0 for J=2 = yuen t-test
B = f .* sum(((1 - (w./repmat(U,1,J))).^2) ./ repmat((h-1),size(Y,1),1),2);

% F and p values
% --------------
F = A./ (1+B);
df = J-1;
dfe = ((3/(J^2-1)) .* (sum(((1 - (w./repmat(U,1,J))).^2) ./ repmat((h-1),size(Y,1),1),2))).^-1;
p = 1 - fcdf(F, df, dfe);





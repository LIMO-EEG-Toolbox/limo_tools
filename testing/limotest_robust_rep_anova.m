function limotest_robust_rep_anova

% unit test under H0 for limo_robust_rep_anova
% beside testing for the input/output formats, it also tests that the
% outputs are valid (1) using 0% trimming compared to limo_rep_anova and
% (2) under H0. This allows ensuring valid F values and accurate control 
% the type 1 error rate 

%% One sample repeated measure
% ----------------------------

N             = 20;                                        % 20 subjects
gp            = ones(N,1);                                 % 1 group
factors       = 3;                                         % N repeated meaures
data          = randn(2,N,factors);                        % data 
result        = limo_rep_anova(data,gp,factors);           % standard ANOVA
robust_result = limo_robust_rep_anova(data,gp,factors,0);  % robust using 0% trimming
disp('one way repeated measure ANOVA')
compare_results(result,robust_result)

% test under H0, to check the type 1 error rate
nboot = 1000;
P1    = NaN(nboot,1);
P2    = NaN(nboot,1);
mu    = zeros(1,factors);
go = 0;
while go==0
    try
        SIGMA = rand(factors);
        SIGMA(SIGMA==diag(SIGMA)) = 1;
        SIGMA = SIGMA - tril(SIGMA,-1) + triu(SIGMA,1)';
        r = mvnrnd(mu,SIGMA,N); go=1;
    end
end

parfor n=1:nboot
    r = mvnrnd(mu,SIGMA,N); data = nan(2,N,factors);
    data(1,:,:) = r; data(2,:,:) = r;
    result        = limo_rep_anova(data,gp,factors);
    robust_result = limo_robust_rep_anova(data,gp,factors); % default 20% trimmed mean
    P1(n) = result.p(1);
    P2(n) = robust_result.p(1);
end

type1_error =  mean([P1 P2]<0.05);

end

%% output
function compare_results(result,robust_result)

if ~any(single(result.F)==single(robust_result.F))
    error('F values differ between standad and robust ANOVA with no trimming')
else
    disp('outputs for F values seem valid')
end

if ~any(single(result.p)==single(robust_result.p))
    error('p values differ between standad and robust ANOVA with no trimming')
else
    disp('outputs for p values seem valid')
end

end


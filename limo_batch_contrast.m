function limo_batch_contrast(limo_file,C)

% routine to compute T contrasts from the batch mode
% -----------------------------
% Copyright (C) LIMO Team 2014

% Cyril Pernet June 2014

load(limo_file);
cd(LIMO.dir)

%% 1st check the contrasts 
disp('contrast checking ..')
for j=1:size(C,1)
    out(j,:) = limo_contrast_checking(LIMO.dir, LIMO.design.X, C(j,:));
    go = limo_contrast_checking(out(j,:),LIMO.design.X);
    if go == 0
        error(sprintf('contrast %g is invalid',j));
    end
end

%% do the analysis

fprintf('loading data - and computing contrast(s) in %s \n',pwd)

load Yr; load Betas;
if isfield(LIMO,'contrast')
    previous_con = size(LIMO.contrast,2);
else
    previous_con = 0;
end

for i=1:size(C,1)  % for each new contrast
    LIMO.contrast{previous_con+i}.V = 'T';
    LIMO.contrast{previous_con+i}.C = out(i,:); save LIMO LIMO
    limo_contrast(Yr, Betas, LIMO, 0,1);
end



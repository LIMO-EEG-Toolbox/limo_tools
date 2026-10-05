function s = limo_glm_boot_chunk(Yr, LIMO, X, channel, boot_table)

% One CHANNEL of the GLM bootstrap under H0 (non time-frequency): runs
% limo_glm_boot once for this channel (all nboot, exactly as limo_glm_handling
% does) and returns that channel's H0 sub-arrays. The bootstrap is chunked by
% channel rather than by bootstrap because limo_glm_boot re-draws the null
% (limo_glm_null uses randperm) on every call, so it must be called exactly
% once per channel to match the non-chunked result.
%
% Returns a struct s with per-channel arrays (only the fields present in the
% design), with the channel dimension dropped:
%   s.R2    [frames x 3 x nboot]
%   s.Betas [frames x nparam x nboot]
%   s.Cond  [frames x ncond x 2 x nboot]   (if conditions)
%   s.Int   [frames x nint  x 2 x nboot]   (if full factorial)
%   s.Cov   [frames x ncont x 2 x nboot]   (if continuous regressors)
% ------------------------------
%  Copyright (C) LIMO Team 2026

nframes = size(Yr,2); nparam = size(X,2);
method  = LIMO.design.method;
nboot   = LIMO.design.bootstrap;
hascond = prod(LIMO.design.nb_conditions) ~= 0;
hasint  = (LIMO.design.fullfactorial == 1);
hascov  = (LIMO.design.nb_continuous ~= 0);

% ---- per-channel limo_glm_boot (all nboot), exactly as limo_glm_handling ----
if LIMO.Level == 2
    Y = squeeze(Yr(channel,:,:));
    index = find(~isnan(Y(1,:)));
    if strcmp(method,'WLS') || strcmp(method,'OLS')
        Weights = squeeze(LIMO.design.weights(channel,index))';
    else
        Weights = squeeze(LIMO.design.weights(channel,:,index));
    end
    model = limo_glm_boot(squeeze(Y(:,index))', X(index,:), Weights, ...
        LIMO.design.nb_conditions, LIMO.design.nb_interactions, LIMO.design.nb_continuous, ...
        method, LIMO.Analysis, boot_table{channel});
else % LIMO.Level == 1
    L = LIMO;
    if strcmp(method,'WLS') || strcmp(method,'OLS')
        L.Weights = squeeze(LIMO.design.weights(channel,:))';
    else
        L.Weights = squeeze(LIMO.design.weights(channel,:,:));
    end
    model = limo_glm_boot(squeeze(Yr(channel,:,:))', L, boot_table);
end

% ---- write-back into per-channel arrays (identical to limo_glm_handling) ----
s.R2    = NaN(nframes,3,nboot);
s.Betas = NaN(nframes,nparam,nboot);
if hascond, s.Cond = NaN(nframes,length(LIMO.design.nb_conditions),2,nboot); end
if hasint,  s.Int  = NaN(nframes,length(LIMO.design.nb_interactions),2,nboot); end
if hascov,  s.Cov  = NaN(nframes,LIMO.design.nb_continuous,2,nboot); end

for B = 1:nboot
    s.Betas(:,:,B) = model.betas{B};
    s.R2(:,1,B)    = model.R2_univariate{B};
    s.R2(:,2,B)    = model.F{B};
    s.R2(:,3,B)    = model.p{B};

    if hascond
        if length(LIMO.design.nb_conditions) == 1
            s.Cond(:,1,1,B) = model.conditions.F{B};
            s.Cond(:,1,2,B) = model.conditions.p{B};
        else
            for i = 1:length(LIMO.design.nb_conditions)
                s.Cond(:,i,1,B) = model.conditions.F{B}(i,:);
                s.Cond(:,i,2,B) = model.conditions.p{B}(i,:);
            end
        end
    end

    if hasint
        if length(LIMO.design.nb_interactions) == 1
            s.Int(:,1,1,B) = model.interactions.F{B};
            s.Int(:,1,2,B) = model.interactions.p{B};
        else
            for i = 1:length(LIMO.design.nb_interactions)
                s.Int(:,i,1,B) = model.interactions.F{B}(i,:);
                s.Int(:,i,2,B) = model.interactions.p{B}(i,:);
            end
        end
    end

    if hascov
        if LIMO.design.nb_continuous == 1
            s.Cov(:,1,1,B) = model.continuous.F{B};
            s.Cov(:,1,2,B) = model.continuous.p{B};
        else
            for i = 1:LIMO.design.nb_continuous
                if all(size(squeeze(s.Cov(:,i,1,B))) == size(squeeze(model.continuous.F{B}(:,i)))) || ...
                        all(size(squeeze(s.Cov(:,i,1,B))) == size(squeeze(model.continuous.F{B}(:,i))'))
                    s.Cov(:,i,1,B) = model.continuous.F{B}(:,i);
                    s.Cov(:,i,2,B) = model.continuous.p{B}(:,i);
                else
                    s.Cov(:,i,1,B) = model.continuous.F{B}(i,:);
                    s.Cov(:,i,2,B) = model.continuous.p{B}(i,:);
                end
            end
        end
    end
end
end

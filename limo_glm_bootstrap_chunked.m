function limo_glm_bootstrap_chunked(LIMO, Yr, X, array, boot_table, nboot, subname)

% Channel-streamed, resumable, memory-bounded GLM bootstrap under H0 (non
% time-frequency). Each channel is bootstrapped with a single limo_glm_boot
% call (exactly as limo_glm_handling), checkpointed to disk, and its slice
% written into the H0 files with a low-memory matfile writer -- so the full
% [channels x frames x ... x nboot] arrays are never held in RAM. The
% bootstrap is chunked by CHANNEL (not by bootstrap) because limo_glm_boot
% re-draws the null on every call (limo_glm_null uses randperm); calling it
% once per channel reproduces the non-chunked behaviour, and a crashed run
% can resume from the last completed channel.
%
% Writes the same H0 files as the non-chunked path: R2 and Betas, plus one
% file per condition / interaction / continuous effect.
% ------------------------------
%  Copyright (C) LIMO Team 2026

nchan = size(Yr,1); nframes = size(Yr,2); nparam = size(X,2);
H0dir     = fullfile(LIMO.dir,'H0');
chunk_dir = fullfile(H0dir,'glm_chunks');
if ~exist(chunk_dir,'dir'), mkdir(chunk_dir); end

% reuse the saved boot_table if present (so a resumed run keeps the same
% resampling and the completed channel chunks stay valid)
btfile = fullfile(H0dir,[subname 'boot_table.mat']);
if exist(btfile,'file')
    bt_saved   = load(btfile); boot_table = bt_saved.boot_table;
else
    save(btfile,'boot_table');
end

hascond = prod(LIMO.design.nb_conditions) ~= 0;
hasint  = (LIMO.design.fullfactorial == 1);
hascov  = (LIMO.design.nb_continuous ~= 0);

isvalid = false(nchan,1); isvalid(array) = true;
ndone   = sum(arrayfun(@(c) exist(local_cp(chunk_dir,c),'file')>0, 1:nchan));
if ndone>0, fprintf('  resuming GLM bootstrap: %d/%d channels on disk\n', ndone, nchan); end

% --- compute the missing channels (one limo_glm_boot call each) ---
for ch = 1:nchan
    if exist(local_cp(chunk_dir,ch),'file'), continue; end
    if isvalid(ch)
        fprintf('  bootstrapping channel %g/%g\n', ch, nchan);
        s = limo_glm_boot_chunk(Yr, LIMO, X, ch, boot_table);
    else
        s = struct();
    end
    s.valid = isvalid(ch); %#ok<STRNU>
    save(local_cp(chunk_dir,ch), '-struct', 's', '-v7.3');
    clear s
end

% --- output filenames (level 1 uses the subname prefix) ---
if LIMO.Level == 1
    f_R2 = fullfile(H0dir,[subname 'R2H0.mat']);
    f_B  = fullfile(H0dir,[subname 'BetasH0.mat']);
else
    f_R2 = fullfile(H0dir,'R2_desc-H0.mat');
    f_B  = fullfile(H0dir,'Betas_desc-H0.mat');
end
if exist(f_R2,'file'), delete(f_R2); end
if exist(f_B,'file'),  delete(f_B);  end
mR2 = matfile(f_R2,'Writable',true);
mB  = matfile(f_B ,'Writable',true);

% per-effect files
condf = {}; intf = {}; covf = {};
if hascond
    nE = length(LIMO.design.nb_conditions); condf = cell(1,nE); condm = cell(1,nE);
    for i=1:nE
        if LIMO.Level==1, nm=sprintf('%sCondition_effect_%gH0',subname,i); else, nm=sprintf('Condition_effect_%g_desc-H0',i); end
        condf{i}=fullfile(H0dir,[nm '.mat']); if exist(condf{i},'file'), delete(condf{i}); end
        condm{i}=matfile(condf{i},'Writable',true);
    end
end
if hasint
    nE = length(LIMO.design.nb_interactions); intf = cell(1,nE); intm = cell(1,nE);
    for i=1:nE
        if LIMO.Level==1, nm=sprintf('%sInteraction_effect_%gH0',subname,i); else, nm=sprintf('Interaction_effect_%g_desc-H0',i); end
        intf{i}=fullfile(H0dir,[nm '.mat']); if exist(intf{i},'file'), delete(intf{i}); end
        intm{i}=matfile(intf{i},'Writable',true);
    end
end
if hascov
    nE = LIMO.design.nb_continuous; covf = cell(1,nE); covm = cell(1,nE);
    for i=1:nE
        if LIMO.Level==1, nm=sprintf('%sdesc-Covariate_effect_%gH0',subname,i); else, nm=sprintf('Covariate_effect_%g_desc-H0',i); end
        covf{i}=fullfile(H0dir,[nm '.mat']); if exist(covf{i},'file'), delete(covf{i}); end
        covm{i}=matfile(covf{i},'Writable',true);
    end
end

% --- merge channels into the H0 files (one channel in RAM at a time) ---
for ch = 1:nchan
    S = load(local_cp(chunk_dir,ch));
    if S.valid
        r2 = reshape(S.R2,    [1 nframes 3 nboot]);
        be = reshape(S.Betas, [1 nframes nparam nboot]);
    else
        r2 = NaN(1,nframes,3,nboot);
        be = NaN(1,nframes,nparam,nboot);
    end
    if ch==1, mR2.H0_R2 = r2; mB.H0_Betas = be;
    else,     mR2.H0_R2(ch,:,:,:) = r2; mB.H0_Betas(ch,:,:,:) = be; end

    if hascond, local_weffect(condm, S, 'Cond', 'H0_Condition_effect',        ch, S.valid, nframes, nboot); end
    if hasint,  local_weffect(intm,  S, 'Int',  'H0_Interaction_effect',       ch, S.valid, nframes, nboot); end
    if hascov,  local_weffect(covm,  S, 'Cov',  'H0_Covariate_effect',         ch, S.valid, nframes, nboot); end
    clear S
end
end

% -------------------------------------------------------------------------
function local_weffect(mc, S, field, var, ch, valid, nframes, nboot)
for i = 1:numel(mc)
    if valid, slice = reshape(S.(field)(:,i,:,:), [1 nframes 2 nboot]);
    else,     slice = NaN(1,nframes,2,nboot); end
    if ch==1, mc{i}.(var) = slice; else, mc{i}.(var)(ch,:,:,:) = slice; end
end
end

function p = local_cp(d,c), p = fullfile(d, sprintf('glm_chunk_%04d.mat', c)); end

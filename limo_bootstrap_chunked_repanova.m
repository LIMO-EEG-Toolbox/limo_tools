function limo_bootstrap_chunked_repanova(H0dir, dims, type, nC, b_compute, rep_files, irep_files, gp_file, opts)

% Chunked, resumable, memory-bounded bootstrap for the repeated-measures ANOVA
% under H0. Bootstraps are processed in chunks (one chunk in memory), each
% chunk checkpointed to disk, then merged into the per-effect H0 files with a
% low-memory matfile writer. Identical result to the non-chunked accumulation
% for the same boot_table.
%
% H0dir      the analysis H0 directory
% dims       [nchan nframes nboot]
% type       1..4 (limo repeated-measures design type)
% nC         number of repeated-measures contrasts (size of the effect dim)
% b_compute  handle: [mainSub,gpSub,intSub] = b_compute(B) for one bootstrap B
% rep_files  cellstr (length nC) of output files for the main effects
%            (variable H0_Rep_ANOVA, [nchan nframes 2 nboot])
% irep_files cellstr (length nC) of interaction-with-group files (type 3/4),
%            else {} (variable H0_Rep_ANOVA_Interaction_with_gp)
% gp_file    group-effect file (type 3/4), else '' (variable H0_Rep_ANOVA_Gp_effect)
% opts       .chunk_size, .chunk_dir (optional)
% ------------------------------
%  Copyright (C) LIMO Team 2026

nchan = dims(1); nframes = dims(2); nboot = dims(3);
if nargin < 9 || isempty(opts), opts = struct; end
withgp = (type == 3 || type == 4);

if ~isfield(opts,'chunk_dir') || isempty(opts.chunk_dir), opts.chunk_dir = fullfile(H0dir,'rep_chunks'); end
if ~isfield(opts,'chunk_size') || isempty(opts.chunk_size)
    bytes_per_boot = nchan*nframes*nC*2*8*(1 + 2*withgp);
    opts.chunk_size = max(1, floor(512*1024^2 / max(bytes_per_boot,1)));
end
opts.chunk_size = max(1, min(round(opts.chunk_size), nboot));
if ~exist(opts.chunk_dir,'dir'), mkdir(opts.chunk_dir); end

starts   = 1:opts.chunk_size:nboot;
ranges   = arrayfun(@(s) s:min(s+opts.chunk_size-1,nboot), starts, 'UniformOutput', false);
n_chunks = numel(ranges);
done     = false(1,n_chunks);
for c = 1:n_chunks, if exist(local_cpath(opts.chunk_dir,c),'file'), done(c) = true; end, end
if any(done), fprintf('  resuming rep-ANOVA bootstrap: %d/%d chunks on disk\n', nnz(done), n_chunks); end

% --- compute the missing chunks ---
for c = find(~done)
    br = ranges{c}; nbc = numel(br);
    fprintf('  [chunk %d/%d] bootstraps %d-%d\n', c, n_chunks, br(1), br(end));
    mains = cell(1,nbc); gps = cell(1,nbc); ints = cell(1,nbc);
    parfor bi = 1:nbc
        [mains{bi}, gps{bi}, ints{bi}] = b_compute(br(bi));
    end
    cMain = NaN(nchan,nframes,nC,2,nbc);
    if withgp, cGp = NaN(nchan,nframes,2,nbc); cInt = NaN(nchan,nframes,nC,2,nbc); else, cGp = []; cInt = []; end
    for bi = 1:nbc
        cMain(:,:,:,:,bi) = mains{bi};
        if withgp, cGp(:,:,:,bi) = gps{bi}; cInt(:,:,:,:,bi) = ints{bi}; end
    end
    save(local_cpath(opts.chunk_dir,c), 'cMain','cGp','cInt','br','-v7.3');
    clear cMain cGp cInt mains gps ints
end

% --- merge chunks into the per-effect H0 files (one chunk in RAM at a time) ---
mRep = cell(1,nC); mInt = cell(1,nC);
for i = 1:nC
    f = local_ext(rep_files{i}); if exist(f,'file'), delete(f); end
    mRep{i} = matfile(f,'Writable',true);
    if withgp
        f = local_ext(irep_files{i}); if exist(f,'file'), delete(f); end
        mInt{i} = matfile(f,'Writable',true);
    end
end
if withgp
    f = local_ext(gp_file); if exist(f,'file'), delete(f); end
    mGp = matfile(f,'Writable',true);
end

for c = 1:n_chunks
    S = load(local_cpath(opts.chunk_dir,c)); br = S.br; nbc = numel(br);
    for i = 1:nC
        slice = reshape(S.cMain(:,:,i,:,:), [nchan nframes 2 nbc]);
        if c == 1, mRep{i}.H0_Rep_ANOVA = slice; else, mRep{i}.H0_Rep_ANOVA(:,:,:,br) = slice; end
        if withgp
            islice = reshape(S.cInt(:,:,i,:,:), [nchan nframes 2 nbc]);
            if c == 1, mInt{i}.H0_Rep_ANOVA_Interaction_with_gp = islice;
            else,      mInt{i}.H0_Rep_ANOVA_Interaction_with_gp(:,:,:,br) = islice; end
        end
    end
    if withgp
        if c == 1, mGp.H0_Rep_ANOVA_Gp_effect = S.cGp; else, mGp.H0_Rep_ANOVA_Gp_effect(:,:,:,br) = S.cGp; end
    end
    clear S
end
end

% -------------------------------------------------------------------------
function p = local_cpath(d,c), p = fullfile(d, sprintf('rep_chunk_%04d.mat', c)); end
function f = local_ext(f), if ~endsWith(f,'.mat'), f = [f '.mat']; end, end

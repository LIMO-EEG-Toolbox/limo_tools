function limo_bootstrap_chunked(H0_file, var_name, dims, array, nboot, chanfun, opts)

% Memory-bounded, resumable, extendable bootstrap accumulation for LIMO.
%
% Computes a bootstrap H0 array of size [nchan nframes 2 nboot] by processing
% the bootstrap dimension in CHUNKS: each chunk is computed, written to disk,
% and only one chunk is ever held in memory. The chunks are then merged into
% the final H0 file with a low-memory (matfile) writer. Because each
% H0(channel,:,:,b) is a deterministic function of the (centred) data and the
% resampling indices boot_table{channel}(:,b), the merged result is IDENTICAL
% to the non-chunked limo_random_robust bootstrap when the same boot_table is
% used -- chunking only changes the order/granularity of accumulation, never
% the numbers.
%
% Why: for high-density data (e.g. 256 channels x 525 frames x 1000 boots the
% full H0 array is ~2 GB) the non-chunked path holds the whole array in RAM
% and saves once at the very end, so it is memory-heavy and a crash loses all
% progress. This routine bounds memory to one chunk and checkpoints every
% chunk, so a run can be resumed from the last completed chunk after a crash,
% and more bootstraps can be appended later.
%
% FORMAT limo_bootstrap_chunked(H0_file, var_name, dims, array, nboot, chanfun)
%        limo_bootstrap_chunked(..., opts)
%
% INPUTS
%   H0_file  full path to the final H0 .mat to write (its folder must exist)
%   var_name name of the variable to store, e.g. 'H0_paired_samples'
%   dims     [nchan nframes 2 nboot] final dimensions
%   array    vector of channel indices that are actually computed; other
%            channels stay NaN (matching limo_random_robust)
%   nboot    total number of bootstraps (= dims(4))
%   chanfun  function handle [tmat,pmat] = chanfun(channel, b_range) returning
%            the t and p maps for the requested bootstrap columns, each of
%            size [nframes x numel(b_range)]. (chanfun owns the per-channel
%            data extraction and the parfor over bootstraps, so its result is
%            byte-for-byte the limo_random_robust per-iteration computation.)
%   opts     optional struct:
%            .chunk_size  bootstraps per chunk (default: auto, ~512 MB/chunk)
%            .chunk_dir   chunk/checkpoint folder (default: <H0dir>/chunks)
%            .keep_chunks keep chunk files after merge (default: true)
%            .verbose     print progress (default: true)
%
% See also LIMO_RANDOM_ROBUST, LIMO_CREATE_BOOT_TABLE
% ------------------------------
%  Copyright (C) LIMO Team 2026

if nargin < 7 || isempty(opts), opts = struct; end
nchan = dims(1); nframes = dims(2); nstat = dims(3);

H0_dir = fileparts(H0_file);
if ~isfield(opts,'chunk_dir')  || isempty(opts.chunk_dir),  opts.chunk_dir  = fullfile(H0_dir,'chunks'); end
if ~isfield(opts,'keep_chunks'),                            opts.keep_chunks = true;  end
if ~isfield(opts,'verbose'),                                opts.verbose     = true;  end
if ~isfield(opts,'chunk_size') || isempty(opts.chunk_size)
    bytes_per_boot = nchan*nframes*nstat*8;
    opts.chunk_size = max(1, floor(512*1024^2 / max(bytes_per_boot,1)));
end
opts.chunk_size = max(1, min(round(opts.chunk_size), nboot));

if ~exist(opts.chunk_dir,'dir'), mkdir(opts.chunk_dir); end
meta_file = fullfile(opts.chunk_dir,[var_name '_chunks_meta.mat']);

% contiguous chunk plan over 1:nboot
starts   = 1:opts.chunk_size:nboot;
b_ranges = arrayfun(@(s) s:min(s+opts.chunk_size-1,nboot), starts, 'UniformOutput', false);
n_chunks = numel(b_ranges);

% resume: a chunk counts as done if its file exists on disk
done = false(1,n_chunks);
for c = 1:n_chunks
    if exist(local_chunk_path(opts.chunk_dir,var_name,c),'file'), done(c) = true; end
end
if opts.verbose && any(done)
    fprintf('  resuming: %d/%d chunks already on disk\n', nnz(done), n_chunks);
end

% compute the missing chunks
for c = find(~done)
    b_range = b_ranges{c};
    if opts.verbose
        fprintf('  [chunk %d/%d] bootstraps %d-%d (size %d)\n', ...
            c, n_chunks, b_range(1), b_range(end), numel(b_range));
    end
    H0_chunk = NaN(nchan, nframes, nstat, numel(b_range));
    for e = 1:numel(array)
        channel = array(e);
        [tmat,pmat] = chanfun(channel, b_range);          % [nframes x nb]
        H0_chunk(channel,:,1,:) = reshape(tmat, 1, nframes, 1, numel(b_range));
        H0_chunk(channel,:,2,:) = reshape(pmat, 1, nframes, 1, numel(b_range));
    end
    save(local_chunk_path(opts.chunk_dir,var_name,c), 'H0_chunk', 'b_range', '-v7.3');
    chunk_size = opts.chunk_size; %#ok<NASGU>
    save(meta_file, 'n_chunks', 'b_ranges', 'chunk_size', 'nboot', '-v7.3');
    clear H0_chunk
end

% merge chunks -> final H0 file (one chunk in RAM at a time)
if opts.verbose, fprintf('  merging %d chunks -> %s\n', n_chunks, H0_file); end
if exist(H0_file,'file'), delete(H0_file); end
m = matfile(H0_file,'Writable',true);
for c = 1:n_chunks
    S = load(local_chunk_path(opts.chunk_dir,var_name,c));   % H0_chunk, b_range
    if c == 1
        m.(var_name) = S.H0_chunk;                           % creates+sizes the variable
    else
        m.(var_name)(:,:,:,S.b_range) = S.H0_chunk;          % grows contiguously
    end
    clear S
end

if ~opts.keep_chunks
    for c = 1:n_chunks, delete(local_chunk_path(opts.chunk_dir,var_name,c)); end
    if exist(meta_file,'file'), delete(meta_file); end
end
end

% -------------------------------------------------------------------------
function p = local_chunk_path(d,v,c)
p = fullfile(d, sprintf('%s_chunk_%04d.mat', v, c));
end

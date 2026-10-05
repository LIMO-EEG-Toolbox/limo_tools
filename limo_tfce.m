function [tfce_score,thresholded_maps] = limo_tfce(varargin)

% LIMO_TFCE Threshold-Free Cluster Enhancement (Smith & Nichols, 2009),
% computed EXACTLY via the eTFCE algorithm (default) or with the classic
% discretised approximation (method = 'classic').
%
% The TFCE integral
%
%       TFCE(v) = integral_{0}^{h_v} extent(h)^E * h^H dh
%
% is, by default, evaluated EXACTLY rather than approximated by a fixed number
% of equally-spaced thresholds (the 'dh' loop of the classic method). Since
% extent(h) is piecewise constant in h (it changes only when h crosses a data
% value), the integral has a closed form on each interval and is obtained in a
% single pass using a disjoint-set (union-find) cluster-retrieval forest:
% supra-threshold voxels are processed in non-increasing order, each voxel's
% cluster extent at its own height is the union-find component size, and
% TFCE(v) is accumulated along each voxel's path to its forest root. This is
% exact (no step-size parameter, no discretisation bias), is typically ~50x
% faster than the classic dh loop, and needs no Image Processing / SPM
% clustering routine.
%
% ATTRIBUTION
% -----------
% The exact-TFCE (eTFCE) algorithm implemented here was NOT developed by the
% LIMO team. It is the method of:
%
%   Chen, X., Weeda, W. D., Nichols, T. E., & Goeman, J. J. (2026).
%   eTFCE: Exact Threshold-Free Cluster Enhancement via Fast Cluster
%   Retrieval. arXiv:2603.03004.   https://arxiv.org/abs/2603.03004
%
% building on the cluster-retrieval algorithm of Chen, Goeman, Krebs, Meijer
% & Weeda (2023). This function is only a MATLAB implementation of that
% published method for LIMO; please cite the paper above when using the
% exact (default) option.
%
% INPUT  tfce_score = limo_tfce(type,data,channeighbstructmat)
%        tfce_score = limo_tfce(type,data,channeighbstructmat,updatebar,E,H,dh,minnbchan,method)
%
%        type 1 = 1D (one channel ERP/Power), 2 = 2D (channels x time, or a
%             single freq x time map), 3 = 3D (ERSP, channels x freq x time).
%             The first dimension must be channels.
%        data  = a single map of t/F values, or a set of maps under H0
%             (stacked along the last dimension).
%        channeighbstructmat = channel neighbourhood matrix (n_chan x n_chan,
%             1 = neighbours). If empty for type 2, a 4-connected freq x time
%             grid is used.
%        updatebar = flag (default 1) to show a waitbar for stacks.
%        E,H   = TFCE exponents (defaults 0.5 and 2).
%        dh    = step size. Used ONLY by method = 'classic'; ignored by the
%             exact method (which has no step size).
%        minnbchan = minimum number of supra-threshold CHANNEL neighbours a
%             point must have to be kept (FieldTrip criterion), iterated to a
%             fixed point before clustering. DEFAULT 2, which reproduces the
%             current LIMO cluster definition (limo_findcluster is called with
%             minnbchan = 2). Set 0 for plain connected components, i.e. the
%             canonical Smith & Nichols / eTFCE-paper TFCE. Applied for type 2
%             and type 3 with a non-empty neighbour matrix; ignored for type 1
%             and empty-neighbour maps (matching the classic code).
%        method = 'exact' (default) for eTFCE, or 'classic' to dispatch to the
%             original discretised algorithm (limo_tfce_classic) for exact
%             reproducibility of pre-eTFCE results.
%
% OUTPUT tfce_score = map of TFCE scores (same spatial dims as the input;
%             bootstrap dimension preserved for stacks).
%        thresholded_maps = the per-dh maps for method = 'classic'; [] for the
%             exact method (which has no discrete dh levels). It is only used
%             as a diagnostic by limo_tfce_handling.
%
% NOTES
%   * minnbchan = 2 (default) keeps the SAME cluster definition as the classic
%     LIMO TFCE, so the only change from the classic result is the exact vs
%     discretised integral. On real data the two agree very closely
%     (correlation ~0.999); a few boundary voxels can differ because the
%     classic integral is a <=200-step approximation. Use method = 'classic'
%     to reproduce previous results bit-for-bit.
%   * Signed maps are split into positive and negative parts, each scored and
%     the two summed -- same convention as the classic implementation.
%   * Both the observed map and the H0 maps go through the same transform, so
%     permutation/bootstrap inference is internally consistent.
%
% References
%   Smith, S. M., & Nichols, T. E. (2009). Threshold-free cluster
%     enhancement. NeuroImage, 44(1), 83-98.
%   Pernet, C., Latinus, M., Nichols, T. E., & Rousselet, G. A. (2015).
%     Cluster-based computational methods for mass univariate analyses of
%     event-related brain potentials/fields. J Neurosci Methods, 250, 85-93.
%   Chen, Weeda, Nichols & Goeman (2026), arXiv:2603.03004 (eTFCE; see above).
%
% See also LIMO_TFCE_CLASSIC, LIMO_TFCE_HANDLING, LIMO_FINDCLUSTER
% ------------------------------
%  eTFCE engine added 2026 (exact method: Chen, Weeda, Nichols & Goeman 2026)
%  Copyright (C) LIMO Team 2026

%% check / parse input
if nargin < 3
    error('limo_tfce:notEnoughArgs','not enough arguments')
elseif nargin > 9
    error('limo_tfce:tooManyArgs','too many arguments')
end

type                = varargin{1};
data                = varargin{2};
channeighbstructmat = varargin{3};

if nargin >= 4 && ~isempty(varargin{4}); updatebar = varargin{4}; else; updatebar = 1;       end
if nargin >= 5 && ~isempty(varargin{5}); E         = varargin{5}; else; E         = 0.5;     end
if nargin >= 6 && ~isempty(varargin{6}); H         = varargin{6}; else; H         = 2;       end
if nargin >= 7 && ~isempty(varargin{7}); dh        = varargin{7}; else; dh        = 0.1;     end
if nargin >= 8 && ~isempty(varargin{8}); minnbchan = varargin{8}; else; minnbchan = 2;       end
if nargin >= 9 && ~isempty(varargin{9}); method    = varargin{9}; else; method    = 'exact'; end
% dh is used ONLY by method = 'classic'; the exact method has no step size.
% minnbchan DEFAULT 2 reproduces the classic LIMO cluster definition; set 0
% for the canonical no-pruning (Smith & Nichols / eTFCE-paper) TFCE. See header.

% Dispatch to the original discretised algorithm for exact reproducibility.
if strcmpi(method,'classic')
    [tfce_score,thresholded_maps] = limo_tfce_classic(type,data,channeighbstructmat,updatebar,E,H,dh);
    return
end

thresholded_maps = [];
if isempty(channeighbstructmat); channeighbstructmat = []; end  % normalise []

%% determine single-map size (sz) and number of stacked maps (nb)
switch type
    case 1
        if isvector(data)
            sz = [1 numel(data)]; nb = 1;
        else
            [x,nb] = size(data); sz = [1 x];
        end
    case 2
        if size(data,3) <= 1
            sz = [size(data,1) size(data,2)]; nb = 1;
        else
            sz = [size(data,1) size(data,2)]; nb = size(data,3);
        end
    case 3
        if size(data,4) <= 1
            sz = [size(data,1) size(data,2) size(data,3)]; nb = 1;
        else
            sz = [size(data,1) size(data,2) size(data,3)]; nb = size(data,4);
        end
    otherwise
        error('limo_tfce:badType','type must be 1, 2 or 3')
end

Nlin = prod(sz);

%% optional minnbchan pruning context (k-core entry-height reduction)
% Applied per spatial slice on the channel graph, exactly where the classic
% code prunes: type 2 (slices = time points) and type 3 (slices = freq*time),
% both only with a non-empty neighbour matrix. Silently skipped otherwise
% (type 1, or empty neighbours), matching classic which does not prune those.
prune = [];
if minnbchan > 0 && (type == 2 || type == 3) && ~isempty(channeighbstructmat)
    prune = struct('k', minnbchan, 'neighb', channeighbstructmat, ...
                   'C', sz(1), 'nslice', Nlin/sz(1));
end

%% build the (value-independent) adjacency once for this geometry
edges = local_build_edges(type, sz, channeighbstructmat);

%% flatten the data to [Nlin x nb] and score each map
dataflat = reshape(double(data), Nlin, nb);
out      = zeros(Nlin, nb);

showbar = (updatebar == 1) && (nb > 1);
if showbar
    try, f = waitbar(0,'percentage of maps analyzed','name','eTFCE'); catch, showbar = false; end
end

for boot = 1:nb
    out(:,boot) = local_etfce_map(dataflat(:,boot), edges, Nlin, E, H, prune);
    if showbar; try, waitbar(boot/nb); end; end %#ok<TRYNC>
end

if showbar; try, close(f); end; end %#ok<TRYNC>

%% reshape to expected output layout (trailing singleton dropped for nb==1)
tfce_score = reshape(out, [sz nb]);

end % limo_tfce


% =====================================================================
%  Core: exact TFCE for one map (handles signed data via pos/neg split)
% =====================================================================
function score = local_etfce_map(M, edges, Nlin, E, H, prune)

if nargin < 6, prune = []; end
M(~isfinite(M)) = 0;                 % NaN/Inf -> background (non supra-threshold)
pospart = max( M, 0);
negpart = max(-M, 0);

% Optional minnbchan pruning: replace each part by its per-time-slice k-core
% "entry height" h_eff. Running standard eTFCE on h_eff is exactly equivalent
% to integrating TFCE over the minnbchan-pruned supra-threshold sets, because
% the k-core is monotone in the threshold (pruned set at level h == {h_eff>h}).
if ~isempty(prune)
    pospart = reshape(local_compute_heff(reshape(pospart, prune.C, prune.nslice), prune.neighb, prune.k), [], 1);
    negpart = reshape(local_compute_heff(reshape(negpart, prune.C, prune.nslice), prune.neighb, prune.k), [], 1);
end

score = local_etfce_pass(pospart, edges, Nlin, E, H) ...
      + local_etfce_pass(negpart, edges, Nlin, E, H);

end


% =====================================================================
%  minnbchan reduction: per-time-slice k-core "entry height".
%  For each voxel (channel c, time t) returns h_eff(c,t) = the largest
%  threshold h at which c is still in the k-core (>= k supra-threshold
%  channel neighbours, iterated) of the channels with value > h at time t.
%  h_eff = 0 for voxels never in the k-core. Equivalent to the iterated
%  FieldTrip/limo_findcluster minnbchan deletion, evaluated at every level.
% =====================================================================
function heff = local_compute_heff(V, neighb, k)

[C,T] = size(V);
A = (neighb ~= 0); A = A | A';            % symmetric channel adjacency
A(1:C+1:end) = false;                      % no self-loops
nbrs = cell(C,1);
for c = 1:C, nbrs{c} = find(A(:,c))'; end  % row vectors of neighbour ids

heff = zeros(C,T);
for t = 1:T
    v   = V(:,t);
    on  = v > 0;
    if ~any(on), continue; end

    deg = zeros(C,1);
    onv = find(on)';
    for c = onv, deg(c) = sum(on(nbrs{c})); end   % degree within the on-set

    removed = ~on;                 % off channels are "removed" background
    hloc    = zeros(C,1);
    current_h = 0;

    % peel order: on-channels by ascending value (threshold removals)
    [~,o] = sort(v(onv),'ascend');
    ord   = onv(o);
    ptr   = 1;

    % cascade queue: on-channels already below the k-neighbour requirement
    q  = find(on & deg < k)';
    qi = 1;

    while true
        % drain degree-cascade (these leave the core at the current height)
        while qi <= numel(q)
            x = q(qi); qi = qi + 1;
            if removed(x), continue; end
            removed(x) = true; hloc(x) = current_h;
            for y = nbrs{x}
                if ~removed(y)
                    deg(y) = deg(y) - 1;
                    if deg(y) < k, q(end+1) = y; end %#ok<AGROW>
                end
            end
        end
        % next threshold removal (lowest-value surviving channel)
        while ptr <= numel(ord) && removed(ord(ptr)), ptr = ptr + 1; end
        if ptr > numel(ord), break; end
        x = ord(ptr); ptr = ptr + 1;
        current_h = v(x);
        removed(x) = true; hloc(x) = current_h;
        for y = nbrs{x}
            if ~removed(y)
                deg(y) = deg(y) - 1;
                if deg(y) < k, q(end+1) = y; end %#ok<AGROW>
            end
        end
    end
    heff(:,t) = hloc;
end
end


% =====================================================================
%  Exact TFCE for the positive part of one map.
%  vals >= 0; only strictly positive entries are nodes (supra-threshold
%  for some h > 0). Integration lower limit is h0 = 0.
% =====================================================================
function T_full = local_etfce_pass(vals, edges, Nlin, E, H)

T_full = zeros(Nlin,1);
active = vals > 0;
N = nnz(active);
if N == 0, return; end

% map active linear indices -> node ids 1..N
activeIdx          = find(active);
nodeid             = zeros(Nlin,1);
nodeid(activeIdx)  = 1:N;
nodeval            = vals(activeIdx);          % node statistic values (>0)

% keep only edges whose two endpoints are both active, in node-id space
ea   = nodeid(edges(:,1));
eb   = nodeid(edges(:,2));
keep = (ea > 0) & (eb > 0);
ea   = ea(keep);
eb   = eb(keep);

% neighbour adjacency in CSR form (store both directions)
src = [ea; eb];
dst = [eb; ea];
[src, ord_e] = sort(src);
dst          = dst(ord_e);
cnt      = accumarray(src, 1, [N 1]);          % #neighbours per node (0 if none)
startp   = cumsum([1; cnt]);                    % node i -> dst(startp(i):startp(i+1)-1)

% process nodes in non-increasing value order (rank 1 = largest)
[~, order] = sort(nodeval, 'descend');
rnk        = zeros(N,1);
rnk(order) = 1:N;

% union-find state
ufp   = (1:N)';      % disjoint-set parent (path-compressed; for fast find)
csize = ones(N,1);   % component size keyed by current representative
croot = (1:N)';      % hierarchy root of a component, keyed by representative
hpar  = (1:N)';      % forest (hierarchy) parent of each node; self == root
ext   = ones(N,1);   % e_v(h_v): cluster extent of node when it is added

for p = 1:N
    i  = order(p);
    s  = startp(i);
    e  = startp(i+1) - 1;

    % representative of i's (currently singleton) component
    ri = i; while ufp(ri) ~= ri, ri = ufp(ri); end
    a = i; while ufp(a) ~= ri, t = ufp(a); ufp(a) = ri; a = t; end   % compress

    for k = s:e
        j = dst(k);
        if rnk(j) >= p, continue; end           % neighbour not yet processed
        rj = j; while ufp(rj) ~= rj, rj = ufp(rj); end
        a = j; while ufp(a) ~= rj, t = ufp(a); ufp(a) = rj; a = t; end  % compress
        if rj == ri, continue; end              % already same component

        cr        = croot(rj);                  % current root of absorbed subtree
        hpar(cr)  = i;                          % i absorbs it -> i becomes its parent

        if csize(ri) >= csize(rj)               % union by size
            ufp(rj)    = ri;
            csize(ri)  = csize(ri) + csize(rj);
            croot(ri)  = i;
        else
            ufp(ri)    = rj;
            csize(rj)  = csize(ri) + csize(rj);
            croot(rj)  = i;
            ri         = rj;                     % representative changed
        end
    end

    ext(i) = csize(ri);                          % extent of i's cluster at height nodeval(i)
end

% Algorithm 2: walk each node's path to its forest root (reverse rank order
% guarantees a node's parent is computed before the node itself).
T = zeros(N,1);
c = 1/(H+1);
for p = N:-1:1
    i = order(p);
    u = hpar(i);
    if u == i
        T(i) = c * (ext(i)^E) * (nodeval(i)^(H+1));
    else
        T(i) = T(u) + c * (ext(i)^E) * (nodeval(i)^(H+1) - nodeval(u)^(H+1));
    end
end

T_full(activeIdx) = T;

end


% =====================================================================
%  Adjacency builder (geometry only, independent of the data values).
%  Returns an [M x 2] list of linear-index pairs into a single map of
%  size sz, matching the connectivity used by limo_findcluster:
%    type 1            : 1D temporal chain
%    type 2 (+neighb)  : per-channel temporal chain + same-time channel edges
%    type 2 (empty)    : 4-connected 2D (freq*time) grid
%    type 3 (+neighb)  : per-channel 4-conn (freq*time) + same-(f,t) channel edges
%    type 3 (empty)    : per-channel 4-conn (freq*time), channels independent
% =====================================================================
function edges = local_build_edges(type, sz, channeighbstructmat)

edges     = zeros(0,2);
hasNeighb = ~isempty(channeighbstructmat);

switch type

    case 1
        T = prod(sz);
        if T >= 2
            edges = [(1:T-1)', (2:T)'];
        end

    case 2
        if hasNeighb
            C = sz(1); Y = sz(2);
            te = zeros(0,2); ce = zeros(0,2);
            % per-channel temporal chain: (c,t)-(c,t+1)
            if Y >= 2
                [cc,tt] = ndgrid(1:C, 1:Y-1);
                a = cc(:) + (tt(:)-1)*C;
                te = [a, a + C];
            end
            % same-time channel edges from the neighbour matrix
            [c1,c2] = local_neighbour_pairs(channeighbstructmat);
            if ~isempty(c1)
                [pp,tt] = ndgrid(1:numel(c1), 1:Y);
                off = (tt(:)-1)*C;
                ce  = [c1(pp(:)) + off, c2(pp(:)) + off];
            end
            edges = [te; ce];
        else
            % single freq*time map: 4-connected grid
            edges = local_grid4(sz(1), sz(2));
        end

    case 3
        C = sz(1); F = sz(2); Tt = sz(3);
        we = zeros(0,2);
        % per-channel frequency edges: (c,f,t)-(c,f+1,t)
        if F >= 2
            [cc,ff,tt] = ndgrid(1:C, 1:F-1, 1:Tt);
            a  = cc(:) + (ff(:)-1)*C + (tt(:)-1)*C*F;
            we = [we; a, a + C];
        end
        % per-channel temporal edges: (c,f,t)-(c,f,t+1)
        if Tt >= 2
            [cc,ff,tt] = ndgrid(1:C, 1:F, 1:Tt-1);
            a  = cc(:) + (ff(:)-1)*C + (tt(:)-1)*C*F;
            we = [we; a, a + C*F];
        end
        % same-(f,t) channel edges from the neighbour matrix
        if hasNeighb
            [c1,c2] = local_neighbour_pairs(channeighbstructmat);
            if ~isempty(c1)
                [pp,ff,tt] = ndgrid(1:numel(c1), 1:F, 1:Tt);
                off = (ff(:)-1)*C + (tt(:)-1)*C*F;
                we  = [we; c1(pp(:)) + off, c2(pp(:)) + off];
            end
        end
        edges = we;
end

end


function [c1,c2] = local_neighbour_pairs(channeighbstructmat)
% unique unordered channel pairs flagged as neighbours
nb        = channeighbstructmat ~= 0;
nb        = triu(nb | nb', 1);           % symmetrise, keep upper triangle
[c1, c2]  = find(nb);
c1 = c1(:); c2 = c2(:);
end


function E = local_grid4(X, Y)
% 4-connected adjacency on an X-by-Y grid (column-major linear indices)
E = zeros(0,2);
if X >= 2
    [xx,yy] = ndgrid(1:X-1, 1:Y);
    a = xx(:) + (yy(:)-1)*X;
    E = [E; a, a + 1];
end
if Y >= 2
    [xx,yy] = ndgrid(1:X, 1:Y-1);
    a = xx(:) + (yy(:)-1)*X;
    E = [E; a, a + X];
end
end

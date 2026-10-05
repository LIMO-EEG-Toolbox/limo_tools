function test_chunked_bootstrap
% Self-contained correctness tests for the chunked LIMO bootstrap.
%
%   * equivalence : the chunked H0 is IDENTICAL to the original (non-chunked)
%                   accumulation for the same boot_table -- proven against an
%                   inline copy of the previous limo_random_robust loop.
%   * granularity : the result does not depend on chunk size (cases 1, 2, 3).
%   * resume      : after losing chunk files, only the missing chunks are
%                   recomputed and the merged H0 is unchanged.
%   * add-more    : extending nboot leaves the earlier bootstraps untouched.
%
% Uses the LIMO default method (Trimmed Mean). Requires the LIMO toolbox on
% the path (limo_random_robust, limo_yuend_ttest, limo_ttest, limo_yuen_ttest,
% limo_trimci, limo_trimmed_mean, limo_create_boot_table).

rng(7);
nchan = 12; nframes = 30; nboot = 60; meth = 'Trimmed Mean';
allok = true;

% ---- CASE 3 (paired): equivalence to the original non-chunked loop ----
d1 = randn(nchan,nframes,18); d2 = d1 + 0.25*randn(nchan,nframes,18);
dc = tempname; mkdir(dc);
L  = mk(dc,nboot,meth); L.design.bootstrap_chunk = 20;     % -> 3 chunks
limo_random_robust(3, d1, d2, [1 2], L);
bn3 = 'Paired_Samples_Ttest_parameter_1_2_desc-H0';
Hc  = ld(fullfile(dc,'H0',[bn3 '.mat']),'H0_paired_samples');
bt  = ld(fullfile(dc,'H0','boot_table.mat'),'boot_table');
array = intersect(find(~isnan(d1(:,1,1))), find(~isnan(d2(:,1,1))));
Href = ref_paired(d1, d2, bt, array, nboot, true);          % inline original algorithm
allok = report('case3 chunked == original non-chunked loop', isequaln(Hc,Href), maxdiff(Hc,Href)) & allok;

% resume: drop the merged H0 + 2 of 3 chunk files, rerun
delete(fullfile(dc,'H0',[bn3 '.mat']));
cf = dir(fullfile(dc,'H0','chunks','H0_paired_samples_chunk_*.mat'));
for k = 1:min(2,numel(cf)), delete(fullfile(cf(k).folder,cf(k).name)); end
limo_random_robust(3, d1, d2, [1 2], L);
Hr = ld(fullfile(dc,'H0',[bn3 '.mat']),'H0_paired_samples');
allok = report('case3 resume after losing chunks == full', isequaln(Hc,Hr), 0) & allok;

% add-more: extend 60 -> 100, earlier bootstraps unchanged
L2 = mk(dc,100,meth); L2.design.bootstrap_chunk = 20;
limo_random_robust(3, d1, d2, [1 2], L2);
Hm = ld(fullfile(dc,'H0',[bn3 '.mat']),'H0_paired_samples');
allok = report('case3 add-more 60->100 prefix unchanged', size(Hm,4)==100 && isequaln(Hm(:,:,:,1:nboot),Hr), 0) & allok;

% ---- granularity: chunk size must not change the result (cases 1,2,3) ----
allok = gran_paired(d1,d2,nboot,meth) & allok;
allok = gran_onesample(randn(nchan,nframes,20),nboot,meth) & allok;
allok = gran_twosample(randn(nchan,nframes,20),randn(nchan,nframes,17),nboot,meth) & allok;

fprintf('-----------------------------------------------\n');
if allok, disp('ALL PASS'); else, error('test_chunked_bootstrap: FAIL'); end
end

% ---------- inline reference = previous (non-chunked) paired loop ----------
function H0 = ref_paired(data1, data2, boot_table, array, nboot, trimmed)
if trimmed
    d1c = data1 - repmat(limo_trimmed_mean(data1),[1 1 size(data1,3)]);
    d2c = data2 - repmat(limo_trimmed_mean(data2),[1 1 size(data2,3)]);
else
    d1c = data1 - repmat(nanmean(data1,3),[1 1 size(data1,3)]);
    d2c = data2 - repmat(nanmean(data2,3),[1 1 size(data2,3)]);
end
H0 = NaN(size(data1,1), size(data1,2), 2, nboot);
for e = 1:numel(array)
    ch = array(e);
    tmp = d1c(ch,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:))));
    tmp = d2c(ch,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:))));
    for b = 1:nboot
        if trimmed
            [t,~,~,~,p,~,~] = limo_yuend_ttest(Y1(1,:,boot_table{ch}(:,b)),Y2(1,:,boot_table{ch}(:,b)));
        else
            [~,~,~,~,~,t,p] = limo_ttest(1,Y1(1,:,boot_table{ch}(:,b)),Y2(1,:,boot_table{ch}(:,b)));
        end
        H0(ch,:,1,b) = t; H0(ch,:,2,b) = p;
    end
end
end

% ---------- granularity checks (chunk=small vs chunk=nboot) ----------
function ok = gran_paired(d1,d2,nboot,meth)
A = runchunk(@(L)limo_random_robust(3,d1,d2,[1 2],L), 'Paired_Samples_Ttest_parameter_1_2_desc-H0','H0_paired_samples',nboot,meth,17,   4243);
B = runchunk(@(L)limo_random_robust(3,d1,d2,[1 2],L), 'Paired_Samples_Ttest_parameter_1_2_desc-H0','H0_paired_samples',nboot,meth,nboot,4243);
ok = report('case3 chunk-size invariance', isequaln(A,B), maxdiff(A,B));
end
function ok = gran_onesample(d,nboot,meth)
A = runchunk(@(L)limo_random_robust(1,d,1,L), 'One_Sample_Ttest_parameter_1_desc-H0','H0_one_sample',nboot,meth,17,   4244);
B = runchunk(@(L)limo_random_robust(1,d,1,L), 'One_Sample_Ttest_parameter_1_desc-H0','H0_one_sample',nboot,meth,nboot,4244);
ok = report('case1 chunk-size invariance', isequaln(A,B), maxdiff(A,B));
end
function ok = gran_twosample(g1,g2,nboot,meth)
A = runchunk(@(L)limo_random_robust(2,g1,g2,[1 2],L), 'Two_Samples_Ttest_parameter_1_2_desc-H0','H0_two_samples',nboot,meth,17,   4245);
B = runchunk(@(L)limo_random_robust(2,g1,g2,[1 2],L), 'Two_Samples_Ttest_parameter_1_2_desc-H0','H0_two_samples',nboot,meth,nboot,4245);
ok = report('case2 chunk-size invariance', isequaln(A,B), maxdiff(A,B));
end
function H = runchunk(callfun, boot_name, var, nboot, meth, csz, seed)
% Seed the RNG identically for the small-chunk and full-chunk runs so both build
% the SAME boot_table (limo_create_boot_table draws from the global RNG, once,
% before chunking). With a shared seed the two runs compare the same resampling
% at different chunk sizes -- which is exactly the invariance being tested.
% (Without this, each run drew an independent unseeded boot_table, so the check
% compared two different resamplings and spuriously failed.)
if nargin >= 7 && ~isempty(seed), rng(seed); end
d = tempname; mkdir(d); L = mk(d,nboot,meth); L.design.bootstrap_chunk = csz;
callfun(L); H = ld(fullfile(d,'H0',[boot_name '.mat']),var);
end

% ---------- helpers ----------
function L = mk(d,nboot,meth)
L = struct('dir',d,'Analysis','Time','Type','Channels');
L.design.bootstrap = nboot; L.design.tfce = 0; L.design.method = meth;
end
function v = ld(f,name), S = load(f); v = S.(name); end
function d = maxdiff(a,b), m=~isnan(a)&~isnan(b); if any(m(:)), d=max(abs(a(m)-b(m))); else, d=0; end, end
function ok = report(nm,ok,md)
if ok, fprintf('  [PASS] %-42s (max|diff|=%.2e)\n',nm,md);
else,  fprintf('  [FAIL] %-42s (max|diff|=%.2e)\n',nm,md); end
end

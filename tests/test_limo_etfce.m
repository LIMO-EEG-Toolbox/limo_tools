function test_limo_etfce
% Validation of the exact eTFCE limo_tfce.m against an INDEPENDENT
% brute-force exact integrator.
%
% The reference recomputes, at every distinct data value a_k, the connected
% components of {val >= a_k} with MATLAB's graph/conncomp and adds the exact
% closed-form interval contribution
%       (1/(H+1)) * extent^E * (a_k^(H+1) - a_{k-1}^(H+1))
% to every node in the component. This shares no code with the union-find
% forest in limo_tfce.m, so agreement is a genuine cross-check.
%
% Run:  matlab -batch "cd('/Users/devon7y/Downloads/limo_tools-3.4'); test_limo_etfce"

rng(7);
E = 0.5; H = 2; tol = 1e-9;
fprintf('=== eTFCE validation (forest vs independent brute-force exact) ===\n');

npass = 0; nfail = 0;

% ---- type 1: 1D signed ERP ----------------------------------------------
T = 40;
m = randn(1,T)*2;
A = adj_chain(T);
[npass,nfail] = check('type1 signed', limo_tfce(1,m,[]), refmap(m(:),A,E,H), tol, npass, nfail);

% ---- type 2 channels x time, with channel neighbours --------------------
C = 12; Y = 30;
neighb = rand_neighb(C);
m = randn(C,Y)*1.5;
A = adj_type2(C,Y,neighb);
% minnbchan=0 -> paper-faithful no-pruning path (validated against plain CC ref)
[npass,nfail] = check('type2 signed +neighb', limo_tfce(2,m,neighb,0,0.5,2,0.1,0), reshape(refmap(m(:),A,E,H),C,Y), tol, npass, nfail);

% ---- type 2 all-positive (F-like) ---------------------------------------
mf = abs(randn(C,Y))*3 + 0.1;
[npass,nfail] = check('type2 F-like (>0)', limo_tfce(2,mf,neighb,0,0.5,2,0.1,0), reshape(refmap(mf(:),A,E,H),C,Y), tol, npass, nfail);

% ---- type 2 empty neighbours -> 4-connected freq*time grid --------------
X = 10; Yg = 16;
mg = randn(X,Yg)*1.2;
Ag = adj_grid4(X,Yg);
[npass,nfail] = check('type2 grid (empty neighb)', limo_tfce(2,mg,[]), reshape(refmap(mg(:),Ag,E,H),X,Yg), tol, npass, nfail);

% ---- type 3 channels x freq x time, with channel neighbours -------------
C3 = 8; F = 6; Tt = 10;
neighb3 = rand_neighb(C3);
m3 = randn(C3,F,Tt)*1.3;
A3 = adj_type3(C3,F,Tt,neighb3);
[npass,nfail] = check('type3 signed +neighb', limo_tfce(3,m3,neighb3,0,0.5,2,0.1,0), reshape(refmap(m3(:),A3,E,H),C3,F,Tt), tol, npass, nfail);

% ---- stack consistency: each slice == single-map call (default minnbchan) --
stack = randn(C,Y,4)*1.5;
sc = limo_tfce(2, stack, neighb, 0);
okstack = true;
for b = 1:4
    s1 = limo_tfce(2, stack(:,:,b), neighb, 0);
    if max(abs(s1(:) - reshape(sc(:,:,b),[],1))) > tol, okstack = false; end
end
[npass,nfail] = report('stack==per-slice', okstack, npass, nfail);

% ---- non-negativity & zeros where expected ------------------------------
sc2 = limo_tfce(2,m,neighb,0,0.5,2,0.1,0);
[npass,nfail] = report('scores finite & >=0', all(isfinite(sc2(:))) && all(sc2(:) >= -tol), npass, nfail);

fprintf('-----------------------------------------------\n');
fprintf('PASSED %d / %d tests\n', npass, npass+nfail);
if nfail > 0
    error('test_limo_etfce: %d test(s) FAILED', nfail);
end
end

% ====================== helpers ======================
function [np,nf] = check(name, got, want, tol, np, nf)
got = got(:); want = want(:);
d = max(abs(got - want));
scale = max(1, max(abs(want)));
ok = d/scale < tol;
[np,nf] = report(sprintf('%-26s (max|err|=%.2e)', name, d), ok, np, nf);
end

function [np,nf] = report(name, ok, np, nf)
if ok, fprintf('  [PASS] %s\n', name); np = np+1;
else,  fprintf('  [FAIL] %s\n', name); nf = nf+1; end
end

function T = refmap(vals, A, E, H)
% independent exact integrator on positive + negative parts
vals(~isfinite(vals)) = 0;
T = refpass(max(vals,0),A,E,H) + refpass(max(-vals,0),A,E,H);
end

function T = refpass(vals, A, E, H)
N = numel(vals); T = zeros(N,1);
av = sort(unique(vals(vals>0)), 'ascend');
prev = 0; c = 1/(H+1);
for k = 1:numel(av)
    a = av(k);
    idx = find(vals >= a);
    if ~isempty(idx)
        comp = conncomp(graph(A(idx,idx)));
        cs   = accumarray(comp(:),1);
        csz  = cs(comp(:));
        T(idx) = T(idx) + c * (csz.^E) * (a^(H+1) - prev^(H+1));
    end
    prev = a;
end
end

function A = adj_chain(T)
i = (1:T-1)'; j = (2:T)';
A = sparse([i;j],[j;i],1,T,T); A = double(A>0);
end

function A = adj_grid4(X,Y)
ii=[];jj=[]; lin=@(x,y) x+(y-1)*X;
for y=1:Y, x=(1:X-1)'; ii=[ii;lin(x,y)]; jj=[jj;lin(x+1,y)]; end
for y=1:Y-1, x=(1:X)'; ii=[ii;lin(x,y)]; jj=[jj;lin(x,y+1)]; end
A=sparse([ii;jj],[jj;ii],1,X*Y,X*Y); A=double(A>0);
end

function A = adj_type2(C,Y,neighb)
ii=[];jj=[]; lin=@(c,t) c+(t-1)*C;
for t=1:Y-1, c=(1:C)'; ii=[ii;lin(c,t)]; jj=[jj;lin(c,t+1)]; end
nb=triu((neighb|neighb'),1); [p,q]=find(nb);
for t=1:Y, ii=[ii;lin(p,t)]; jj=[jj;lin(q,t)]; end
A=sparse([ii;jj],[jj;ii],1,C*Y,C*Y); A=double(A>0);
end

function A = adj_type3(C,F,Tt,neighb)
ii=[];jj=[]; lin=@(c,f,t) c+(f-1)*C+(t-1)*C*F;
for t=1:Tt, for f=1:F-1, c=(1:C)'; ii=[ii;lin(c,f,t)]; jj=[jj;lin(c,f+1,t)]; end, end
for t=1:Tt-1, for f=1:F, c=(1:C)'; ii=[ii;lin(c,f,t)]; jj=[jj;lin(c,f,t+1)]; end, end
nb=triu((neighb|neighb'),1); [p,q]=find(nb);
for t=1:Tt, for f=1:F, ii=[ii;lin(p,f,t)]; jj=[jj;lin(q,f,t)]; end, end
A=sparse([ii;jj],[jj;ii],1,C*F*Tt,C*F*Tt); A=double(A>0);
end

function nb = rand_neighb(C)
% random symmetric channel adjacency with no isolated structure issues
nb = rand(C) < 0.3; nb = triu(nb,1); nb = nb | nb';
end

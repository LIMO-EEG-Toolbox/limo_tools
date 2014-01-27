function [eigen_vectors,eigen_values] = limo_decomp(E,H)

% this function is used to decompose inv(E)*H
%
% Following Rencher 2002 (Methods of multivariate
% analysis - Wiley) we note that eig(inv(E)*H) =
% eig((E^1/2)*H*inv(E^1/2)) = eig(inv(U')*H*inv(U))
%
% E^1/2 is the square root matrix of E and U'U = E
% (Cholesky factorization). Using the Cholesky 
% factorisation, we return positve eigen values
% from inv(U')*H*inv(U) which is positive semidefinite.
% If this procedre fails (E is not positive definite) 
% we then use an SVD decomposition
%
% Cyril Pernet v2 29-05-2009
% -----------------------------
%  Copyright (C) LIMO Team 2010

try
    U = chol(E);
    [eigen_vectors D] = eig(inv(U')*H*inv(U));
    eigen_values = diag(D);

catch
    y = (pinv(E)*H);
    [m n]   = size(y);
    if m > n
        [v s v] = svd(y*y');
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
        eigen_vectors = v;
    else
        [u s u] = svd(y'*y);
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
        eigen_vectors = u;
    end
    d  = sign(sum(v)); u = u*d;
    eigen_values  = u*sqrt(s(1)/n);
end

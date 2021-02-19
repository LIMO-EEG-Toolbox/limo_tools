function [eigen_vectors,eigen_values] = limo_decomp(E,H,type)

% FORMAT: [eigen_vectors,eigen_values] = limo_decomp(E,H,type)
%
% INPUT E and H are matrices, typically square symmetric Sum of Squares and
%       Cross Products
%       type is 'Chol' (default) or 'SVD')
%
% OUTPUT: the eigen vectors and values of the decomposition of inv(E)*H
%
% Following Rencher 2002 (Methods of multivariate analysis - Wiley) we note
% that eig(inv(E)*H) = % eig((E^1/2)*H*inv(E^1/2)) = eig(inv(U')*H*inv(U))
% and E^1/2 is the square root matrix of E and U'U = E (Cholesky factorization).
% Using the Cholesky factorisation, we return positve eigen values from
% inv(U')*H*inv(U) which is positive semidefinite. If this procedre fails
% (E is not positive definite) we then use an eigen value decomposition of pinv(E)*H
% It is also possible to procede using an SVD decomposition using the argument
% type ('SVD')
%
% Cyril Pernet 2009
% Cyril Pernet and Iege Bassez 2017
% -----------------------------
%  Copyright (C) LIMO Team 2010

% check input
if nargin < 2
    error('not enough arguments in')
elseif nargin == 2
    type = 'Chol';
end

% proceede
if strcmpi(type,'chol')
    try
        U = chol(E);
        [eigen_vectors, D] = eig(inv(U')*H*inv(U));
        eigen_values = diag(D);
        
    catch
        [eigen_vectors, eigen_values] = eig((pinv(E)*H));
    end
    
    
elseif strcmpi(type,'SVD')
    y = (pinv(E)*H);
    [m, n]   = size(y);
    if m > n
        [v,s,v] = svd(y*y');
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
        eigen_vectors = v;
    else
        [u, s,u] = svd(y'*y);
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
        eigen_vectors = u;
    end
    d  = sign(sum(v)); u = u*d;
    eigen_values  = u*sqrt(s(1)/n);
end


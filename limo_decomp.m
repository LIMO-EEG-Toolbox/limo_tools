function [eigen_vectors,eigen_values] = limo_decomp(varargin)

% FORMAT: [eigen_vectors,eigen_values] = limo_decomp(E,H,method)
%
% INPUT E and H are matrices, typically square symmetric Sum of Squares and
%               Cross Products
%       method is the decomposition method for pinv(E)*H
%              - pseudo (default)
%              - Cholesky
%              - SVD
%
% OUTPUT: the eigen vectors and values of the decomposition of inv(E)*H
%
% It is often the case that E is not positive definite, thus the default is 
% to use the eigen value decomposition: eig((pinv(E)*H))
% It is however also possible to use two other methods: 
% - Cholesky factorization eig(inv(E)*H) = eig((E^1/2)*H*inv(E^1/2)) = eig(inv(U')*H*inv(U))
%   --> returns positve eigen values from inv(U')*H*inv(U)
% - SVD to find first the eigen vectors and then their values
%
% Cyril Pernet & Iege Bassez
% -----------------------------
%  Copyright (C) LIMO Team 2018

%% inputs
E = varargin{1};
H = varargin{2};

if nargin == 2
    cov_method = 'pseudo'; % the default
else 
    cov_method = varargin{3};
end

% decompose
if strcmpi(cov_method, 'pseudo')
    [vec, D] = eig((pinv(E)*H));

    % sort eigenvalues and then sort eigenvectors in order of decreasing eigenvalues
    [e,ei] = sort(diag(D)); 
    ordered_eigenvalues = flipud(e);
    vec = vec(:,flipud(ei));

    % validate if correct eigenvalues and eigenvectors of matrix pinv(E)*H:
    if round((pinv(E)*H) * vec, 4) == round(vec * diag(ordered_eigenvalues), 4)
        eigen_vectors = vec;
        eigen_values = ordered_eigenvalues;
    else
        error('limo_decomp could not decompose the cross-product matrix (find eigenvalues/vectors), try regularized covariance if not done so already');
    end
    
elseif strcmp(method, 'cholesky')
        U = chol(E);
        [eigen_vectors, D] = eig(inv(U')*H*inv(U));
        eigen_values = diag(D);

elseif strcmpi(method,'SVD')
    y = (pinv(E)*H);
    [m, n]   = size(y);
    if m > n
        [v,s,v] = svd(y*y');
        s       = diag(s);
        v       = v(:,1);
        u       = y*v/sqrt(s(1));
        eigen_vectors = v;
    else
        [u,s,u] = svd(y'*y);
        s       = diag(s);
        u       = u(:,1);
        v       = y'*u/sqrt(s(1));
        eigen_vectors = u;
    end
    d  = sign(sum(v)); u = u*d;
    eigen_values  = u*sqrt(s(1)/n);
end    
    
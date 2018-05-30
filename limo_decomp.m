function [eigen_vectors,eigen_values] = limo_decomp(varargin)

% FORMAT: [eigen_vectors,eigen_values] = limo_decomp(E,H)
%
% INPUT E and H are matrices, typically square symmetric Sum of Squares and
%       Cross Products
%
% OUTPUT: the eigen vectors and values of the decomposition of inv(E)*H
%
% Cyril Pernet 2009
% Cyril Pernet and Iege Bassez 2017
% -----------------------------
%  Copyright (C) LIMO Team 2010

E = varargin{1};
H = varargin{2};

if nargin == 2
    cov_method = 'pseudo'; % the default
else 
    cov_method = varargin{3};
end

if strcmp(cov_method, 'pseudo')
    % procede
    [vec, D] = eig((pinv(E)*H));

    % sort eigenvalues and then sort eigenvectors in order of decreasing eigenvalues
    [e,ei] = sort(diag(D));  % eigenvalues of inv(U')*H*inv(U) == eigenvalues of inv(E)*H
    ordered_eigenvalues = flipud(e);
    vec = vec(:,flipud(ei));

    % validate if correct eigenvalues and eigenvectors of matrix pinv(E)*H:
    if round((pinv(E)*H) * vec, 4) == round(vec * diag(ordered_eigenvalues), 4)
        eigen_vectors = vec;
        eigen_values = ordered_eigenvalues;
    else
        error('this method could not find the correct eigenvalues or eigenvectors')
    end
elseif strcmp(cov_method, 'regularized')
%     n = size(Y,1);
%     [RegularizedCovariance, ~] = cov1para(Y);
%     RegularizedE = RegularizedCovariance .* (n-k);
  
    % procede
    [vec, D] = eig((inv(E)*H));

    % sort eigenvalues and then sort eigenvectors in order of decreasing eigenvalues
    [e,ei] = sort(diag(D));  % eigenvalues of inv(U')*H*inv(U) == eigenvalues of inv(E)*H
    ordered_eigenvalues = flipud(e);
    vec = vec(:,flipud(ei));

    % validate if correct eigenvalues and eigenvectors of matrix pinv(E)*H:
    if round((inv(E)*H) * vec, 4) == round(vec * diag(ordered_eigenvalues), 4)
        eigen_vectors = vec;
        eigen_values = ordered_eigenvalues;
    else
        error('this method could not find the correct eigenvalues or eigenvectors')
    end
end

end
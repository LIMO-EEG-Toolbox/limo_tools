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

% If this procedre fails (E not of full rank) we 
% use the decomposition of pinv(E)*H

if rank(E) == size(E,1)
   U = chol(E);
   [b D] = eig(inv(U')*H*inv(U)); % b: eigenvectors, D: diagonal matrix with eigenvalues
    
    % validate if correct eigenvalues and eigenvectors of matrix inv(U')*H*inv(U):
    if round((inv(U')*H*inv(U)) * b, 4) ~= round(b * D, 4) % needs to be equal 
        errordlg('something went wrong in the decomposition inv(U`)*H*inv(U)');
    else
        % adjustment to find eigenvectors of matrix inv(E)*H (page 279 Rencher) 
        a = inv(U)*b; % these are the eigenvectors of inv(E)*H
        
        % sort eigenvectors 
        [e,ei] = sort(diag(D));
        a = a(:,flipud(ei)); % in order of increasing eigenvalues
        ordered_eigenvalues = flipud(e);
        
        % validate if correct eigenvalues and eigenvectors of matrix inv(E)*H:
        if round((inv(E)*H) * a, 4) == round(a * diag(ordered_eigenvalues), 4) % needs to be one (equal)
            % check if A * v = Eigenvalue * V
            for i=1:rank(rank(inv(E)*H))
                round((inv(E)*H) * a(:,i), 4) == round(ordered_eigenvalues(i) * a(:,i), 4)
            end 
            eigen_vectors = a;
            eigen_values = ordered_eigenvalues; % deccreasing order
        end   
    end
    
else    
    A = (pinv(E)*H);
    [vec, val] = eig(A);
    [val2,ei] = sort(diag(val));
    vec = vec(:,flipud(ei)); % sorting the eigenvectors by decreasing eigenvalues
    val = flipud(val2); % decreasing eigenvalues
    left =  round((pinv(E)*H) * vec, 4);
    right = round(vec * diag(val), 4);
    % validate if correct eigenvalues and eigenvectors of matrix pinv(E)*H:
    if round(A * vec, 4) == round(vec * diag(val), 4);
        eigen_vectors = vec;
        eigen_values = val;
    end 
    %round(A * vec(:,1), 10) == round(val(1) * vec(:,1),10)
    %round(A * vec(:,2), 10) == round(val(2) * vec(:,2) ,10)
      
end
end

 
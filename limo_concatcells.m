function out = limo_concatcells(in)

% routine for recursive concatenation of a cell array of matrices
% FORMAT out = limo_concatcells(in)
%
% INPUT in is a cell arrays of matrices to concatenate in the last dimension
% OUTPUT out is the concatenate matrix
%
% Cyril Pernet Nov 2014
% -----------------------------
% Copyright (C) LIMO Team 2014

%% check data in
if isnumeric(in)
    disp('single matrix, nothing to concatenate')
    try 
        out = cell2mat(in);
    catch
        out = in; 
    end
    return
end


%% check consistency of matrices
n = length(in);
ref = size(in{1});
last_dim = numel(ref);
for data=2:n
    if sum(size(in{data}) == ref) ~= (last_dim-1)
        error('data to concatenate are of different sizes')
    end
end

%% now do it
out = cell2mat(in(1));
for concat = 2:n
    disp('concatenating matrices of data ...')
    B = cell2mat(in(concat));
    out = cat(last_dim,out,B);
end


function out = limo_concatcells(in,type)

% routine for recursive concatenation of a cell array of matrices
% FORMAT out = limo_concatcells(in,type)
%
% INPUT in is a cell array of matrices to concatenate in the last dimension
%          or a cell array of names to concatenate into a single cell
%       type is 'num' or 'char'
% OUTPUT out is a concatenate matrix or a concatenated cell of names
%
% Cyril Pernet
% ------------------------------
%  Copyright (C) LIMO Team 2021

% default is to deal with matrices
if nargin == 1
    type = 'num';
end

if strcmpi(type,'num')
    % check data in
    if isnumeric(in)
        disp('single matrix, nothing to concatenate')
        try
            out = cell2mat(in);
        catch
            out = in;
        end
        return
    end
    
    % check consistency of matrices
    n = length(in);
    ref = size(in{1});
    last_dim = numel(ref);
    for data=2:n
        if sum(size(in{data}) == ref) ~= (last_dim-1)
            error('data to concatenate are of different sizes')
        end
    end
    
    % now do it
    out = cell2mat(in(1));
    for concat = 2:n
        disp('concatenating matrices of data ...')
        B = cell2mat(in(concat));
        out = cat(last_dim,out,B);
    end
    
else % type is char
    out = in{1};
    for concat = 2:length(in)
        out = cat(2,out,in{concat});
    end
end
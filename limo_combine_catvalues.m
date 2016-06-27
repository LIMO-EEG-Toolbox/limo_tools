function out = limo_combine_catvalues(varargin)

% simple routine to combine columns from a categorical variable file
%
% FORMAT out = limo_combine_catvalue(cat_file, options)
%
% INPUT cat_file is simply the file name of your text file that include
%       multiple colimns
%       option: 'columns' followed by the columns you want to combine
%                        (default is all)
%               'format' is 'txt' to create a new text file or 'var' for
%               simply returning the new combined values (default)
%
% OUTPUT out is either the name of the new text file created or the new 
%        combined values (see optons) 
%
% EXEMPLE 
% new_values = limo_combine_catvalue('./mydocs/myfile.txt', 'column', [1 2 3])
% new_file = limo_combine_catvalue('./mydocs/myfile.txt', 'format', 'txt')
%
% Cyril Pernet 27 June 2016
% ------------------------------------------
% Copyright (C) LIMO Team 2016

%% Check data input

% defaults
columns = 'all';
format = 'val';

if nargin == 0
    
elseif nargin == 1
    
else
    
end

combine = 

%% Combine specified columns


%% output 


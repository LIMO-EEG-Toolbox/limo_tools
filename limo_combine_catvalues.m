function [single_col,final_matrix] = limo_combine_catvalues(varargin)

% simple routine to combine columns from a categorical variable file
%
% FORMAT out = limo_combine_catvalue(cat_file, options)
%
% INPUT cat_file is simply the file name of your text file that include
%                multiple columns coding each for levels in a factor
%       option: 'columns' followed by the columns you want to combine
%                        (default is all)
%               'format' is 'txt' to create a new text file or 'var' for
%               simply returning the new combined values (default)
%
% OUTPUT single_col is either the name of the new text file created or the new
%        combined values (see optons)
%        final_matrix is a matrix form of 1 and 0
%
% EXEMPLES
% new_file = limo_combine_catvalues('./mydocs/myfile.txt', 'columns', [1 2 3],'format', 'txt'))
% new_values = limo_combine_catvalues('./mydocs/myfile.txt', 'format', 'var')
%
% Cyril Pernet 27 June 2016
% ------------------------------------------
% Copyright (C) LIMO Team 2016

%% Check data input

% defaults
columns = 'all';
format = 'val';

if nargin == 0
    [filename, pathname] = uigetfile('*.txt', 'Pick a text file for your categorical variables');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel'); return
    else
        disp(['User selected ', fullfile(pathname, filename)])
        varargin{1} = fullfile(pathname, filename);
    end
    cat_values = load(filename);
elseif nargin == 1
    cat_values = load(varargin{1});
else
    cat_values = load(varargin{1});
    for n=1:(nargin-1)
        if strcmpi(varargin{n},'columns')
            columns = varargin{n+1};
        elseif strcmpi(varargin{n},'format')
            format = varargin{n+1};
        end
    end
end

if ischar(columns)
    combine = [1:size(cat_values,2)];
else
    combine = columns;
end


%% Combine specified columns

[x,nb_conditions]= limo_quick_design(cat_values(:,combine));
[tmpX interactions] = limo_make_interactions(x, nb_conditions);
N = size(tmpX,2) - interactions(end) + 1;
final_matrix = tmpX(:,N:end);

%% output
single_col = sum((repmat([1:interactions(end)],[size(final_matrix,1),1]) .* final_matrix),2);
if strcmpi(format,'txt')
    name = [varargin{1}(1:end-4) '_combined.txt'];
    save(name,'single_col','-ascii');
    single_col = name;
end

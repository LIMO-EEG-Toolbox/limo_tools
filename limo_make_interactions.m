function [tmpX interactions] = limo_make_interactions(x, nb_conditions)

% internal routine copied from limo_design_matrix to get the interaction terms
%
% FORMAT: [tmpX interactions] = limo_make_interactions(x, nb_conditions)
% 
% INPUT: x the design matrix
%        nb_conditions = vector with levels per factor e.g. [3 2 2]
% OUTPUT tmpX: the interaction matrix
%        interactions = the vector of levels per interaction e.g. [6 6 4 12]
%
% Cyril Pernet v1 - 21-06-2013
% v2 fixed bug for designs > 3 factors - 05-07-2013
% -------------------------------------------------
%  Copyright (C) LIMO Team 2015

if nb_conditions == 0
    disp('number of condition = 0')
    tmpX = []; interactions = [];
    return
end

% get each part of x for the right factors
nb_factors = size(nb_conditions,2);
F{1} = x(:,1:nb_conditions(1)); index = nb_conditions(1)+1;
for f = 2:nb_factors
    F{f} = x(:,index:(index+nb_conditions(f)-1));
    index = index+nb_conditions(f);
end

% look for which factors to combine
index = 1;
for n=2:nb_factors
    interaction{index} = nchoosek([1:nb_factors],n);
    index = index + 1; 
end

% combine those factors
tmpX = x;
index = 1;
for i=1:size(interaction,2) % for each interaction levels
    for j=1:size(interaction{i},1) % for each interaction in this level
        combination = interaction{i}(j,:);
        if size(combination,2) == 2 % 2 factors
            I = [];
            a = F{combination(1)};
            b = F{combination(2)};
            for m=1:size(a,2)
                tmp = repmat(a(:,m),1,size(b,2)).*b;
                I = [I tmp];
            end
            I = I(:,find(sum(I))); % removes the silly zero columns
            
        else % > 2 factors
            l = 1;
            while l < size(combination,2)
                if l == 1
                    I = [];
                    a = F{combination(l)};
                    b = F{combination(l+1)};
                    for m=1:size(a,2)
                        tmp = repmat(a(:,m),1,size(b,2)).*b;
                        I = [I tmp];
                    end
                    I = I(:,find(sum(I)));
                    C{l} = I; l = l+1;
                else
                    I = [];
                    a = C{l-1};
                    b = F{combination(l+1)};
                    for m=1:size(a,2)
                        tmp = repmat(a(:,m),1,size(b,2)).*b;
                        I = [I tmp];
                        clear tmp
                    end
                    I = I(:,find(sum(I)));
                    C{l} = I; l = l+1;
                end
            end
            I = C{end};
            clear C
        end
        interactions(index) = size(I,2);
        tmpX = [tmpX I];
        index = index +1;
    end
end


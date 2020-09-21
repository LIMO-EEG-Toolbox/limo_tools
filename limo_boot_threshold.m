function [mask,pvalues]= limo_boot_threshold(values,H0_values,p)
% computes the cell-wise p value using bootrap 
% 
% FORMAT: [mask,values]= limo_boot_threshold(values,H0_values,p)
% 
% INPUTs: values    = observed data (e.g. F values)
%         H0_values = distrutions of null values (e.g. bootstrap F)
%         p value   = threshold the apply for significance
%   
% OUTPUTS: mask     = significant p values
%          pvalues  = the p values from comparing values to H0_values
%
% Cyril Pernet
% ------------------------------------------------------------------
%  Copyright (C) LIMO Team 2019

tmp = NaN(size(values,1),size(values,2));
sorted_values = sort(H0_values,3); clear H0_values;
if all(sorted_values(:))>=0 % i.e. F values
    U    = round((1-p)*size(sorted_values,3));
    mask = (values >= sorted_values(:,:,U));
    for row = 1:size(values,1)
        for column = 1:size(values,2)
            tmp(row,column) = sum(values(row,column)>squeeze(sorted_values(row,column,:)));
        end
    end
    pvalues = 1- (tmp ./ size(sorted_values,3)) ; % p values
else % i.e. T values
    low  = round(p*size(sorted_values,3)/2);
    high = size(sorted_values,3) - low;
    mask = (values <= sorted_values(:,:,low))+(values >= sorted_values(:,:,high));
    for row = 1:size(values,1)
        for column = 1:size(values,2)
            tmp(row,column) = sum(values(row,column)>squeeze(sorted_values(row,column,:)));
        end
    end
    pvalues = min((tmp ./ size(sorted_values,3)), 1- (tmp ./ size(sorted_values,3))) ; % p values
end
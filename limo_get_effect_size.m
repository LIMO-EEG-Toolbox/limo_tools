function [name,clusters]=limo_get_effect_size(file,mask)

% simple routine to compute effect sizes from a result file
%
% FORMAT limo_get_effect_size(file,mask)
%
% INPUTS file: is a result file like a t-test or ANOVA
%        mask: is optional ([] by default) and is a N-ary matrix of clusters
%
% OUTPUTS name: is the name of the file created
%              [file_name]_effectsize.mat is created with Cohen's d or patial
%              eta square values at each cell
%         clusters: if mask is provided as input, the 1st eigen value if
%               provided for each cluster in mask (if some cells are outlying
%               values the 1st eigen value weights them down, otherwise it
%               the same as the mean)
%
% Cyril Pernet 2022
% ------------------------------
%  Copyright (C) LIMO Team 2022

%% reminder of standardized effect sizes
% Cohen's d = (u1-u2) / std
% Partial eta^2 np^2 = (df*F) / (df*F+dfe)
% Mahalanobid distance D = (u1-u2)*inv(S)*(u1-u2)

%% check inputs
% e.g. file = 'Rep_ANOVA_Factor_1.mat';
[filepath,filename]=fileparts(file);
if isempty(filepath)
    filepath = pwd;
end

if ~exist(fullfile(filepath,'LIMO.mat'))
    error('cannot find a LIMO.mat in the same filder as this file, this is required for this function to work')
else
    LIMO = load(fullfile(filepath,'LIMO.mat'));
    LIMO = LIMO.LIMO;
end

%% compute effect sizes based on design

if contains(FileName,'one_sample','IgnoreCase',true) || contains(FileName,'two_samples','IgnoreCase',true) || ...
        contains(FileName,'paired_samples','IgnoreCase',true) || contains(FileName,'con_','IgnoreCase',true) || ...
        contains(FileName,'ess_','IgnoreCase',true)
    
elseif contains(LIMO.design.name,'regression','IgnoreCase',true) && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
        contains(LIMO.design.name,'ANOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
        contains(LIMO.design.name,'ANCOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
    
elseif contains(LIMO.design.name,'Repeated','IgnoreCase',true)   % All stuffs for repeated measures ANOVA
    F = load(filename);
    F = F.(cell2mat(fieldnames(F)));
    if numel(size(F)) == 3
        F = squeeze(F(:,:,1));
    else
        F = squeeze(F(:,:,:,1));
    end
    
    if ~contains(FileName,'Rep_ANOVA_Interaction') && ~contains(FileName,'Rep_ANOVA_Gp')
        T2 = (df/dfe)*F;
        D  = (N*T2) / n1*n2;
    elseif contains(FileName,'Rep_ANOVA_Gp')
        
    elseif contains(FileName,'Rep_ANOVA_Interaction')
        
    end
    
end



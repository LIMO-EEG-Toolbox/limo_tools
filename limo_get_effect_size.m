function [name,clusters] = limo_get_effect_size(file,mask)

% simple routine to compute effect sizes from a result file
%
% FORMAT [name,clusters] = limo_get_effect_size(file,mask)
%
% INPUTS file: is a result file like a t-test or ANOVA
%        mask: is optional ([] by default) and is a N-ary matrix of clusters
%
% OUTPUTS name: is the name of the file created
%              [file_name]_effectsize.mat is created with Cohen's d or patial
%              eta square values at each cell
%         clusters: if mask is provided as input, it returns a structure with
%               summary statistics per cluster for the computed effect size
%               (not for the T/F statistics, use limo_get_summary.m)
%               - eigenmode is the 'direction' of the effect size (usually,
%                        clusters are right skewed and thus it represents
%                        better the 'average' effect size)
%               - median provided as a comparison point to eigen mode
%               - mean provided as a comparison point to eigen mode
%               - min and max for completeness
%
% Cyril Pernet 2022
% ------------------------------
%  Copyright (C) LIMO Team 2022

name     = [];
clusters = [];

%% reminder of standardized effect sizes
% Cohen's d = (u1-u2) / std
% Partial eta^2 np^2 = (df*F) / (df*F+dfe)
% Mahalanobis distance D = (u1-u2)*inv(S)*(u1-u2) = % Hotelling Tsquare / N

%% check inputs

if nargin == 0
    file = uigetfile('.mat','select a LIMO stat file');
    if isempty(file)
        return
    end
end

[filepath,filename,ext]=fileparts(file);
if isempty(filepath)
    filepath = pwd;
end
filename = [filename ext];
assert(exist(fullfile(filepath,filename),'file'), 'file %s not found', filename)

if ~exist(fullfile(filepath,'LIMO.mat'),'file')
    error('cannot find a LIMO.mat in the same filder as this file, this is required for this function to work')
else
    LIMO = load(fullfile(filepath,'LIMO.mat'));
    LIMO = LIMO.LIMO;
end

%% compute effect sizes based on design

if contains(filename,'one_sample','IgnoreCase',true) || contains(filename,'two_samples','IgnoreCase',true) || ...
        contains(filename,'paired_samples','IgnoreCase',true)
    T = load(fullfile(filepath,filename));
    T = T.(cell2mat(fieldnames(T)));
    if numel(size(T)) == 3
        mu  = squeeze(T(:,:,1));        
        se  = squeeze(T(:,:,2));
        df  = squeeze(T(:,:,3));
    else
        mu  = squeeze(T(:,:,:,1));
        se  = squeeze(T(:,:,:,2));
        df  = squeeze(T(:,:,:,3));
    end
    n           = df+1;
    effect_size = mu ./ (se.*sqrt(n));
    name        = fullfile(filepath,[filename(1:end-4) '_Cohensd.mat']);
    
elseif contains(LIMO.design.name,'regression','IgnoreCase',true) && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
        contains(LIMO.design.name,'ANOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true) || ...
        contains(LIMO.design.name,'ANCOVA') && ~contains(LIMO.design.name,'Repeated','IgnoreCase',true)
    stats = load(fullfile(filepath,filename));
    stats = stats.(cell2mat(fieldnames(stats)));
    if numel(size(stats)) == 3
        if size(stats,3) == 2
            stats = squeeze(stats(:,:,1));
        else
            if contains(filename,'R2')
                R2    = squeeze(stats(:,:,1));
            else
                mu    = squeeze(stats(:,:,1));
                se    = squeeze(stats(:,:,2));
                df    = squeeze(stats(:,:,3));
                stats = squeeze(stats(:,:,4));
            end
        end
    else
        if size(stats,4) == 2
            stats = squeeze(stats(:,:,:,1));
        else
            if contains(filename,'R2')
                R2    = squeeze(stats(:,:,:,1));
            else
                mu    = squeeze(stats(:,:,:,1));
                se    = squeeze(stats(:,:,:,2));
                df    = squeeze(stats(:,:,:,3));
                stats = squeeze(stats(:,:,:,4));
            end
        end
    end
    
    if contains(filename,'con')
        n           = df+1;
        effect_size = mu ./ (se.*sqrt(n));
        name        = fullfile(filepath,[filename(1:end-4) '_Cohensd.mat']);
    else
        if contains(LIMO.design.name,'regression','IgnoreCase',true) && ...
                contains(filename,'Covariate')
            A   = LIMO.model.continuous_df(1)*stats;
            B   = (A+repmat(LIMO.model.continuous_df(2),size(A,1),size(A,2)));
            effect_size = A ./B ;
            name        = fullfile(filepath,[filename(1:end-4) '_PartialEta2.mat']);
        elseif contains(LIMO.design.name,'regression','IgnoreCase',true) && ...
                contains(filename,'R2')
            effect_size = R2 ./ (1-R2);
            name        = fullfile(filepath,[filename(1:end-4) '_Cohensf2.mat']);
        else
            if contains(LIMO.design.name,'ANOVA')
                A       = stats.*repmat(LIMO.design.df,1,size(stats,2));
                B       = A+LIMO.design.dfe;
            elseif contains(LIMO.design.name,'ANCOVA')
                if contains(filename,'condition','IgnoreCase',true)
                    A       = LIMO.model.conditions_df(1).*stats;
                    B       = A+repmat(LIMO.model.conditions_df(2),size(A,1),size(A,2));
                elseif contains(filename,'covariate','IgnoreCase',true)
                    A       = LIMO.model.continuous_df(1).*stats;
                    B       = A+repmat(LIMO.model.continuous_df(2),size(A,1),size(A,2));
                end
            end
            effect_size = A ./B ;
            name        = fullfile(filepath,[filename(1:end-4) '_PartialEta2.mat']);
        end
    end
    
elseif contains(LIMO.design.name,'Repeated','IgnoreCase',true)   % All stuffs for repeated measures ANOVA
    F = load(fullfile(filepath,filename));
    F = F.(cell2mat(fieldnames(F)));
    if numel(size(F)) == 3
        if size(F,3) == 2
            F   = squeeze(F(:,:,1));
        else
            df  = squeeze(F(:,:,3));
            F   = squeeze(F(:,:,4));
        end
    else
        if size(F,4) == 2
            F   = squeeze(F(:,:,:,1));
        else
            df  = squeeze(F(:,:,:,3));
            F   = squeeze(F(:,:,:,4));
        end
    end
    
    if ~contains(filename,'Rep_ANOVA_Interaction') && ~contains(filename,'Rep_ANOVA_Gp')
        if contains(filename,'Main_effect','IgnoreCase',true)
            index1     = strfind(filename,'Main_effect')+length('Main_effect')+1;
            index2     = max(strfind(filename,'_'))-1;
            effect_nb  = eval(filename(index1:index2));
        elseif contains(filename,'Interaction','IgnoreCase',true)
            index1     = strfind(filename,'Interaction')+length('Interaction')+1;
            index2     = max(strfind(filename,'_'))-1;
            effect_nb  = eval(filename(index1:index2));
        else
            index1     = strfind(filename,'Factor')+length('Factor')+1;
            effect_nb  = eval(filename(index1:end));
        end
        df          = squeeze(LIMO.design.df(:,effect_nb));
        dfe         = squeeze(LIMO.design.dfe(:,effect_nb));
        T2          = F.*repmat((df./dfe),1,size(F,2));
        N           = size(LIMO.design.X,1)/size(LIMO.design.C{effect_nb},2);
        effect_size = T2 ./ N;
        name        = fullfile(filepath,[filename(1:end-4) '_MahalanobisD.mat']);
    elseif contains(filename,'Rep_ANOVA_Gp')
        A           = (LIMO.design.group.df'.*F);
        B           = (A+repmat(LIMO.design.group.dfe',1,size(A,2)));
        effect_size = A ./B ;
        name        = fullfile(filepath,[filename(1:end-4) '_PartialEta2.mat']);
    elseif contains(filename,'Rep_ANOVA_Interaction')
        effect_nb   = str2double(filename(strfind(filename,'Factor_')+7:strfind(filename,'Factor_')+6+strfind(filename(strfind(filename,'Factor_')+6:end),'_')));
        df          = squeeze(LIMO.design.interaction.df(:,effect_nb));
        dfe         = squeeze(LIMO.design.interaction.dfe(:,effect_nb));
        T2          = F.*repmat((df./dfe),1,size(F,2));
        N           = size(LIMO.design.X,1)/size(LIMO.design.C{effect_nb},2);
        effect_size = T2 ./ N;
        name        = fullfile(filepath,[filename(1:end-4) '_MahalanobisD.mat']);
    elseif contains(filename,'ess')
        if contains(filename,'gp_interaction')
            effect_nb   = str2double(filename(strfind(filename,'ess_gp_interaction_')+19:end-4));
        else
            effect_nb   = str2double(filename(strfind(filename,'ess_')+4:end-4));
        end
        N           = size(LIMO.design.X,1)/size(LIMO.contrast{effect_nb}.C,2);
        T2          = F.*(df./(N-df));
        effect_size = T2 ./ N;
    end
end

if exist('name','var')
    save(name,'effect_size');
else
   error('effect size not computed, likely filename not handled') 
end

%% deal with clusters
% --------------------
if exist('mask','var')
    % quickly make if N-ary if binary
    if length(unique(mask)) == 2
        mask = limo_findcluster(mask,LIMO.data.neighbouring_matrix,2);
    end
    
    num = unique(mask);
    num(num==0) = [];
    for c = size(num,2):-1:1
        data                  = effect_size(mask == num(c));
        clusters(c).eigenmode = sqrt(eig(data'*data)/length(data));
        clusters(c).median    = median(data);
        clusters(c).mean      = mean(data);
        clusters(c).min       = min(data);
        clusters(c).max       = max(data);
    end
end


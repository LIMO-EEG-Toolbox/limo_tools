function tfce_score = limo_tfce_select(filename,LIMO)

% data handling function for limo_tfce, ie fiegure out what file it
% is, and what dimension to compute tfce on -
%
% FORMAT limo_tfce_sect(filename,LIMO)
%
% INPUT filename is a file genrate by LIMO EEG (eg one_sample_ttest_parameters_X.mat)
%       LIMO is the LIMO structure
%
% OUTPUT tfce_score is the score matrix computed in limo_tfce
%
% Cyril Pernet v1 08-02-2016
% --------------------------------------
% Copyright (C) LIMO Team 2016

%% check data structure
if isfield(LIMO.data,'neighbouring_matrix') && LIMO.design.bootstrap ~=0
    if exist('TFCE','dir')
        if strcmp(questdlg('TFCE directory detected, overwrite?','data check','Yes','No','No'),'No');
            return
        end
    end
    fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Computing TFCE GLM takes a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
    mkdir TFCE; PCT_test = ver('distcomp');
elseif ~isfield(LIMO.data,'neighbouring_matrix')
    disp('No TFCE performed, neighbourhood matrix missing')
elseif  LIMO.design.bootstrap ==0
    disp('No TFCE performed, since there was no bootstraps computed')
end

%% GLM cases
if LIMO.level == 1; % || LIMO.level == 2 && LIMO.design.name
    cd(LIMO.data.dir)
    
    % R2
    load R2.mat; fprintf('Creating R2 TFCE scores \n'); cd('TFCE');
    if size(R2,1) == 1
        tfce_score(1,:) = limo_tfce(1, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
    else
        tfce_score = limo_tfce(2, squeeze(R2(:,:,2)),LIMO.data.neighbouring_matrix);
    end
    save('tfce_R2','tfce_score'); clear R2; cd ..;
    
    cd('H0'); fprintf('Thresholding H0_R2 using TFCE \n'); load H0_R2;
    if size(H0_R2,1) == 1
        if ~isempty(PCT_test)
            tfce_H0_score = NaN(1,size(H0_R2,2),LIMO.design.bootstrap);
            parfor b=1:nboot
                tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
            end
        else
            tfce_H0_score(1,:,:) = limo_tfce(1, squeeze(H0_R2(:,2,:)),LIMO.data.neighbouring_matrix);
        end
    else
        if ~isempty(PCT_test)
            tfce_H0_score = NaN(size(H0_R2,1),size(H0_R2,2),LIMO.design.bootstrap);
            parfor b=1:nboot
                tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_R2(:,:,2,b)),LIMO.data.neighbouring_matrix,0);
            end
        else
            tfce_H0_score = limo_tfce(2, squeeze(H0_R2(:,:,2,:)),LIMO.data.neighbouring_matrix);
        end
    end
    save('tfce_H0_R2','tfce_H0_score'); clear H0_R2; cd ..;
    
    % conditions
    if prod(LIMO.design.nb_conditions) ~=0
        for i=1:length(LIMO.design.nb_conditions)
            name = sprintf('Condition_effect_%g.mat',i); load(name);
            cd('TFCE'); fprintf('Creating Condition %g TFCE scores \n',i)
            if size(Condition_effect,1) == 1
                tfce_score(1,:) = limo_tfce(1, squeeze(Condition_effect(:,:,1)),LIMO.data.neighbouring_matrix);
                tfce_score = limo_tfce(2, squeeze(Condition_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            end
            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
            clear Condition_effect tfce_score; cd ..
        end
        
        cd('H0'); fprintf('Creating H0 Condition(s) TFCE scores \n');
        for i=1:length(LIMO.design.nb_conditions)
            name = sprintf('H0_Condition_effect_%g.mat',i); load(name);
            if size(H0_Condition_effect,1) == 1
                if exist('parfor','file')
                    tfce_H0_score = NaN(1,size(H0_Condition_effect,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Condition_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score(1,:,:) = limo_tfce(1,squeeze(H0_Condition_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                end
            else
                if exist('parfor','file')
                    tfce_H0_score = NaN(size(H0_Condition_effect,1),size(H0_R2,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Condition_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score = limo_tfce(2,squeeze(H0_Condition_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                end
            end
            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_H0_score');
            clear H0_Condition_effect tfce_H0_score;
        end
        cd ..
    end
    
    % interactions
    if LIMO.design.fullfactorial == 1
        for i=1:length(LIMO.design.fullfactorial)
            name = sprintf('Interaction_effect_%g.mat',i); load(name);
            cd('TFCE'); fprintf('Creating Interaction %g TFCE scores \n',i)
            if size(Interaction_effect,1) == 1
                tfce_score(1,:) = limo_tfce(1,squeeze(Interaction_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            else
                tfce_score = limo_tfce(2,squeeze(Interaction_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            end
            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
            clear Interaction_effect tfce_score; cd ..
        end
        
        cd('H0'); fprintf('Creating H0 Interaction(s) TFCE scores \n');
        for i=1:length(LIMO.design.fullfactorial)
            name = sprintf('H0_Interaction_effect_%g.mat',i); load(name);
            if size(H0_Interaction_effect,1) == 1
                if exist('parfor','file')
                    tfce_H0_score = NaN(1,size(H0_Interaction_effect,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Interaction_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score(1,:,:) = limo_tfce(1,squeeze(H0_Interaction_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                end
            else
                if exist('parfor','file')
                    tfce_H0_score = NaN(size(H0_Interaction_effect,1),size(H0_Interaction_effect,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score = limo_tfce(2,squeeze(H0_Interaction_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                end
            end
            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_H0_score');
            clear H0_Interaction_effect tfce_H0_score;
        end
        cd ..
    end
    
    % covariates / continuous regressors
    if LIMO.design.nb_continuous ~=0
        for i=1:LIMO.design.nb_continuous
            name = sprintf('Covariate_effect_%g.mat',i); load(name);
            cd('TFCE'); fprintf('Creating Covariate %g TFCE scores \n',i);
            if size(Covariate_effect,1) == 1
                tfce_score(1,:) = limo_tfce(1,squeeze(Covariate_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            else
                tfce_score = limo_tfce(2,squeeze(Covariate_effect(:,:,1)),LIMO.data.neighbouring_matrix);
            end
            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_score');
            clear Covariate_effect tfce_score; cd ..
        end
        
        cd('H0'); fprintf('Creating H0 Covariate(s) TFCE scores \n');
        for i=1:LIMO.design.nb_continuous
            name = sprintf('H0_Covariate_effect_%g.mat',i); load(name);
            if size(H0_Covariate_effect,1) == 1
                if ~isempty(PCT_test)
                    tfce_H0_score = NaN(1,size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(1,:,b) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score(1,:,:) = limo_tfce(1,squeeze(H0_Covariate_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                end
            else
                if ~isempty(PCT_test)
                    tfce_H0_score = NaN(size(H0_Covariate_effect,1),size(H0_Covariate_effect,2),LIMO.design.bootstrap);
                    parfor b=1:nboot
                        tfce_H0_score(:,:,b) = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,b)),LIMO.data.neighbouring_matrix,0);
                    end
                else
                    tfce_H0_score = limo_tfce(2,squeeze(H0_Covariate_effect(:,:,1,:)),LIMO.data.neighbouring_matrix);
                end
            end
            full_name = sprintf('tfce_%s',name); save(full_name,'tfce_H0_score');
            clear H0_Covariate_effect tfce_H0_score
        end
        cd ..
    end
end

%% 2nd level T-tests and 1way ANOVA

%% 2nd level rep. measure ANOVA



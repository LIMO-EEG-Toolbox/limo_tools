function limo_contrast_execute(LIMOfile, handles)

% Execute contrast analysis
%
% Nicolas Chauveau, Arnaud Delorme & Cyril Pernet 29-04-2009 v2
% updated for 1st level bootstrap + some fix 20-06-2013
% ------------------------------
%  Copyright (C) LIMO Team 2019

if ischar(LIMOfile)
    LIMO = load('-mat', LIMOfile);
    LIMO = LIMO.LIMO;
else
    LIMO = LIMOfile;
end

disp('executing contrast')

if LIMO.design.bootstrap ~=0 && exist([LIMO.dir filesep 'H0'],'dir') && ...
        ~contains(LIMO.design.method,'Generalized Welch''s method','IgnoreCase',true)
    choice = questdlg('(re)compute contrast bootstrap?','bootstrap choice','compute bootstrap contrast','don''t compute any bootstraps','compute bootstrap contrast');
else
    choice = 'don''t compute any bootstraps';
end

% ------------------------------------------------
% 1st level contrat & 2st level ANOVA/ANCOVA/Regression
% ------------------------------------------------
if LIMO.Level == 1 || ...
        LIMO.Level == 2 && contains(LIMO.design.name,'regression','IgnoreCase',true) || ...
        LIMO.Level == 2 && contains(LIMO.design.name,'N-ways','IgnoreCase',true) || ...
        LIMO.Level == 2 && contains(LIMO.design.name,'ANCOVA','IgnoreCase',true) || ...
        LIMO.Level == 2 && contains(LIMO.design.name,'ANOVA','IgnoreCase',true) && ...
        ~contains(LIMO.design.name,'Repeated')

    if isfield(LIMO,'contrast')
        previous_con = size(LIMO.contrast,2);
    else
        previous_con = 0;
    end
    index = previous_con+1;

    % update LIMO.mat
    LIMO.contrast{index}.C = handles.C;
    if handles.F == 0
        LIMO.contrast{index}.V = 'T';
    else
        LIMO.contrast{index}.V = 'F';
    end

    if exist(LIMO.dir,'dir')
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
    end

    % -------------------------------------------------------
    if strcmp(LIMO.design.type_of_analysis,'Mass-univariate')
        % -------------------------------------------------------

        if contains(LIMO.design.name,'ANOVA','IgnoreCase',true) && ...
                contains(LIMO.design.method,'Generalized Welch''s method','IgnoreCase',true)
            if handles.F == 1
                warndlg(sprintf('there is no F contrast possible for Generalized Welch''s method ANOVA'),'Robust ANOVA info')
            else
                warndlg(sprintf('no T contrasts for Generalized Welch''s method ANOVA,\nswitching to robust t-tests for sub-groups comparison'),'Robust ANOVA info')
                if exist(LIMO.dir,'dir')
                    data = load(fullfile(LIMO.dir,'Yr.mat'));
                else
                    warning('%s doesn''t exist, pulling data from the local dir',LIMO.dir)
                    data = load(fullfile(pwd,'Yr.mat')); LIMO.dir = pwd;
                end

                if strcmp(LIMO.Analysis ,'Time-Frequency')
                    limo_random_robust(2,data.Yr(:,:,:,find(LIMO.design.X(:,handles.C == 1))),...
                        data.Yr(:,:,:,find(LIMO.design.X(:,handles.C == -1))),...
                        find(handles.C ~= 0),LIMO);
                else
                    limo_random_robust(2,data.Yr(:,:,find(LIMO.design.X(:,handles.C == 1))),...
                        data.Yr(:,:,find(LIMO.design.X(:,handles.C == -1))),...
                        find(handles.C ~= 0),LIMO);
                end
            end
        else % standard GLM type ANOVA/ANCOVA

            if ~exist(LIMO.dir,'dir')
                LIMO.dir = pwd;
            end
            Yfile = dir(fullfile(LIMO.dir,'*Yr.mat'));
            Bfile = dir(fullfile(LIMO.dir,'*Betas.mat'));
            limo_contrast(fullfile(LIMO.dir,Yfile.name), fullfile(LIMO.dir,Bfile.name), LIMO, handles.F,1);

            if LIMO.design.bootstrap ~= 0 && strcmpi(choice,'compute bootstrap contrast')
                Yr       = load(fullfile(LIMO.dir,Yfile.name)); 
                Yr       = Yr.Yr;
                H0Bfile  = dir(fullfile(LIMO.dir,['H0' filesep '*BetasH0.mat']));
                H0_Betas = load(fullfile(H0Bfile.folder,H0Bfile.name)); 
                H0_Betas = H0_Betas.H0_Betas;
                if strcmp(LIMO.Analysis ,'Time-Frequency')
                    disp('preparing Time-Frequency H0 data matrix');
                    tmp = zeros(size(H0_Betas,1), size(H0_Betas,2)*size(H0_Betas,3), size(H0_Betas,4), size(H0_Betas,5));
                    for boot = 1:size(H0_Betas,5)
                        tmp(:,:,:,boot)= limo_tf_4d_reshape(squeeze(H0_Betas(:,:,:,:,boot)));
                    end
                    limo_contrast(limo_tf_4d_reshape(Yr), tmp, LIMO, handles.F,2);
                    clear Yr tmp
                else
                    limo_contrast(Yr, H0_Betas, LIMO, handles.F,2);
                end
                clear Yr tmp
                disp('boostrapped contrasts done ...')
            end
        end

        % -------------------------------------------------------
    elseif strcmp(LIMO.design.type_of_analysis,'Multivariate')
        % -------------------------------------------------------

        LIMO.contrast = handles.F;
        limo_contrast(squeeze(Yr(:,time,:))', squeeze(Betas(:,time,:))', [], LIMO, handles.F,1);
        if exist(LIMO.dir,'dir')
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3');
        end

    end
    clear Yr Betas

    % -------------------------------------------
    %          2nd level Repeated measure ANOVA
    % -------------------------------------------
elseif LIMO.Level == 2 && contains(LIMO.design.name,'Repeated')

    if isfield(LIMO,'contrast')
        previous_con = size(LIMO.contrast,2);
    else
        previous_con = 0;
    end
    index = previous_con+1;

    % update LIMO.mat
    LIMO.contrast{index}.C = handles.C;
    LIMO.contrast{index}.V = 'F'; % always F since we use Hotelling test

    % create ess files and call limo_rep_anova adding C
    Yr = load(fullfile(LIMO.dir,'Yr.mat')); Yr = Yr.Yr;
    if exist(LIMO.dir,'dir')
        save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO','-v7.3')
    end
    limo_contrast(Yr,LIMO,3);

    if strcmpi(choice,'compute bootstrap contrast')
        limo_contrast(Yr, LIMO, 4);
    end

    if LIMO.design.tfce == 1 && strcmpi(choice,'compute bootstrap contrast')
        filename = fullfile(LIMO.dir,['ess_' num2str(index) '.mat']);
        limo_tfce_handling(filename)
        if LIMO.design.nb_conditions ~= 1
            filename = fullfile(LIMO.dir,['ess_gp_interaction_' num2str(index) '.mat']);
            limo_tfce_handling(filename)
        end
    end
    clear Yr LIMO
    disp('contrast evaluation done ...')
end
        

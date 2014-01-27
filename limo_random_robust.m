function limo_random_robust(varargin)

% This function makes the result files for the random effects of various tests
% as well as organizes and makes files for boostrap. Tt is interfaced with
% limo_random_effect which itself interfaces with the user to select and pass
% the data in the appropriate format. Limo_random_robust calls low level
% functions like limo_ttest to perform the actual computation.
%
% FORMAT
%
% limo_random_robust(1,y,parameter number,nboot,tfce)
%                    1                = a one-sample t-test
%                    y                = data (dim electrodes, frames, subjects)
%                    parameter number = describe which parameters is analysed (e.g. 1)
%                    nboot            = nb of resamples
%                    tfce             = 1 run tcfe.
%                   
% limo_random_robust(2,y1,y2,parameter number,nboot,tfce);
%                    2                = two samples t-test
%                    y1               = data (dim electrodes, frames, subjects)
%                    y2               = data (dim electrodes, frames, subjects)
%                    parameter number = describe which parameters is analysed (e.g. 1)
%                    nboot            = nb of resamples
%                    tfce             = 1 run tcfe.
%
% limo_random_robust(3,y1,y2,parameter number,nboot,tfce);
%                    3                = paired t-test
%                    y1               = data (dim electrodes, frames, subjects)
%                    y2               = data (dim electrodes, frames, subjects)
%                    parameter number = describe which parameters were analysed (e.g. 1 2)
%                    nboot            = nb of resamples
%                    tfce             = 1 run tcfe.
%
% limo_random_robust(4,y,X,parameter number,nboot,tfce);
%                    4                = regression analysis
%                    y                = data (dim electrodes, frames, subjects)
%                    X                = continuous regressor(s)
%                    parameter number = describe which parameters is analysed (e.g. 1)
%                    nboot            = nb of resamples
%                    tfce             = 1 run tcfe.
%
% limo_random_robust(5,y,cat,cont,nboot,tfce)
%                    5     = N-way ANOVA/ANCOVA
%                    y     = data (dim electrodes, frames, subjects)
%                    cat   = categorical variable(s)
%                    cont  = continuous regressors (covariates) 
%                    nboot = nb of resamples
%                    tfce  = 1 run tcfe.
%
% limo_random_robust(6,y,gp,factor_levels,nboot,tfce)
%                    6 = Repeated measures ANOVA/ANCOVA using multivariate approach
%                    y = data (dim electrodes, frames, subjects, measures)
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    nboot = nb of resamples
%                    tfce = 1 run tcfe.
%
% OUPUT
%
% write on the disk matrices correponding to the test
% 1 sample t-test: one_parameter_X (electrodes, frames, [mean value, std, nb_subjects, t, p])
%                  boot_one_sample_X (electrodes, frames, [T values under H1, T values under H0, p values under H0], nboot)
%
% See also LIMO_TRIMCI LIMO_YUEN_TTEST LIMO_YUEND_TTEST LIMO TFCE
% LIMO_DESIGN_MATRIX LIMO_GLM1 LIMO_EEG(4) LIMO_REP_ANOVA
% LIMO_CREATE_BOOT_TABLE
% -----------------------------
%  Copyright (C) LIMO Team 2010

% v1: Cyril Pernet and Guillaume Rousselet 24-08-2009
% v2: Cyril Pernet 12-07-2010
% 25-08-2010: GAR fixed bug in one-sample bootstrap + added new NaN check of boot indices
% 29-08-2010:  Cyril Pernet and Guillaume Rousselet ANOVAs sorted out (boot_index per cell)
% 11/12-2012 Marianne Latinus added TFCE computation. Also added a check
% that if a bootstrap file with same characteristics (nboot, and nelec) on
% the same data exists the bootstrap step is skipped.
% 02/06/2013 cleaned up + changed ANOVA to create matrices amd run
% limo_glm(4) - updated repeated measure - boot_tables are made via
% function -- thx to Benedikt Ehinger for spotting a few bugs here and
% there
% Novembre 2013 - Matt Craddock fixes applied on rep ANOVA for single
% electrode analyses


%% inputs checks for tfce
type  = varargin{1};
alpha = .05; % used for a basic computation like in t-test but data aren't thresholded here

if ismember(type, [1])
    tfce = varargin{5};
else
    tfce = varargin{6};
end

if tfce == 1
    files = dir('LIMO.mat');
    if isempty(files)
        [file,dir_path] = uigetfile('LIMO.mat','select a LIMO.mat file');
        if file ==0
            return
        else
            cd(dir_path); load LIMO.mat;
        end
    else
        load('LIMO.mat');
    end
    cd(LIMO.dir);
    if isfield(LIMO.data,'neighbouring_matrix')
        channeighbstructmat = LIMO.data.neighbouring_matrix;
    else
        disp('tfce found, creating neighbouring matrix')
        answer = questdlg('load or compute neighbouring matrix?','channel neighbouring definition','Load','Compute','Load');
        if strcmp(answer,'Load')
            [file,path,whatsup] = uigetfile('*.mat','select neighbourghing matrix (or expected chanloc file)');
            if whatsup == 0
                disp('selection aborted');
                return
            else
                cd(path); load(file); cd(LIMO.dir);
            end
        else
            channeighbstructmat = limo_expected_chanlocs(LIMO.data.data, LIMO.data.data_dir);
        end
        LIMO.data.neighbouring_matrix = channeighbstructmat;
        save LIMO LIMO
    end
end
        
%% start

switch type
    %--------------------------------------------------------------------------
    %            One Sample t-test     //     bootstrap-t method
    %--------------------------------------------------------------------------
    case {1}
        
        data                = varargin{2}; 
        parameter           = varargin{3};
        nboot               = varargin{4};
        tfce                = varargin{5};
        clear varargin
        
        % ------------------------------------------------
        % make a one_sample file per parameter (electrodes, frames, [mean value, se, df, t, p])
        % se and df can be used to build CI and p to threshhold the data
        one_sample = NaN(size(data,1), size(data,2), 5);
        name = sprintf('one_sample_ttest_parameter_%g',parameter);
       
        for electrode = 1:size(data,1) % run per electrode because we have to remove NaNs
            fprintf('analyse parameter %g electrode %g \n',parameter, electrode); 
            tmp = data(electrode,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
            [one_sample(electrode,:,4),one_sample(electrode,:,1),trimci,one_sample(electrode,:,2),one_sample(electrode,:,5),tcrit,one_sample(electrode,:,3)]=limo_trimci(Y);
            % [one_sample(electrode,:,1),one_sample(electrode,:,3),ci,sd,n,one_sample(electrode,:,4),one_sample(electrode,:,5)] = limo_ttest(1,Y,0,5/100);
            % one_sample(electrode,:,2) = sd./sqrt(n);
            clear tmp Y
        end
        save ([name],'one_sample', '-v7.3')

        
        % ------------------------------------------------
        % Bootstrap
        if nboot > 0
            
            bootex = 1;
            boot_name = sprintf('H0_one_sample_ttest_parameter_%g',parameter);
            if exist(['H0', filesep, boot_name, '.mat'], 'file')
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'Yes');
                    bootex = 1;
                else 
                    bootex = 0;
                end
            end
            
            if bootex == 1;
                mkdir H0
                % create a boot one_sample file to store data under H0 and H1
                H0_one_sample = NaN(size(data,1), size(data,2),2,nboot); % stores T and p values for each boot under H0
                % create centered data to estimate H0
                centered_data = data - repmat(nanmean(data,3),[1 1 size(data,3)]);
                % get boot table
                disp('making boot table ...')
                boot_table = limo_create_boot_table(data,nboot);
                save(['H0', filesep, 'boot_table'], 'boot_table')

                % get results under H0
                for electrode = 1:size(data,1)
                    fprintf('bootstrap: electrode %g parameter %g \n',electrode,parameter);
                    tmp = centered_data(electrode,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                    for b=1:nboot
                        [H0_one_sample(electrode,:,1,b),tmdata,trimci,se,H0_one_sample(electrode,:,2,b),tcrit,df]=limo_trimci(Y(1,:,boot_table{electrode}(:,b)));
                        % [m,dfe,ci,sd,n,H0_one_sample(electrode,:,1,b),H0_one_sample(electrode,:,2,b)] = limo_ttest(1,Y(1,:,boot_table{electrode}(:,b)),0,5/100);
                    end
                    clear tmp Y
                end % closes for electrode
                save (['H0', filesep, boot_name],'H0_one_sample','-v7.3');
            end
            
            if tfce ~= 0 && size(data,1) > 1
                mkdir tfce
                tfce_name = sprintf('tfce_one_sample_ttest_parameter_%g',parameter);
                tfce_H0_name = sprintf('tfce_H0_one_sample_ttest_parameter_%g',parameter);
                % do tfce for the current data
                fprintf('Thresholding One Sample T-test using TFCE \n');
                tfce_one_sample = limo_tfce(squeeze(one_sample(:,:,4)),channeighbstructmat);
                save(['tfce', filesep, tfce_name], 'tfce_one_sample', '-v7.3'); clear tfce_one_sample;
                fprintf('Thresholding H0 One Sample T-test using TFCE \n');
                tfce_H0_one_sample = limo_tfce(squeeze(H0_one_sample(:,:,1,:)),channeighbstructmat);
                save(['H0', filesep, tfce_H0_name],'tfce_H0_one_sample', '-v7.3'); clear tfce_H0_one_sample;
            end
            clear H0_one_sample
        end % closes if nboot > 0
        clear one_sample
        disp('one sample t-test done')
        
        
        %--------------------------------------------------------------------------
        %            Two Samples t-test     //     percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {2}
        
        data1               = varargin{2}; 
        data2               = varargin{3}; 
        parameter           = varargin{4};
        nboot               = varargin{5};
        tfce                = varargin{6};
        clear varargin
        
        % ------------------------------------------------ 
        % make a two_samples file per parameter (electrodes, frames, [mean value, se, df, t, p])
        two_samples = NaN(size(data1,1), size(data1,2),5);
        name = sprintf('two_samples_ttest_parameter_%g',parameter);
        
        array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
        for e = 1:size(array,1)
            electrode = array(e);
            fprintf('analyse parameter %g electrode %g',parameter, electrode); disp(' ');
            tmp = data1(electrode,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            tmp = data2(electrode,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            [two_samples(electrode,:,4),two_samples(electrode,:,1),two_samples(electrode,:,2),CI,two_samples(electrode,:,5),tcrit,two_samples(electrode,:,3)]=limo_yuen_ttest(Y1,Y2); clear Y1 Y2
            % [two_samples(electrode,:,1),two_samples(electrode,:,3),ci,sd,n,two_samples(electrode,:,4),two_samples(electrode,:,5)]=limo_ttest(2,Y1,Y2,.05); 
            % sd = sd.^2; a = sd(1,:)./size(Y1,3); b = sd(1,:)./size(Y2,3);
            % two_samples(electrode,:,2) = sqrt(a + b); clear Y1 Y2
        end
        save ([name],'two_samples', '-v7.3')
        
        % ------------------------------------------------
        if nboot > 0

            bootex = 1;
            boot_name = sprintf('H0_two_samples_ttest_parameter_%g',parameter);
            if exist(['H0', filesep, boot_name, '.mat'], 'file')
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'Yes');
                    bootex = 1;
                else 
                    bootex = 0;
                end
            end
            
            if bootex == 1;
                mkdir H0
                % create a boot one_sample file to store data under H0 
                H0_two_samples = NaN(size(data1,1), size(data1,2), 2, nboot); % stores T and p values for each boot
                % create centered data to estimate H0
                data1_centered = data1 - repmat(nanmean(data1,3),[1 1 size(data1,3)]);
                data2_centered = data2 - repmat(nanmean(data2,3),[1 1 size(data2,3)]);
                % get boot table
                disp('making boot tables ...')
                boot_table1 = limo_create_boot_table(data1,nboot);
                boot_table2 = limo_create_boot_table(data2,nboot);
                save(['H0', filesep, 'boot_table1'], 'boot_table1')
                save(['H0', filesep, 'boot_table2'], 'boot_table2')

                % get results under H0
                for e = 1:size(array,1)
                    electrode = array(e);
                    fprintf('bootstrapping electrode %g/%g parameter %g \n',e,size(array,1),parameter);
                    tmp = data1_centered(electrode,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    tmp = data2_centered(electrode,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    for b=1:nboot
                        [H0_two_samples(electrode,:,1,b),diff,se,CI,H0_two_samples(electrode,:,2,b),tcrit,df]=limo_yuen_ttest(Y1(1,:,boot_table1{electrode}(:,b)),Y2(1,:,boot_table2{electrode}(:,b))); 
                        % [m,dfe,ci,sd,n,H0_two_samples(electrode,:,1,b),H0_two_samples(electrode,:,2,b)]=limo_ttest(2,Y1(1,:,boot_table1{electrode}(:,b)),Y2(1,:,boot_table2{electrode}(:,b)),0.05); 
                    end
                    clear Y1 Y2
                end
                save (['H0', filesep, boot_name],'H0_two_samples','-v7.3');
            end

            if (tfce ~= 0 && size(data1,1) > 1) % check if tfce is on and if more than one electrode
                mkdir tfce
                tfce_name = sprintf('tfce_two_samples_ttest_parameter_%g',parameter);
                tfce_H0_name = sprintf('tfce_H0_two_samples_ttest_parameter_%g',parameter);
                fprintf('Thresholding Two Sample T-test using TFCE \n');
                tfce_two_samples = limo_tfce(squeeze(two_samples(:,:,4)),channeighbstructmat);
                save(['tfce', filesep, tfce_name], 'tfce_two_samples', '-v7.3'); clear tfce_two_samples;
                % do tfce for the boostrapped data
                fprintf('Thresholding H0 Two Sample T-test using TFCE \n');
                tfce_H0_two_samples = limo_tfce(squeeze(H0_two_samples(:,:,1,:)),channeighbstructmat);
                save(['H0', filesep, tfce_H0_name],'tfce_H0_two_samples', '-v7.3'); clear tfce_H0_two_samples;
            end
            clear H0_two_samples
        end % closes if nboot > 0
        clear two_samples
        disp('two samples t-test done')
        
        
        %--------------------------------------------------------------------------
        %            Paired t-test     //     percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {3}
        
        data1               = varargin{2}; 
        data2               = varargin{3}; 
        parameter           = varargin{4}; 
        nboot               = varargin{5};
        tfce                = varargin{6};
        clear varargin
        
        % ------------------------------------------------
        % make a paired_samples file per parameter (electrodes, frames, [mean value, se, df, t, p])
        paired_samples = NaN(size(data1,1), size(data1,2),5);
        name = sprintf('paired_samples_ttest_parameter_%s',num2str(parameter')');
        
        array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
        for e = 1:size(array,1)
            electrode = array(e);
            fprintf('analyse parameter %s electrode %g',num2str(parameter')', electrode); disp(' ');
            tmp = data1(electrode,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            tmp = data2(electrode,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            [paired_samples(electrode,:,4),paired_samples(electrode,:,1),paired_samples(electrode,:,2),CI,paired_samples(electrode,:,5),tcrit,paired_samples(electrode,:,3)]=limo_yuend_ttest(Y1,Y2); clear Y1 Y2
            % [paired_samples(electrode,:,1),paired_samples(electrode,:,3),ci,sd,n,paired_samples(electrode,:,4),paired_samples(electrode,:,5)]=limo_ttest(1,Y1,Y2,.05); clear Y1 Y2
            % paired_samples(electrode,:,2) = sd./sqrt(n);
        end
        save ([name],'paired_samples', '-v7.3')
        
        % ------------------------------------------------
        if nboot > 0

            bootex = 1;
            boot_name = sprintf('H0_paired_samples_ttest_parameter_%s',num2str(parameter')');
            if exist(['H0', filesep, boot_name, '.mat'], 'file')
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'Yes');
                    bootex = 1;
                else 
                    bootex = 0;
                end
            end
            
            if bootex == 1;
                mkdir H0
                % create a boot one_sample file to store data under H0 
                H0_paired_samples = NaN(size(data1,1), size(data1,2), 2, nboot); % stores T and p values for each boot
                % create centered data to estimate H0
                data1_centered = data1 - repmat(nanmean(data1,3),[1 1 size(data1,3)]);
                data2_centered = data2 - repmat(nanmean(data2,3),[1 1 size(data2,3)]);
                % get boot table
                disp('making boot tables ...')
                boot_table = limo_create_boot_table(data1,nboot);
                save(['H0', filesep, 'boot_table'], 'boot_table')

                % get results under H0
                for e = 1:size(array,1)
                    electrode = array(e);
                    fprintf('bootstrapping electrode %g/%g parameter %s \n',e,size(array,1),num2str(parameter')');
                    tmp = data1_centered(electrode,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    tmp = data2_centered(electrode,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    for b=1:nboot
                        [H0_paired_samples(electrode,:,1,b),diff,se,CI,H0_paired_samples(electrode,:,2,b),tcrit,df]=limo_yuend_ttest(Y1(1,:,boot_table{electrode}(:,b)),Y2(1,:,boot_table{electrode}(:,b))); 
                        % [m,dfe,ci,sd,n,H0_paired_samples(electrode,:,1,b),H0_paired_samples(electrode,:,2,b)]=limo_ttest(1,Y1(1,:,boot_table{electrode}(:,b)),Y2(1,:,boot_table{electrode}(:,b)),0.05); 
                    end
                    clear Y1 Y2
                end
                save (['H0', filesep, boot_name],'H0_paired_samples','-v7.3');
            end

            if (tfce ~= 0 && size(data1,1) > 1) % check if tfce is on and if more than one electrode
                mkdir tfce
                tfce_name = sprintf('tfce_paired_samples_ttest_parameter_%s',num2str(parameter')');
                tfce_H0_name = sprintf('tfce_H0_paired_samples_ttest_parameter_%s',num2str(parameter')');
                fprintf('Thresholding Paired Sample T-test using TFCE \n');
                tfce_paired_samples = limo_tfce(squeeze(paired_samples(:,:,4)),channeighbstructmat);
                save(['tfce', filesep, tfce_name], 'tfce_paired_samples', '-v7.3'); clear tfce_paired_samples;
                % do tfce for the boostrapped data
                fprintf('Thresholding H0 Paired Sample T-test using TFCE \n');
                tfce_H0_paired_samples = limo_tfce(squeeze(H0_paired_samples(:,:,1,:)),channeighbstructmat);
                save(['H0', filesep, tfce_H0_name],'tfce_H0_paired_samples', '-v7.3'); clear tfce_H0_paired_samples;
            end
            clear H0_paired_samples
        end % closes if nboot > 0
        clear paired_samples
        disp('paired t-test done')
        
        %--------------------------------------------------------------------------
        %                    Regression  // percentile bootstrap under H0 and H1
        %--------------------------------------------------------------------------
    case {4}
        
        data                = varargin{2}; 
        regressors          = varargin{3}; % the predictors across subjects like e.g. age
        parameter           = varargin{4}; % the parameters from 1st level matrices the regression is computed on, coded as e.g. 1
        nboot               = varargin{5};
        tfce                = varargin{6};
        clear varargin
        
        % ------------------------------------------------
        % update the LIMO structure
        load LIMO
        LIMO.design.type_of_analysis  = 'Mass-univariate';
        LIMO.data.Cat                 = 0;
        LIMO.data.Cont                = regressors;
        LIMO.data.data_dir            = pwd;
        LIMO.design.fullfactorial     = 0;
        answer = questdlg('zscore regressor(s)?','Regression option','Yes','No','Yes');
        if isempty(answer)
            return
        elseif strcmp(answer,'Yes')
            LIMO.design.zscore            = 1;
        else
            LIMO.design.zscore            = 0;
        end
        LIMO.design.status            = 'to do';
        LIMO.design.method            = 'OLS'; % change to 'IRLS' for robust (to be validated)
        % make design matrix and files
        [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
            LIMO.design.nb_continuous] = limo_design_matrix(data, LIMO,1);
        
        a = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
        if strcmp(a,'Yes')
            save LIMO LIMO
            clear data regressors files LIMO
            limo_eeg(4); disp('regression analysis done');
        else
            return
        end

       
        %--------------------------------------------------------------------------
        %                    N-ways ANOVA / ANCOVA
        %--------------------------------------------------------------------------
    case {5}
        
        data                = varargin{2}; 
        cat                 = varargin{3};
        cont                = varargin{4};
        nboot               = varargin{5};
        tfce                = varargin{6};
        clear varargin

        % ------------------------------------------------
        % update the LIMO structure
        load LIMO
        LIMO.design.type_of_analysis  = 'Mass-univariate';
        LIMO.data.Cat                 = cat;
        LIMO.data.Cont                = cont;
        LIMO.data.data_dir            = pwd;
        if size(cat,2) > 1
            LIMO.design.fullfactorial     = 1;
        else
            LIMO.design.fullfactorial     = 0;
        end
        LIMO.design.zscore            = 1;
        LIMO.design.status            = 'to do';
        LIMO.design.method            = 'OLS'; % change to 'IRLS' for robust (to be validated)
        
        % make design matrix and files
        [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
            LIMO.design.nb_continuous] = limo_design_matrix(data, LIMO,1);
        
        a = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
        if strcmp(a,'Yes')
            save LIMO LIMO; clear data LIMO
            limo_eeg(4); disp('ANOVA/ANCOVA analysis done');
        else
            return
        end
        
        %----------------------------------------------------------------------------------------------
        %                    Repeated Measure ANOVA (multivariate approach) - bootstrap centering data
        %----------------------------------------------------------------------------------------------
    case {6}
        
        data          = varargin{2}; 
        gp_vector     = varargin{3};
        factor_levels = varargin{4};
        nboot         = varargin{5};
        tfce          = varargin{6}; 
        clear varargin
        
        % ------------------------------------------------
        % update the LIMO structure
        load LIMO
        LIMO.design.type_of_analysis  = 'Mass-univariate';
        LIMO.data.Cat                 = gp_vector;
        LIMO.data.Cont                = 0;
        LIMO.data.data_dir            = pwd;
        LIMO.design.nb_conditions     = length(unique(gp_vector));
        LIMO.design.nb_continuous     = 0;
        LIMO.design.fullfactorial     = 0;
        LIMO.design.zscore            = 0;
        LIMO.design.repeated_measure  = factor_levels;
        LIMO.design.method            = 'M'; % = trimmed means = change to 'TM' for robust standard ANOVA (to be validated)
        save LIMO LIMO
        
        % specific stuff for repeated measures
        % from the input we know which case to handle
        if unique(gp_vector) == 1
            % one sample
            if length(factor_levels) ==1
                type = 1;
            elseif length(factor_levels) >1
                type = 2;
            end
        else
            % k samples
            if length(factor_levels) ==1
                type = 3;
            elseif length(factor_levels) >1
                type = 4;
            end
        end
        
         
        % make files to be stored
        if type == 1 % one factor
            C                     = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            tmp_Rep_ANOVA         = NaN(size(data,1),size(data,2),1,2); % store F and p (means can computed via the central tendency and CI button)
            load LIMO;
            LIMO.design.effects{1} = 'Main effect';
            LIMO.design.C{1} = C;
            x = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X = [x ones(size(x,1),1)];
            
        elseif type == 2 % many factor
            C = limo_OrthogContrasts(factor_levels);
            tmp_Rep_ANOVA = NaN(size(data,1),size(data,2),length(C),2); % store F and p for each within factor and interactions
            load LIMO; LIMO.design.C = C; index = length(factor_levels)+1;
            for i= 1:length(factor_levels); LIMO.design.effects{i} = ['Main effect ' num2str(i)]; end
            for i= 2:length(factor_levels); n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1); LIMO.design.effects{index} = ['Interaction ' num2str(n(j,:))]; index = index+1; end; end
            x = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X = [x ones(size(x,1),1)];
            
        elseif type == 3 % one factor within and one factor between
            gp_values = unique(gp_vector); k = length(gp_values); X = NaN(size(gp_vector,1),k+1);
            for g =1:k; X(:,g) = gp_vector == gp_values(g); end; X(:,end) = 1; % design matrix for gp effects
            C                                 = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            tmp_Rep_ANOVA                     = NaN(size(data,1),size(data,2),1,2);
            Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),1,2);
            load LIMO; LIMO.design.C{1} = C; LIMO.design.effects = 'Main effect';
            x = kron(X(:,1:k),eye(prod(factor_levels)));
            LIMO.design.X = [x sum(x,2)]; % just for display
             
        elseif type == 4 % many factors within and one factor between
            gp_values = unique(gp_vector); k = length(gp_values); X = NaN(size(gp_vector,1),k+1);
            for g =1:k; X(:,g) = gp_vector == gp_values(g); end; X(:,end) = 1; % design matrix for gp effects
            C                                 = limo_OrthogContrasts(factor_levels);
            tmp_Rep_ANOVA                     = NaN(size(data,1),size(data,2),length(C),2);
            Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),length(C),2);
            load LIMO; LIMO.design.C = C;
            for i= 1:length(factor_levels);
                LIMO.design.effects{i} = ['Main effect ' num2str(i)];
            end
            index = length(factor_levels)+1;
            for i= 2:length(factor_levels);
                n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1)
                    LIMO.design.effects{index} = ['Interaction ' num2str(n(j,:))];
                    index = index+1;
                end
            end
            x = kron(X(:,1:k),eye(prod(factor_levels)));
            LIMO.design.X = [x sum(x,2)]; % just for display
        end
        
        % check the design with user
        figure('Name','Design matrix'); set(gcf,'Color','w'); imagesc(LIMO.design.X);
        colormap('gray');  title('ANOVA model','FontSize',16);xlabel('regressors');
        ylabel('subjects'); drawnow;
        go = questdlg('start the analysis?');
        if strcmp(go,'No')
            return
        end
        save LIMO LIMO; 
        
        % do the analysis
        array = find(~isnan(data(:,1,1,1)));
        for e = 1:length(array)
            electrode = array(e);
            fprintf('analyse electrode %g/%g\n ...', electrode,size(data,1)); 
            tmp = squeeze(data(electrode,:,:,:));
            Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
            gp = gp_vector(find(~isnan(tmp(1,:,1))),:);
            if type == 3 || type == 4
                XB = X(find(~isnan(tmp(1,:,1))),:);
            end
            
            if type == 1
                if strcmp(LIMO.design.method,'TM')
                    result = limo_robust_rep_anova(Y,gp,factor_levels,C); % trimmed means
                else
                    result = limo_rep_anova(Y,gp,factor_levels,C); % usual means
                end
                tmp_Rep_ANOVA(electrode,:,1,1)             = result.F;
                tmp_Rep_ANOVA(electrode,:,1,2)             = result.p;
                
            elseif type == 2
                if strcmp(LIMO.design.method,'TM')
                    result = limo_robust_rep_anova(Y,gp,factor_levels,C); % trimmed means
                else
                    result = limo_rep_anova(Y,gp,factor_levels,C); % usual means
                end
                tmp_Rep_ANOVA(electrode,:,:,1)             = result.F';
                tmp_Rep_ANOVA(electrode,:,:,2)             = result.p';
                 
            elseif type == 3
                if strcmp(LIMO.design.method,'TM')
                    result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB); % trimmed means
                else
                    result = limo_rep_anova(Y,gp,factor_levels,C,XB); % usual means
                end
                tmp_Rep_ANOVA(electrode,:,1,1)                     = result.repeated_measure.F;
                tmp_Rep_ANOVA(electrode,:,1,2)                     = result.repeated_measure.p;
                Rep_ANOVA_Gp_effect(electrode,:,1)                 = result.gp.F;
                Rep_ANOVA_Gp_effect(electrode,:,2)                 = result.gp.p;
                tmp_Rep_ANOVA_Interaction_with_gp(electrode,:,1)   = result.interaction.F;
                tmp_Rep_ANOVA_Interaction_with_gp(electrode,:,2)   = result.interaction.p;
                
            elseif type == 4
                if strcmp(LIMO.design.method,'TM')
                    result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB); % trimmed means
                else
                    result = limo_rep_anova(Y,gp,factor_levels,C,XB); % usual means
                end
                tmp_Rep_ANOVA(electrode,:,:,1)                     = result.repeated_measure.F';
                tmp_Rep_ANOVA(electrode,:,:,2)                     = result.repeated_measure.p';
                Rep_ANOVA_Gp_effect(electrode,:,1)                 = result.gp.F;
                Rep_ANOVA_Gp_effect(electrode,:,2)                 = result.gp.p;
                tmp_Rep_ANOVA_Interaction_with_gp(electrode,:,:,1) = result.interaction.F';
                tmp_Rep_ANOVA_Interaction_with_gp(electrode,:,:,2) = result.interaction.p';
            end
            clear tmp Y gp result
        end
        
        
        % save stuff
        % ---------
        for i=1:size(tmp_Rep_ANOVA,3)
            name = sprintf('Rep_ANOVA_Factor_%g',i);
            % save each factor effect as F/p values
            % Rep_ANOVA = squeeze(tmp_Rep_ANOVA(:,:,i,:)); 
            % use reshape instead of squeeze in case there is only 1 electrode
            Rep_ANOVA = reshape(tmp_Rep_ANOVA(:,:,i,:),...
                [size(tmp_Rep_ANOVA,1) size(tmp_Rep_ANOVA,2) size(tmp_Rep_ANOVA,4)]);
            save([name],'Rep_ANOVA', '-v7.3'); 
        end
        
        if type == 3 || type ==4
            for i=1:size(tmp_Rep_ANOVA_Interaction_with_gp,3)
                name = sprintf('Rep_ANOVA_Interaction_gp_Factor_%g',i);
                % save each interaction effect as F/p values
                % Rep_ANOVA_Interaction_with_gp = squeeze(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,:));
                Rep_ANOVA_Interaction_with_gp = reshape(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,:),...
                    [size(tmp_Rep_ANOVA_Interaction_with_gp,1) size(tmp_Rep_ANOVA_Interaction_with_gp,2) size(tmp_Rep_ANOVA_Interaction_with_gp,4)]); 
                save([name],'Rep_ANOVA_Interaction_with_gp', '-v7.3'); clear Rep_ANOVA_Interaction_with_gp;
            end
            save Rep_ANOVA_Gp_effect Rep_ANOVA_Gp_effect -v7.3; % always only 1 effect
        end
        
        
        % ----------------------------------------------------------------
        if nboot > 0
            
            bootex = 1;
            name = 'Rep_ANOVA';
            boot_name = sprintf('H0%sboot_H0_%s', filesep, name);
            if exist(['H0', filesep, boot_name, '.mat'], 'file')
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'Yes');
                    bootex = 1;
                else 
                    bootex = 0;
                end
            end
            
            
            if bootex == 1;
                mkdir H0
                % create files to store bootstrap under H1 and H0
                disp('making bootstrap files ...')
                if type ==1
                    tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),1,2,nboot);
                elseif type == 2
                    tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),length(C),2,nboot);
                elseif type == 3
                    tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),1,2,nboot);
                    boot_H0_Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2,nboot);
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),1,2,nboot);
                else
                    tmp_boot_H0_Rep_ANOVA                     = NaN(size(data,1),size(data,2),length(C),2,nboot);
                    boot_H0_Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2,nboot);
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),length(C),2,nboot);
                end
            
                % the data have to be centered (H0) for each cell
                centered_data = NaN(size(data,1),size(data,2),size(data,3),size(data,4));           
                nb_conditions = prod(factor_levels);

                if type ==1 || type ==2
                    for condition=1:nb_conditions
                        avg = repmat(nanmean(data(:,:,:,condition),3),[1 1 size(data,3)]);
                        centered_data(:,:,:,condition) = data(:,:,:,condition) - avg; clear avg
                    end
                else
                    for gp=1:LIMO.design.nb_conditions % = nb gp
                        gp_index = find(LIMO.data.Cat == gp);
                        for condition=1:nb_conditions
                            avg = repmat(nanmean(data(:,:,gp_index,condition),3),[1 1 length(gp_index)]);
                            centered_data(:,:,gp_index,condition) = data(:,:,gp_index,condition) - avg; clear avg
                        end
                    end
                end
                save(['H0', filesep, 'centered_data'], 'centered_data', '-v7.3');

                % create an index to use across all electrodes and frames
                % (different per gp but identical across conditions)
                disp('making random table...')
                % boot_table = limo_create_boot_table(squeeze(data(:,:,:,1)),nboot);
                boot_table = limo_create_boot_table(reshape(data(:,:,:,1),[size(data,1) size(data,2) size(data,3)]),nboot);
                save(['H0', filesep, 'boot_table'], 'boot_table')
            
                % compute bootstrap under H0 for F and p
                for B=1:nboot
                    fprintf('Repeated Measures ANOVA bootstrap %g \n ...', B);
                    for e = 1:length(array)
                        electrode = array(e);
                        
                        tmp = squeeze(centered_data(electrode,:,boot_table{electrode}(:,B),:));
                        Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                        gp = gp_vector(find(~isnan(tmp(1,:,1))));
                        if type == 3 || type == 4
                            XB = X(find(~isnan(tmp(1,:,1))));
                        end
                        
                        if type == 1
                            if strcmp(LIMO.design.method,'TM')
                                result = limo_robust_rep_anova(Y,gp,factor_levels,C);
                            else
                                result = limo_rep_anova(Y,gp,factor_levels,C);
                            end
                            tmp_boot_H0_Rep_ANOVA(electrode,:,1,1,B)         = result.F;
                            tmp_boot_H0_Rep_ANOVA(electrode,:,1,2,B)         = result.p;
                        elseif type == 2
                            if strcmp(LIMO.design.method,'TM')
                                result = limo_robust_rep_anova(Y,gp,factor_levels,C);
                            else
                                result = limo_rep_anova(Y,gp,factor_levels,C);
                            end
                            tmp_boot_H0_Rep_ANOVA(electrode,:,:,1,B)         = result.F';
                            tmp_boot_H0_Rep_ANOVA(electrode,:,:,2,B)         = result.p';
                        elseif type == 3
                            if strcmp(LIMO.design.method,'TM')
                                result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
                            else
                                result = limo_rep_anova(Y,gp,factor_levels,C,XB);
                            end
                            tmp_boot_H0_Rep_ANOVA(electrode,:,1,1,B)                     = result.repeated_measure.F;
                            tmp_boot_H0_Rep_ANOVA(electrode,:,1,2,B)                     = result.repeated_measure.p;
                            H0_Rep_ANOVA_Gp_effect(electrode,:,1,B)                      = result.gp.F;
                            H0_Rep_ANOVA_Gp_effect(electrode,:,2,B)                      = result.gp.p;
                            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,1,1,B) = result.interaction.F;
                            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,1,2,B) = result.interaction.p;
                        elseif type == 4
                            if strcmp(LIMO.design.method,'TM')
                                result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
                            else
                                result = limo_rep_anova(Y,gp,factor_levels,C,XB);
                            end
                            tmp_boot_H0_Rep_ANOVA(electrode,:,:,1,B)                     = result.repeated_measure.F';
                            tmp_boot_H0_Rep_ANOVA(electrode,:,:,2,B)                     = result.repeated_measure.p';
                            H0_Rep_ANOVA_Gp_effect(electrode,:,1,B)                      = result.gp.F;
                            H0_Rep_ANOVA_Gp_effect(electrode,:,2,B)                      = result.gp.p;
                            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,:,1,B) = result.interaction.F';
                            tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(electrode,:,:,2,B) = result.interaction.p';
                            clear y result
                        end
                        clear XB Y gp tmp
                    end
                end
                
                % save 
                for i=1:size(tmp_boot_H0_Rep_ANOVA,3)
                    name = sprintf('H0_Rep_ANOVA_Factor_%g',i);
                    H0_Rep_ANOVA = squeeze(tmp_boot_H0_Rep_ANOVA(:,:,i,:,:)); % save each factor effect as F/p/nboot values
                    save(['H0', filesep, name],'H0_Rep_ANOVA', '-v7.3'); 
                end
                
                if type == 3 || type ==4
                    for i=1:size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp,3)
                        name = sprintf('H0_Rep_ANOVA_Interaction_gp_Factor_%g',i);
                        H0_Rep_ANOVA_Interaction_with_gp = squeeze(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,i,:)); % save each interaction effect as F/p values
                        save(['H0', filesep, name],'H0_Rep_ANOVA_Interaction_with_gp', '-v7.3'); clear H0_Rep_ANOVA_Interaction_with_gp;
                    end
                    
                    save(['H0', filesep, 'H0_Rep_ANOVA_Gp_effect'], 'H0_Rep_ANOVA_Gp_effect', '-v7.3'); 
                end
            end
            
            % ------------------------- TFCE ---------------
            if (tfce ~= 0 && size(data,1) > 1) % check if tfce is on and if more than one electrode
                mkdir tfce; 
                
                fprintf('Thresholding Rep ANOVA using TFCE \n');
                for i=1:size(tmp_Rep_ANOVA,3)
                    fprintf('analyzing effect %g \n',i)
                    tfce_name = sprintf('tfce%stfce_Rep_ANOVA_Factor_%g',filesep,i);
                    tfce_Rep_ANOVA = limo_tfce(squeeze(tmp_Rep_ANOVA(:,:,i,1)),channeighbstructmat);
                    save(tfce_name, 'tfce_Rep_ANOVA'); clear tfce_Rep_ANOVA 
                end
                
                if type == 3 || type == 4
                    % gp effect
                    fprintf('analyzing gp effect \n')
                    tfce_name = sprintf('tfce%stfce_Rep_ANOVA_Gp_effect',filesep);
                    tfce_Rep_ANOVA_Gp_effect = limo_tfce(squeeze(Rep_ANOVA_Gp_effect(:,:,1)),channeighbstructmat);
                    save(tfce_name, 'tfce_Rep_ANOVA_Gp_effect'); clear tfce_Rep_ANOVA_Gp_effect Rep_ANOVA_Gp_effect
                    % interactions
                    for i=1:size(tmp_Rep_ANOVA_Interaction_with_gp,3)
                        fprintf('analyzing interaction effect %g \n',i)
                        tfce_name = sprintf('tfce%stfce_Rep_ANOVA_Interaction_gp_Factor_%g',filesep,i);
                        tfce_Rep_ANOVA_Interaction_with_gp = limo_tfce(squeeze(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,1)),channeighbstructmat);
                        save(tfce_name, 'tfce_Rep_ANOVA_Interaction_with_gp'); clear tfce_Rep_ANOVA_Interaction_with_gp tmp_Rep_ANOVA_Interaction_with_gp
                    end
                end
                
                % do tfce for the boostrapped data
                % -------------------------------
                fprintf('Thresholding bootstrapped Rep ANOVA using TFCE \n');
                for i=1:size(tmp_Rep_ANOVA,3)
                    fprintf('analyzing effect %g \n',i);
                    tfce_name = sprintf('H0%stfce_H0_Rep_ANOVA_Factor_%g',filesep,i);
                    tfce_H0_Rep_ANOVA = limo_tfce(squeeze(tmp_boot_H0_Rep_ANOVA(:,:,i,1,:)),channeighbstructmat);
                    save(tfce_name, 'tfce_H0_Rep_ANOVA'); clear tfce_H0_Rep_ANOVA 
                end
                
                if type == 3 || type == 4
                    % gp effect
                    fprintf('analyzing gp effect \n')
                    tfce_name = sprintf('H0%stfce_H0_Rep_ANOVA_Gp_effect',filesep);
                    tfce_H0_Rep_ANOVA_Gp_effect = limo_tfce(squeeze(H0_Rep_ANOVA_Gp_effect(:,:,1,:)),channeighbstructmat);
                    save(tfce_name, 'H0_tfce_Rep_ANOVA_Gp_effect'); clear tfce_H0_Rep_ANOVA_Gp_effect H0_Rep_ANOVA_Gp_effect
                    % interactions
                    for i=1:size(tmp_Rep_ANOVA_Interaction_with_gp,3)
                        fprintf('analyzing interaction effect %g \n',i)
                        tfce_name = sprintf('H0%stfce_H0_Rep_ANOVA_Interaction_gp_Factor_%g',filesep,i);
                        tfce_H0_Rep_ANOVA_Interaction_with_gp = limo_tfce(squeeze(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,i,1,:)),channeighbstructmat);
                        save(tfce_name, 'tfce_H0_Rep_ANOVA_Interaction_with_gp'); clear tfce_H0_Rep_ANOVA_Interaction_with_gp 
                    end
                end
            end
        end
        disp('Repeated Measures ANOVA done')
end




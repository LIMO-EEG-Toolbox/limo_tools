function LIMOPath = limo_random_robust(varargin)

% This function makes the result files for the random effects of various tests
% as well as organizes and makes files for boostrap. It is interfaced with
% limo_random_effect which itself interfaces with the user to select and pass
% the data in the appropriate format. Limo_random_robust calls low level
% functions to perform the actual computation and add the trimmed mean or mean
% of parameters for the correponding test (helps for vizualizing effects)
%
% FORMAT LIMOPath = limo_random_robust(test,data,label,LIMO)
%
% INPUTS
%
% limo_random_robust(1,y,parameter number,LIMO)
%                    1 = a one-sample t-test
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                      = the name of the Yr file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(2,y1,y2,parameter number,LIMO)
%                    2 = two samples t-test
%                    y1 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y1r file
%                    y2 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y2r file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(3,y1,y2,parameter number,LIMO)
%                    3 = paired t-test
%                    y1 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y1r file
%                    y2 = data (dim channels, time or freq, subjects)
%                       = data (dim channels, freq, time, subjects)
%                       = the name of the Y2r file
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(4,y,X,parameter number,LIMO)
%                    4 = regression analysis
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                    X = continuous regressor(s)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%
% limo_random_robust(5,y,cat,cont,LIMO,'go',option)
%                    5 = N-way ANOVA/ANCOVA
%                    y = data (dim channels, time or freq, subjects)
%                      = data (dim channels, freq, time, subjects)
%                      = the name of the Yr file
%                    cat = categorical variable(s)
%                    cont = continuous regressors (covariates)
%                    LIMO the basic structure with data, design and channel info
%                    'go' is optional and prompt or not the design
%                         options are 'yes' (prompt, default),
%                         or 'no' (no prompt, usuful for scripting)
%
% limo_random_robust(6,y,gp,factor_levels,LIMO,'go',option)
%                    6 = Repeated measures ANOVA/ANCOVA using multivariate approach
%                    y = data (dim channels, time or freq, subjects, measures)
%                      = data (dim channels, freq, time, subjects, measures)
%                      = the name of the Yr file
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    LIMO the basic structure with data, design and channel info
%                         or the full name of the LIMO file
%                    'go' is optional and prompt or not the pseudo-design
%                         options are 'yes' (prompt, default),
%                         or 'no' (no prompt, usuful for scripting)
%
% OUTPUT
% write on the disk matrices correponding to the test (Yr and LIMO.mat are generated in limo_random_select,
% and for Regression, ANOVA, the LIMO.mat structure is updated)
%
% 1 one_sample_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_one_sample_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)
%
% 2 two_samples_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_two_samples_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)
%
% 3 paired_samples_parameter_X (channels, frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_paired_samples_ttest_parameter_X (channels, frames, [T values under H0, p values under H0], LIMO.design.bootstrap)
%
% 4 R2 (channels, frames [time, freq or freq-time], [F p values])
%   H0_R2 (channels, frames, [F p], LIMO.design.bootstrap)
%   Covariate_effect_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_Covariate_effect_X (channels, frames, [F p], LIMO.design.bootstrap)
%
% 5 Condition_effect_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_Condition_effect_X (channels, frames, [F p], LIMO.design.bootstrap)
%   Covariate_effect_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_Covariate_effect_X (channels, frames, [F p], LIMO.design.bootstrap)
%
% 6 Rep_ANOVA_Factor_X (channels, frames [time, freq or freq-time], [F p values])
%   Rep_ANOVA_Gp_effect (channels, frames [time, freq or freq-time], [F p values])
%   Rep_ANOVA_Interaction_gp_Factor_X (channels, frames [time, freq or freq-time], [F p values])
%   H0_XXXXX same as above, including LIMO.design.bootstrap on the last dimension
%
% LIMOPath = LIMO.dir or [] if failed
%
% See also LIMO_TRIMCI LIMO_YUEN_TTEST LIMO_YUEND_TTEST LIMO_ROBUST_1WAY_ANOVA
% LIMO_GLM1 LIMO_EEG(4) LIMO_EEG_TF(4) LIMO_REP_ANOVA LIMO_CREATE_BOOT_TABLE
% ------------------------------
%  Copyright (C) LIMO Team 2020

%% inputs checks
LIMOPath = [];
if nargin == 0
    help limo_random_robust
    return
else
    type  = varargin{1};
    if type <1 || type > 6
        error('type argument must be between 1 and 6')
    end
end

%% start

switch type
    %--------------------------------------------------------------------------
    % One Sample t-test // bootstrap-t method
    %--------------------------------------------------------------------------
    case {1}
        
        if ischar(varargin{2})
            data      = load(varargin{2});
            data      = data.(cell2mat(fieldnames(data)));
        else
            data  = varargin{2};
        end
        parameter = varargin{3};
        if ischar(varargin{4})
            LIMO = load(varargin{4});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{4})
            LIMO = varargin{4};
        end
        clear varargin
        cd(LIMO.dir);
        
        if strcmp(LIMO.Analysis,'Time-Frequency') 
            data = limo_tf_4d_reshape(data);
        end
        
        % ------------------------------------------------
        % check the data structure
        for e=1:size(data,1)
            tmp = isnan(data(e,1,:));
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects - analysis aborded']);
                return
            end
        end
        clear tmp
        
        % ------------------------------------------------
        % make a one_sample file per parameter (channels, frames, [mean value, se, df, t, p])
        one_sample = NaN(size(data,1), size(data,2), 5);
        name       = sprintf('one_sample_ttest_parameter_%g',parameter);
        
        for channel = 1:size(data,1) % run per channel because we have to remove NaNs
            fprintf('analyse parameter %g channel %g \n',parameter, channel);
            tmp = data(channel,:,:);
            if nansum(tmp(1,:)) == 0
                error('there is at least one empty channel using your expected chanlocs')
            else
                Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
            end
            
            if ~isfield(LIMO.design,'method')
                LIMO.design.method = 'Trimmed Mean';
            end
            
            if strcmpi(LIMO.design.method,'Trimmed Mean')
                [one_sample(channel,:,4),one_sample(channel,:,1),~,one_sample(channel,:,2), ...
                    one_sample(channel,:,5),~,one_sample(channel,:,3)] = limo_trimci(Y);
            elseif strcmpi(LIMO.design.method,'Mean')
                [one_sample(channel,:,1),one_sample(channel,:,3),~,sd,n, ...
                    one_sample(channel,:,4),one_sample(channel,:,5)] = limo_ttest(1,Y,0,5/100);
                one_sample(channel,:,2) = sd./sqrt(n);
            else
                error('unrecognized LIMO.design.method: %s',LIMO.design.method)
            end
            clear tmp Y
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            one_sample = limo_tf_4d_reshape(one_sample);
        end
        save (name,'one_sample', '-v7.3')
        LIMOPath = LIMO.dir;
        
        % ------------------------------------------------
        % Bootstrap
        if LIMO.design.bootstrap > 0
            
            bootex = 1;
            boot_name = sprintf('H0_one_sample_ttest_parameter_%g',parameter);
            if exist(['H0', filesep, boot_name, '.mat'], 'file')
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    bootex = 1;
                else
                    bootex = 0;
                end
            end
            
            if bootex == 1
                mkdir H0
                % create a boot one_sample file to store data under H0 and H1
                H0_one_sample = NaN(size(data,1), size(data,2),2,LIMO.design.bootstrap); % stores T and p values for each boot under H0
                % create centered data to estimate H0
                centered_data = data - repmat(limo_trimmed_mean(data),[1 1 size(data,3)]);
                % centered_data = data - repmat(nanmean(data,3),[1 1 size(data,3)]);
                % get boot table
                disp('making boot table ...')
                boot_table = limo_create_boot_table(data,LIMO.design.bootstrap);
                save(['H0', filesep, 'boot_table'], 'boot_table')
                
                % get results under H0
                for channel = 1:size(data,1)
                    fprintf('bootstrap: channel %g parameter %g \n',channel,parameter);
                    tmp = centered_data(channel,:,:);
                    Y   = tmp(1,:,find(~isnan(tmp(1,1,:))));
                    if strcmpi(LIMO.design.method,'Trimmed Mean')
                        parfor b=1:LIMO.design.bootstrap
                            [t{b},~,~,~,p{b},~,~] = limo_trimci(Y(1,:,boot_table{channel}(:,b)));
                        end
                    elseif strcmpi(LIMO.design.method,'Mean')
                        parfor b=1:LIMO.design.bootstrap
                            [~,~,~,~,~,t{b},p{b}] = limo_ttest(1,Y(1,:,boot_table{channel}(:,b)),0,5/100);
                        end
                    end
                    
                    for b=1:LIMO.design.bootstrap
                        H0_one_sample(channel,:,1,b) = t{b};
                        H0_one_sample(channel,:,2,b) = p{b};
                    end
                    clear tmp Y
                end % closes for channel
                
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    H0_one_sample = limo_tf_5d_reshape(H0_one_sample);
                end
                save (['H0', filesep, boot_name],'H0_one_sample','-v7.3');
            end
        end
        
        if LIMO.design.tfce ~= 0
            limo_tfce_handling(fullfile(LIMO.dir,name));
            LIMO.design.tfce = 1;
            save(fullfile(LIMO.dir,'LIMO.mat'));
        end
        disp('one sample t-test done')
            
        %--------------------------------------------------------------------------
        % Two Samples t-test // percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {2}
        
        if ischar(varargin{2})
            data1  = load(varargin{2});
            data1  = data1.(cell2mat(fieldnames(data1)));
        else
            data1  = varargin{2};
        end
        if ischar(varargin{3})
            data2  = load(varargin{3});
            data2  = data2.(cell2mat(fieldnames(data2)));
        else
            data2  = varargin{3};
        end
        parameter = varargin{4};
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        clear varargin
       
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            data1 = limo_tf_4d_reshape(data1);
            data2 = limo_tf_4d_reshape(data2);
        end
        
        % ------------------------------------------------
        % check the data structure
        if size(data1,1) ~= size(data2,1)
            error(['groups have a different number of ' LIMO.Type])
        end
        
        for e=1:size(data1,1)
            tmp = isnan(data1(e,1,:));
            tmp2 = isnan(data2(e,1,:));
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty in group 1 - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects in group 1 - analysis aborded']);
                return
            elseif length(tmp2) == sum(isnan(tmp2))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty in group 2 - analysis aborded']);
                return
            elseif (length(tmp2) - sum(isnan(tmp2))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects in group 2 - analysis aborded']);
                return
            end
        end
        clear tmp tmp2
        
        % ------------------------------------------------
        % make a two_samples file per parameter (channels, frames, [mean value, se, df, t, p])
        two_samples = NaN(size(data1,1), size(data1,2),5);
        name = sprintf('two_samples_ttest_parameter_%g_%g',parameter);
        
        array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
        for e = 1:size(array,1)
            channel = array(e);
            fprintf('analyse parameter %g channel %g',parameter, channel); disp(' ');
            tmp = data1(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            tmp = data2(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            [two_samples(channel,:,4),two_samples(channel,:,1),two_samples(channel,:,2),...
                CI,two_samples(channel,:,5),tcrit,two_samples(channel,:,3)]=limo_yuen_ttest(Y1,Y2); clear Y1 Y2
            % [two_samples(channel,:,1),two_samples(channel,:,3),ci,sd,n,two_samples(channel,:,4),two_samples(channel,:,5)]=limo_ttest(2,Y1,Y2,.05);
            % sd = sd.^2; a = sd(1,:)./size(Y1,3); b = sd(1,:)./size(Y2,3);
            % two_samples(channel,:,2) = sqrt(a + b); clear Y1 Y2
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            two_samples = limo_tf_4d_reshape(two_samples);
        end
        save (name,'two_samples', '-v7.3')
        LIMOPath = LIMO.dir;
        
        % ------------------------------------------------
        % compute the null
        if LIMO.design.bootstrap > 0
            
            bootex = 1;
            boot_name = sprintf('H0_two_samples_ttest_parameter_%g_%g',parameter);
            if exist(['H0', filesep, boot_name, '.mat'], 'file')
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    bootex = 1;
                else
                    bootex = 0;
                end
            end
            
            if bootex == 1
                mkdir H0
                % create a boot one_sample file to store data under H0
                H0_two_samples = NaN(size(data1,1), size(data1,2), 2, LIMO.design.bootstrap); % stores T and p values for each boot
                % create centered data to estimate H0
                data1_centered = data1 - repmat(limo_trimmed_mean(data1),[1 1 size(data1,3)]);
                data2_centered = data2 - repmat(limo_trimmed_mean(data2),[1 1 size(data2,3)]);
                % data1_centered = data1 - repmat(nanmean(data1,3),[1 1 size(data1,3)]);
                % data2_centered = data2 - repmat(nanmean(data2,3),[1 1 size(data2,3)]);
                % get boot table
                disp('making boot tables ...')
                boot_table1 = limo_create_boot_table(data1,LIMO.design.bootstrap);
                boot_table2 = limo_create_boot_table(data2,LIMO.design.bootstrap);
                save(['H0', filesep, 'boot_table1'], 'boot_table1')
                save(['H0', filesep, 'boot_table2'], 'boot_table2')
                
                % get results under H0
                for e = 1:size(array,1)
                    channel = array(e);
                    fprintf('bootstrapping channel %g/%g \n',e,size(array,1));
                    tmp = data1_centered(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    tmp = data2_centered(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    if exist('parfor','file') ~=0
                        parfor b=1:LIMO.design.bootstrap
                            [t{b},~,~,~,p{b},~,~]=limo_yuen_ttest(Y1(1,:,boot_table1{channel}(:,b)),Y2(1,:,boot_table2{channel}(:,b)));
                        end
                        
                        for b=1:LIMO.design.bootstrap
                            H0_two_samples(channel,:,1,b) = t{b};
                            H0_two_samples(channel,:,2,b) = p{b};
                        end
                        clear t p
                        
                    else
                        for b=1:LIMO.design.bootstrap
                            [H0_two_samples(channel,:,1,b),diff,se,CI,H0_two_samples(channel,:,2,b),tcrit,df]=limo_yuen_ttest(Y1(1,:,boot_table1{channel}(:,b)),Y2(1,:,boot_table2{channel}(:,b)));
                            % [m,dfe,ci,sd,n,H0_two_samples(channel,:,1,b),H0_two_samples(channel,:,2,b)]=limo_ttest(2,Y1(1,:,boot_table1{channel}(:,b)),Y2(1,:,boot_table2{channel}(:,b)),0.05);
                        end
                    end
                    clear Y1 Y2
                end
                
                
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    H0_two_samples = limo_tf_5d_reshape(H0_two_samples);
                end
                save (['H0', filesep, boot_name],'H0_two_samples','-v7.3');
            end
        end % closes if LIMO.design.bootstrap > 0

        if LIMO.design.tfce ~= 0
            limo_tfce_handling(fullfile(LIMO.dir,name));
            LIMO.design.tfce = 1;
            save(fullfile(LIMO.dir,'LIMO.mat'));
        end
        disp('two samples t-test done')
        
        
        %--------------------------------------------------------------------------
        % Paired t-test // percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {3}
        
        if ischar(varargin{2})
            data1  = load(varargin{2});
            data1  = data1.(cell2mat(fieldnames(data1)));
        else
            data1  = varargin{2};
        end
        if ischar(varargin{3})
            data2  = load(varargin{3});
            data2  = data2.(cell2mat(fieldnames(data2)));
        else
            data2  = varargin{3};
        end
        parameter = varargin{4};
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        clear varargin
        
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            data1 = limo_tf_4d_reshape(data1);
            data2 = limo_tf_4d_reshape(data2);
        end
        clear varargin
        
        % ------------------------------------------------
        % check the data structure
        if size(data1,1) ~= size(data2,1)
            error(['samples have a different number of ' LIMO.Type ',not a paired t-tests'])
        end
        
        for e=1:size(data1,1)
            tmp = isnan(data1(e,1,:));
            tmp2 = isnan(data2(e,1,:));
            if length(tmp) ~= length(isnan(tmp2))
                errordlg([LIMO.Type ' ' num2str(e) ' has unpaired data - analysis aborded, not a paired t-test']);
                return
            elseif length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects - analysis aborded']);
                return
            end
        end
        clear tmp tmp2
        
        % ------------------------------------------------
        % make a paired_samples file per parameter (channels, frames, [mean value, se, df, t, p])
        paired_samples = NaN(size(data1,1), size(data1,2),5);
        name = sprintf('paired_samples_ttest_parameter_%s',num2str(parameter')');
        
        array = intersect(find(~isnan(data1(:,1,1))),find(~isnan(data2(:,1,1))));
        for e = 1:size(array,1)
            channel = array(e);
            fprintf('analyse parameter %s channel %g',num2str(parameter')', channel); disp(' ');
            tmp = data1(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            tmp = data2(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
            [paired_samples(channel,:,4),paired_samples(channel,:,1),paired_samples(channel,:,2),...
                CI,paired_samples(channel,:,5),tcrit,paired_samples(channel,:,3)]=limo_yuend_ttest(Y1,Y2); clear Y1 Y2
            % [paired_samples(channel,:,1),paired_samples(channel,:,3),ci,sd,n,paired_samples(channel,:,4),paired_samples(channel,:,5)]=limo_ttest(1,Y1,Y2,.05); clear Y1 Y2
            % paired_samples(channel,:,2) = sd./sqrt(n);
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
            paired_samples = limo_tf_4d_reshape(paired_samples);
        end
        save (name,'paired_samples', '-v7.3')
        LIMOPath = LIMO.dir; 
        
        % ------------------------------------------------
        if LIMO.design.bootstrap > 0
            
            bootex = 1;
            boot_name = sprintf('H0_paired_samples_ttest_parameter_%s',num2str(parameter')');
            if exist(['H0', filesep, boot_name, '.mat'], 'file')
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'Yes')
                    bootex = 1;
                else
                    bootex = 0;
                end
            end
            
            if bootex == 1
                mkdir H0
                % create a boot one_sample file to store data under H0
                H0_paired_samples = NaN(size(data1,1), size(data1,2), 2, LIMO.design.bootstrap); % stores T and p values for each boot
                % create centered data to estimate H0
                data1_centered = data1 - repmat(limo_trimmed_mean(data1),[1 1 size(data1,3)]);
                data2_centered = data2 - repmat(limo_trimmed_mean(data2),[1 1 size(data2,3)]);
                % data1_centered = data1 - repmat(nanmean(data1,3),[1 1 size(data1,3)]);
                % data2_centered = data2 - repmat(nanmean(data2,3),[1 1 size(data2,3)]);
                % get boot table
                disp('making boot table ...')
                boot_table = limo_create_boot_table(data1,LIMO.design.bootstrap);
                save(['H0', filesep, 'boot_table'], 'boot_table')
                
                % get results under H0
                for e = 1:size(array,1)
                    channel = array(e);
                    fprintf('bootstrapping channel %g/%g parameter %s \n',e,size(array,1),num2str(parameter')');
                    tmp = data1_centered(channel,:,:); Y1 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    tmp = data2_centered(channel,:,:); Y2 = tmp(1,:,find(~isnan(tmp(1,1,:)))); clear tmp
                    if exist('parfor','file') ~=0
                        parfor b=1:LIMO.design.bootstrap
                            [t{b},~,~,~,p{b},~,~]=limo_yuend_ttest(Y1(1,:,boot_table{channel}(:,b)),Y2(1,:,boot_table{channel}(:,b)));
                        end
                        
                        for b=1:LIMO.design.bootstrap
                            H0_paired_samples(channel,:,1,b) = t{b};
                            H0_paired_samples(channel,:,2,b) = p{b};
                        end
                        clear t p
                        
                    else
                        for b=1:LIMO.design.bootstrap
                            [H0_paired_samples(channel,:,1,b),diff,se,CI,H0_paired_samples(channel,:,2,b),tcrit,df]=limo_yuend_ttest(Y1(1,:,boot_table{channel}(:,b)),Y2(1,:,boot_table{channel}(:,b)));
                            % [m,dfe,ci,sd,n,H0_paired_samples(channel,:,1,b),H0_paired_samples(channel,:,2,b)]=limo_ttest(1,Y1(1,:,boot_table{channel}(:,b)),Y2(1,:,boot_table{channel}(:,b)),0.05);
                        end
                    end
                    clear Y1 Y2
                end
                
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    H0_paired_samples = limo_tf_5d_reshape(H0_paired_samples);
                end
                save (['H0', filesep, boot_name],'H0_paired_samples','-v7.3');
            end
        end % closes if LIMO.design.bootstrap > 0
        
        if LIMO.design.tfce ~= 0
            limo_tfce_handling(fullfile(LIMO.dir,name));
            LIMO.design.tfce = 1;
            save(fullfile(LIMO.dir,'LIMO.mat'));
        end
        disp('Paired t-test done')
        
        %------------------------------------------------------------------
        % Regression // percentile bootstrap under H0
        %------------------------------------------------------------------
    case {4}
        
        if ischar(varargin{2})
            data = load(varargin{2});
            data = data.cell2mat(fieldnames(data));
        else
            data = varargin{2};
        end
        regressors = varargin{3}; % the predictors across subjects like e.g. age
        parameter  = varargin{4}; % the parameters from 1st level matrices the regression is computed on (just for name)
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        
        if nargin >5
            for in = 6:2:nargin
                if strcmpi(varargin{in},'zscore')
                    answer = varargin{in+1};
                elseif strcmpi(varargin{in},'go')
                    go = varargin{in+1};
                end
            end
        else
            go = 'No';
        end
        clear varargin
        cd(LIMO.dir);
        

        % ------------------------------------------------
        % check the data structure
        for e=1:size(data,1)
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp = isnan(data(e,1,1,:));
            else
                tmp = isnan(data(e,1,:));
            end
            
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 3
                errordlg([LIMO.Type ' ' num2str(e) ' has less than 3 subjects - analysis aborded']);
                return
            elseif (length(tmp) - sum(isnan(tmp))) < 6
                warndlg([LIMO.Type ' ' num2str(e) ' has less than 6 subjects - regression results will likely be biased']);
            end
        end
        
        % ------------------------------------------------
        % update the LIMO structure
        LIMO.data.Cat                = 0;
        LIMO.data.Cont               = regressors;
        LIMO.data.data_dir           = pwd;
        LIMO.design.type_of_analysis = 'Mass-univariate';
        LIMO.design.fullfactorial    = 0;
        LIMO.design.status           = 'to do';
        LIMO.design.method           = 'OLS' ; % 'IRLS' beter but H0 boot too conservative
        
        if ~exist('answer','var') 
            answer = questdlg('zscore regressor(s)?','Regression option','Yes','No','Yes');
        end
        
        if isempty(answer)
            return
        elseif strcmpi(answer,'Yes')
            LIMO.design.zscore = 1;
        else
            LIMO.design.zscore = 0;
        end
        save LIMO LIMO
        
        % make design matrix and files
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix_tf(data, LIMO,1);
        else
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix(data, LIMO,1);
        end
        
        % ------------------------------------------------
        % do the analysis
        if strcmpi(go,'no')
            go = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
            close('LIMO design')
        end
                
        if strcmpi(go,'Yes')
            save LIMO LIMO
            if nargout ~= 0, LIMOPath = fullfile(pwd,'LIMO.mat'); end
            clear data regressors files; limo_eeg(4); disp('regression analysis done');
        else
            return
        end
                
        %--------------------------------------------------------------------------
        % N-ways ANOVA / ANCOVA
        %--------------------------------------------------------------------------
    case {5}
        
        data  = varargin{2};
        if ischar(data)
            data = load(data);
            data = data.(cell2mat(filenames(data)));
        end
        cat   = varargin{3}; 
        cont  = varargin{4}; 
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        cd(LIMO.dir);
        
        if nargin ==7
            if strcmpi(varargin{6},'go')
                go = varargin{7};
            end
        end
        clear varargin
        cd(LIMO.dir);        
        
        % ------------------------------------------------
        % check the data structure
        for e=1:size(data,1)
            if strcmp(LIMO.Analysis,'Time-Frequency')
                tmp = isnan(data(e,1,1,:));
            else
                tmp = isnan(data(e,1,:));
            end
            
            if length(tmp) == sum(isnan(tmp))
                errordlg([LIMO.Type ' ' num2str(e) ' is empty - analysis aborded']);
                return
            end
        end
        
        % ------------------------------------------------
        % update the LIMO structure
        LIMO.design.type_of_analysis  = 'Mass-univariate';
        LIMO.data.Cat                 = cat;
        LIMO.data.Cont                = cont;
        LIMO.data.data_dir            = pwd;
        LIMO.design.zscore            = 1;
        if size(cat,2) > 1
            LIMO.design.fullfactorial = 1;
        else
            LIMO.design.fullfactorial = 0;
        end
        
        % make design matrix and files
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix_tf(data, LIMO,1);
        else
            [LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions,...
                LIMO.design.nb_continuous] = limo_design_matrix(data, LIMO,1);
        end
        save LIMO LIMO
        
        % ------------------------------------------------
        % do the analysis
        if ~exist('go','var')
            go = questdlg('run the analysis?','Start GLM analysis','Yes','No','Yes');
            close('LIMO design')
        end
        
        if strcmpi(go,'Yes')
            if LIMO.design.fullfactorial == 0 && LIMO.design.nb_continuous == 0
                if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                    data = limo_tf_4d_reshape(data);
                end
                Yhat             = NaN(size(data));
                Condition_effect = NaN(size(data,1),size(data,2),2);
                array            = find(~isnan(data(:,1,1)));
                for e=1:size(array,1)
                    channel = array(e); fprintf('processing channel %g \n,',channel);
                    [Condition_effect(channel,:,1), Condition_effect(channel,:,2),Yhat(channel,:,:)] = ...
                        limo_robust_1way_anova(squeeze(data(channel,:,:)),LIMO.design.X(:,1:end-1),20); % no intercept in this model
                end
                delete(fullfile(LIMO.dir,'Betas.mat')); % no betas here
                delete(fullfile(LIMO.dir,'R2.mat'));    % no R2
                
                if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                    Condition_effect = limo_tf_4d_reshape(Condition_effect);
                    Yhat = limo_tf_4d_reshape(Yhat); % these are the trimmed mean
                    data = limo_tf_4d_reshape(data);
                end
                save('Condition_effect_1.mat','Condition_effect', '-v7.3');
                clear  Condition_effect_1
                save('Yhat.mat','Yhat', '-v7.3');                                          
                Res = data - Yhat;
                clear Yhat
                save('Res.mat','Res', '-v7.3');                                            
                clear Res
                LIMO.design.status = 'done';
                LIMO.design.method = 'Generalized Welch''s method';
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            else
                LIMO.design.method = 'IRLS';
                LIMO.design.status = 'to do';
                save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
                clear data LIMO
                LIMOPath = LIMO.dir; 
                limo_eeg(4); return
            end
        else
            disp('Analysis aborted')
            return
        end

        
        % do the bootsrap for the 1-way robust ANOVA  -- done in limo_eeg(4) for other designs
        % ------------------------------------------------------------------------------------
        if LIMO.design.bootstrap ~= 0 &&  ...
                LIMO.design.fullfactorial == 0 && LIMO.design.nb_continuous == 0
            mkdir(fullfile(LIMO.dir,'H0'));
            
            if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                data = limo_tf_4d_reshape(data);
            end
            
            for c=1:(size(LIMO.design.X,2)-1) % size(LIMO.data.data,2) % center data
                index = find(LIMO.design.X(:,c));
                data(:,:,index) = data(:,:,index) - repmat(limo_trimmed_mean(data(:,:,index)),[1 1 length(index)]);
            end
            boot_table            = limo_create_boot_table(data,LIMO.design.bootstrap);
            H0_Condition_effect   = NaN(size(data,1),size(data,2),2,LIMO.design.bootstrap);
            
            array = find(~isnan(data(:,1,1))); % skip empty channels
            for b=1:LIMO.design.bootstrap
                fprintf('computing boostrap %g/%g\n',b,LIMO.design.bootstrap);
                for channel=1:size(array,1)
                    e     = array(channel);
                    index = find(~isnan(squeeze(data(e,1,:))));
                    X     = LIMO.design.X(index,1:end-1);
                    if sum(sum(X) == 0) ==0
                        [H0_Condition_effect(e,:,1,b), H0_Condition_effect(e,:,2,b)] = limo_robust_1way_anova(squeeze(data(e,:,boot_table{e}(:,b))),X,20); % no intercept in this model
                    end
                end
            end
            
            if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
                H0_Condition_effect = limo_tf_5d_reshape(H0_Condition_effect);
            end
            save(['H0' filesep 'H0_Condition_effect_1'],'H0_Condition_effect', '-v7.3');
            save(['H0' filesep 'boot_table'],'boot_table', '-v7.3');
            clear data H0_Condition_effect ;
        end
        
        % do TFCE
        if LIMO.design.tfce ~= 0 % check if tfce is on and if more than one channel
            fprintf('Thresholding bootstrapped N-way ANOVA using TFCE \n');
            limo_tfce_handling(fullfile(LIMO.dir,'Condition_effect_1.mat'));
            LIMO.design.tfce = 1;
            save(fullfile(LIMO.dir,'LIMO.mat'));
        end

        if strcmp(LIMO.design.method,'Generalized Welch''s method')
            disp('Robust Gp ANOVA for trimmed means (generalization of Welch''s method) done')
        else
            disp('ANOVA/ANCOVA analysis done');
        end
        LIMOPath = LIMO.dir;
        
        
        %----------------------------------------------------------------------------------------------
        % Repeated Measure ANOVA (multivariate approach) - bootstrap centering data
        %----------------------------------------------------------------------------------------------
    case {6}
        
        if ischar(varargin{2})
            data = load(varargin{2});
            data = data.(cell2mat((fieldnames(data))));
        else
            data = varargin{2}; % e,f,subjects,measures
        end
        
        gp_vector         = varargin{3}; % length of data, indices groups
        factor_levels     = varargin{4}; % vector eg [2 3]
        if ischar(varargin{5})
            LIMO = load(varargin{5});
            LIMO = LIMO.LIMO;
        elseif isstruct(varargin{5})
            LIMO = varargin{5};
        end
        
        if strcmp(LIMO.Analysis,'Time-Frequency') || strcmp(LIMO.Analysis,'ITC')
            tmp = NaN(size(data,1), size(data,2)*size(data,3),size(data,4),size(data,5));
            for measure = 1:size(data,5)
                if size(data,1) == 1
                    tmp(:,:,:,measure) = limo_tf_4d_reshape(data(1,:,:,:,measure));
                else
                    tmp(:,:,:,measure) = limo_tf_4d_reshape(squeeze(data(:,:,:,:,measure)));
                end
            end
            clear data; data=tmp; clear tmp;
            if size(data,3) ~= length(gp_vector)
                error('gp vector is not commensurate to the data (dimension 4)')
            end
            if size(data,4) ~= prod(factor_levels)
                error('factor_levels are not commensurate to the data (dimension 5)')
            end
        else
            if size(data,3) ~= length(gp_vector)
                error('gp vector is not commensurate to the data (dimension 3)')
            end
            if size(data,4) ~= prod(factor_levels)
                error('factor_levels are not commensurate to the data (dimension 4)')
            end
        end
        
        if nargin > 5
            if strcmpi(varargin{6},'go')
                go = varargin{7};
            end
        else
            go = 'no';
        end
        clear varargin
        
        % ------------------------------------------------
        % update the LIMO structure
        if isfield(LIMO,'dir')
            cd(LIMO.dir)
        end
        LIMO.data.Cat                = gp_vector;
        LIMO.data.Cont               = 0;
        LIMO.data.data_dir           = pwd;
        LIMO.design.type_of_analysis = 'Mass-univariate';
        LIMO.design.nb_conditions    = length(unique(gp_vector));
        LIMO.design.nb_interactions  = 0;
        LIMO.design.nb_continuous    = 0;
        LIMO.design.fullfactorial    = 0;
        LIMO.design.zscore           = 0;
        LIMO.design.repeated_measure = factor_levels;
        save('LIMO.mat','LIMO')
        
        % specific stuff for repeated measures
        % from the input we know which case to handle
        if unique(gp_vector) == 1
            % one sample
            if length(factor_levels) ==1
                type = 1; % simple repeated measure ANOVA (1 effect)
            elseif length(factor_levels) >1
                type = 2; % multiple factors (effects and interactions)
            end
        else
            % k samples
            if length(factor_levels) ==1
                type = 3; % gp * repeated measure (2 effects + interactions)
            elseif length(factor_levels) >1
                type = 4; % gp * multiple factors (effects and interactions)
            end
        end
        
        
        % make files to be stored
        % -----------------------
        if type == 1 % one factor
            LIMO.design.method     = 'Mean'; % change to Trimmed Mean for robust ANOVA
            C                      = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            tmp_Rep_ANOVA          = NaN(size(data,1),size(data,2),1,2); % store F and p
            LIMO.design.effects{1} = 'Main effect';
            LIMO.design.C{1}       = C;
            x                      = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X          = [x ones(size(x,1),1)];
            
        elseif type == 2 % many factors
            LIMO.design.method = 'Mean';
            C                  = limo_OrthogContrasts(factor_levels);
            tmp_Rep_ANOVA      = NaN(size(data,1),size(data,2),length(C),2); % store F and p for each within factor and interactions
            LIMO.design.C      = C;
            index = length(factor_levels)+1;
            
            for i= 1:length(factor_levels)
                LIMO.design.effects{i} = ['Main effect ' num2str(i)];
            end
            
            for i= 2:length(factor_levels)
                n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1)
                    LIMO.design.effects{index} = ['Interaction ' num2str(n(j,:))]; index = index+1;
                end
            end
            x             = kron(eye(prod(factor_levels)),ones(size(data,3),1));
            LIMO.design.X = [x ones(size(x,1),1)];
            
        elseif type == 3 % one factor within and one factor between
            LIMO.design.method                = 'Mean';
            gp_values                         = unique(gp_vector);
            k                                 = length(gp_values);
            X                                 = NaN(size(gp_vector,1),k+1);
            for g =1:k
                X(:,g) = gp_vector == gp_values(g);
            end
            X(:,end) = 1; % design matrix for gp effects
            C                                 = [eye(size(data,4)-1) ones(size(data,4)-1,1).*-1]; % contrast
            tmp_Rep_ANOVA                     = NaN(size(data,1),size(data,2),1,2);
            Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),1,2);
            LIMO.design.C{1}                  = C;
            LIMO.design.effects{1}            = 'Main effect';
            x                                 = kron(X(:,1:k),eye(prod(factor_levels)));
            LIMO.design.X                     = [x sum(x,2)]; % just for display
            LIMO.design.nb_interactions       = length(LIMO.design.C);
           
        elseif type == 4 % many factors within and one factor between
            LIMO.design.method                = 'Mean';
            gp_values                         = unique(gp_vector);
            k                                 = length(gp_values);
            X                                 = NaN(size(gp_vector,1),k+1);
            for g =1:k
                X(:,g) = gp_vector == gp_values(g);
            end
            X(:,end)                          = 1; % design matrix for gp effects
            C                                 = limo_OrthogContrasts(factor_levels);
            tmp_Rep_ANOVA                     = NaN(size(data,1),size(data,2),length(C),2);
            Rep_ANOVA_Gp_effect               = NaN(size(data,1),size(data,2),2);
            tmp_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),length(C),2);
            LIMO.design.C                     = C;
            for i= 1:length(factor_levels)
                LIMO.design.effects{i} = ['Main effect ' num2str(i)];
            end
            index                             = length(factor_levels)+1;
            for i= 2:length(factor_levels)
                n = nchoosek([1:length(factor_levels)],i);
                for j=1:size(n,1)
                    LIMO.design.effects{index} = ['Interaction ' num2str(n(j,:))];
                    index = index+1;
                end
            end
            x                                  = kron(X(:,1:k),eye(prod(factor_levels)));
            LIMO.design.X                      = [x sum(x,2)]; % just for display
            LIMO.design.nb_interactions        = length(LIMO.design.C);
       end
        
        if isempty(dir('Rep_ANOVA_Factor*.mat'))
            
            % check the design with user
            % --------------------------
            if ~strcmpi(go,'Yes')
                figure('Name','LIMO design'); set(gcf,'Color','w');
                imagesc(LIMO.design.X); colormap('gray');
                title('ANOVA model','FontSize',16);xlabel('regressors');
                ylabel('subjects'); drawnow;
                go = questdlg('start the analysis?');
                close('LIMO design')
                if ~strcmpi(go,'Yes')
                    return
                end
            end
            save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO');
            
            % do the analysis
            % ---------------
            array = find(~isnan(data(:,1,1,1)));
            for e = 1:length(array)
                channel = array(e);
                fprintf('analyse channel %g/%g\n ...', channel,size(data,1));
                tmp = squeeze(data(channel,:,:,:));
                if size(data,2) == 1
                    Y = ones(1,size(tmp,1),size(tmp,2)); Y(1,:,:) = tmp;
                    gp = gp_vector(find(~isnan(Y(1,:,1))),:);
                    Y = Y(:,find(~isnan(Y(1,:,1))),:);
                else
                    Y = tmp(:,find(~isnan(tmp(1,:,1))),:);
                    gp = gp_vector(find(~isnan(tmp(1,:,1))),:);
                end
                
                if type == 3 || type == 4
                    XB = X(find(~isnan(tmp(1,:,1))),:);
                end
                
                if type == 1
                    if strcmp(LIMO.design.method,'Trimmed Mean')
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,1,1) = result.F;
                    tmp_Rep_ANOVA(channel,:,1,2) = result.p;
                elseif type == 2
                    if strcmp(LIMO.design.method,'Trimmed Mean')
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,:,1) = result.F';
                    tmp_Rep_ANOVA(channel,:,:,2) = result.p';
                elseif type == 3
                    if strcmp(LIMO.design.method,'Trimmed Mean')
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C,XB); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,1,1)                   = result.repeated_measure.F;
                    tmp_Rep_ANOVA(channel,:,1,2)                   = result.repeated_measure.p;
                    Rep_ANOVA_Gp_effect(channel,:,1)               = result.gp.F;
                    Rep_ANOVA_Gp_effect(channel,:,2)               = result.gp.p;
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,1) = result.interaction.F;
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,2) = result.interaction.p;
                elseif type == 4
                    if strcmp(LIMO.design.method,'Trimmed Mean')
                        result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB); % trimmed means
                    else
                        result = limo_rep_anova(Y,gp,factor_levels,C,XB); % usual means
                    end
                    tmp_Rep_ANOVA(channel,:,:,1)                     = result.repeated_measure.F';
                    tmp_Rep_ANOVA(channel,:,:,2)                     = result.repeated_measure.p';
                    Rep_ANOVA_Gp_effect(channel,:,1)                 = result.gp.F;
                    Rep_ANOVA_Gp_effect(channel,:,2)                 = result.gp.p;
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,:,1) = result.interaction.F';
                    tmp_Rep_ANOVA_Interaction_with_gp(channel,:,:,2) = result.interaction.p';
                end
                nb_effects = size(tmp_Rep_ANOVA,3);
                clear tmp Y gp result
            end
            
            % save stuff
            % ---------
            Rep_filenames = cell(1,nb_effects);
            for i=1:nb_effects
                if contains(LIMO.design.effects{i},'Main effect')
                    if isfield(LIMO.design,'factor_names')
                        Rep_filenames{i} = sprintf('Rep_ANOVA_Main_effect_%g_%s.mat',i,LIMO.design.factor_names{i});
                    else
                        Rep_filenames{i} = sprintf('Rep_ANOVA_Main_effect_%g.mat',i);
                    end
                else
                    Interaction = LIMO.design.effects{i}(length('Interaction')+1:end);
                    Interaction(isspace(Interaction)) = [];
                    Rep_filenames{i} = sprintf('Rep_ANOVA_Interaction_Factors_%s.mat',Interaction);
                end
                
                % save each factor effect as F/p values
                % use reshape instead of squeeze in case there is only 1 channel
                Rep_ANOVA = reshape(tmp_Rep_ANOVA(:,:,i,:),...
                    [size(tmp_Rep_ANOVA,1) size(tmp_Rep_ANOVA,2) size(tmp_Rep_ANOVA,4)]);
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    Rep_ANOVA = limo_tf_4d_reshape(Rep_ANOVA);
                end
                save(Rep_filenames{i},'Rep_ANOVA', '-v7.3');
                if nargout ~= 0, LIMOPath = [fullfile(pwd,Rep_filenames{i}),'.mat']; end
            end
            
            if type == 3 || type ==4
                IRep_filenames = cell(1,nb_effects);
                for i=1:size(tmp_Rep_ANOVA_Interaction_with_gp,3)
                    if contains(LIMO.design.effects{i},'Main effect')
                        if isfield(LIMO.design,'factor_names')
                            IRep_filenames{i} = sprintf('Rep_ANOVA_Interaction_gp_Factor_%g_%s.mat',i,LIMO.design.factor_names{i});
                        else
                            IRep_filenames{i} = sprintf('Rep_ANOVA_Interaction_gp_Factor_%g.mat',i);
                        end
                    else
                        Interaction = LIMO.design.effects{i}(length('Interaction')+1:end);
                        Interaction(isspace(Interaction)) = [];
                        IRep_filenames{i} = sprintf('Rep_ANOVA_Interaction_gp_Factors_%s.mat',Interaction);
                    end
                    
                    % save each interaction effect as F/p values
                    Rep_ANOVA_Interaction_with_gp = reshape(tmp_Rep_ANOVA_Interaction_with_gp(:,:,i,:),...
                        [size(tmp_Rep_ANOVA_Interaction_with_gp,1) size(tmp_Rep_ANOVA_Interaction_with_gp,2) size(tmp_Rep_ANOVA_Interaction_with_gp,4)]);
                    if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                        Rep_ANOVA_Interaction_with_gp = limo_tf_4d_reshape(Rep_ANOVA_Interaction_with_gp);
                    end
                    save(IRep_filenames{i},'Rep_ANOVA_Interaction_with_gp', '-v7.3');
                    clear Rep_ANOVA_Interaction_with_gp;
                    if nargout ~= 0, LIMOPath = [fullfile(pwd,IRep_filenames{i}),'.mat']; end
                end
                
                % Main group effet
                save('Rep_ANOVA_Gp_effect.mat','Rep_ANOVA_Gp_effect','-v7.3');
                if nargout ~= 0, LIMOPath = fullfile(pwd,'Rep_ANOVA_Gp_effect.mat'); end
            end
        end
        
        % if skipping the above
        if ~exist('nb_effects','var')
            nb_effects = length(LIMO.design.C);
        end
        
        % clear up all tmp files
        clear tmp_Rep_ANOVA
        if type == 3 || type == 4
            clear Rep_ANOVA_Gp_effect tmp_Rep_ANOVA_Interaction_with_gp
        end
        
        % ----------------------------------------------------------------
        % now do the bootstrap
        % ---------------------
        if LIMO.design.bootstrap > 0
            boot_files = dir(fullfile(LIMO.dir,'H0'));
            if ~isempty(boot_files)
                answer = questdlg('a boostrap file already exist - overwrite?','data check','Yes','No','Yes');
                if strcmp(answer,'No')
                    return
                end
            end
            
            % create files to store bootstrap under H1 and H0
            mkdir(fullfile(LIMO.dir,'H0'))
            disp('making bootstrap files ...')
            if type ==1
                tmp_boot_H0_Rep_ANOVA = NaN(size(data,1),size(data,2),1,2,LIMO.design.bootstrap);
                X = [];
            elseif type == 2
                tmp_boot_H0_Rep_ANOVA = NaN(size(data,1),size(data,2),length(C),2,LIMO.design.bootstrap);
                X = [];
            elseif type == 3
                tmp_boot_H0_Rep_ANOVA = NaN(size(data,1),size(data,2),1,2,LIMO.design.bootstrap);
                H0_Rep_ANOVA_Gp_effect = NaN(size(data,1),size(data,2),2,LIMO.design.bootstrap);
                tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),1,2,LIMO.design.bootstrap);
            else
                tmp_boot_H0_Rep_ANOVA = NaN(size(data,1),size(data,2),length(C),2,LIMO.design.bootstrap);
                H0_Rep_ANOVA_Gp_effect = NaN(size(data,1),size(data,2),2,LIMO.design.bootstrap);
                tmp_boot_H0_Rep_ANOVA_Interaction_with_gp = NaN(size(data,1),size(data,2),length(C),2,LIMO.design.bootstrap);
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
            save(fullfile(LIMO.dir,['H0', filesep, 'centered_data']), 'centered_data', '-v7.3');
            
            % create an index to use across all channels and frames
            % (different per gp but identical across conditions)
            disp('making random table...')
            if LIMO.design.bootstrap == 1; LIMO.design.bootstrap = 1000; end
            boot_table = limo_create_boot_table(squeeze(data(:,:,:,1)),LIMO.design.bootstrap);
            save(fullfile(LIMO.dir,['H0', filesep, 'boot_table']), 'boot_table', '-v7.3');
            
            % compute bootstrap under H0 for F and p
            fprintf('Bootstrapping Repeated Measures ANOVA\n');
            parfor B=1:LIMO.design.bootstrap
                array = find(~isnan(data(:,1,1,1)));
                
                % preallocation for parfor
                if type ==1
                    tmp_boot_H0_Rep_ANOVA_sub = NaN(size(data,1),size(data,2),1,2);
                elseif type == 2
                    tmp_boot_H0_Rep_ANOVA_sub = NaN(size(data,1),size(data,2),length(C),2);
                elseif type == 3
                    tmp_boot_H0_Rep_ANOVA_sub = NaN(size(data,1),size(data,2),1,2);
                    H0_Rep_ANOVA_Gp_effect_sub = NaN(size(data,1),size(data,2),2);
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp_sub = NaN(size(data,1),size(data,2),1,2);
                else
                    tmp_boot_H0_Rep_ANOVA_sub = NaN(size(data,1),size(data,2),length(C),2);
                    H0_Rep_ANOVA_Gp_effect_sub = NaN(size(data,1),size(data,2),2);
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp_sub = NaN(size(data,1),size(data,2),length(C),2);
                end
 
                for e = 1:length(array)
                    channel = array(e);
                    if e == 1
                        fprintf('parallel boot %g channel %g',B,channel);
                    elseif e==length(array)
                        fprintf(' %g\n',channel);
                    else
                        fprintf(' %g',channel);
                    end
                    % get data per channel
                    tmp = squeeze(centered_data(channel,:,boot_table{channel}(:,B),:));
                    if size(centered_data,2) == 1
                        Y  = ones(1,size(tmp,1),size(tmp,2)); Y(1,:,:) = tmp;
                        gp = gp_vector(find(~isnan(Y(1,:,1))),:);
                        Y  = Y(:,find(~isnan(Y(1,:,1))),:);
                    else
                        Y  = tmp(:,find(~isnan(tmp(1,:,1))),:);
                        gp = gp_vector(find(~isnan(tmp(1,:,1))));
                    end
                    
                    if type == 3 || type == 4
                        XB = X(find(~isnan(tmp(1,:,1))));
                    else
                        XB = [];
                    end
                    
                    if type == 1
                        if strcmp(LIMO.design.method,'Trimmed Mean')
                            result = limo_robust_rep_anova(Y,gp,factor_levels,C);
                        else
                            result = limo_rep_anova(Y,gp,factor_levels,C);
                        end
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,1,1) = result.F;
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,1,2) = result.p;
                    elseif type == 2
                        if strcmp(LIMO.design.method,'Trimmed Mean')
                            result = limo_robust_rep_anova(Y,gp,factor_levels,C);
                        else
                            result = limo_rep_anova(Y,gp,factor_levels,C);
                        end
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,1,1) = result.F';
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,1,2) = result.p';
                    elseif type == 3
                        if strcmp(LIMO.design.method,'Trimmed Mean')
                            result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
                        else
                            result = limo_rep_anova(Y,gp,factor_levels,C,XB);
                        end
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,1,1) = result.repeated_measure.F;
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,1,2) = result.repeated_measure.p;
                        H0_Rep_ANOVA_Gp_effect_sub(channel,:,1) = result.gp.F;
                        H0_Rep_ANOVA_Gp_effect_sub(channel,:,2) = result.gp.p;
                        tmp_boot_H0_Rep_ANOVA_Interaction_with_gp_sub(channel,:,1,1) = result.interaction.F;
                        tmp_boot_H0_Rep_ANOVA_Interaction_with_gp_sub(channel,:,1,2) = result.interaction.p;
                    elseif type == 4
                        if strcmp(LIMO.design.method,'Trimmed Mean')
                            result = limo_robust_rep_anova(Y,gp,factor_levels,C,XB);
                        else
                            result = limo_rep_anova(Y,gp,factor_levels,C,XB);
                        end
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,:,1) = result.repeated_measure.F';
                        tmp_boot_H0_Rep_ANOVA_sub(channel,:,:,2) = result.repeated_measure.p';
                        H0_Rep_ANOVA_Gp_effect_sub(channel,:,1)  = result.gp.F;
                        H0_Rep_ANOVA_Gp_effect_sub(channel,:,2)  = result.gp.p;
                        tmp_boot_H0_Rep_ANOVA_Interaction_with_gp_sub(channel,:,:,1) = result.interaction.F';
                        tmp_boot_H0_Rep_ANOVA_Interaction_with_gp_sub(channel,:,:,2) = result.interaction.p';
                    end
                end
                tmp_boot_H0_Rep_ANOVA(:,:,:,:,B)  = tmp_boot_H0_Rep_ANOVA_sub;
                if type == 3 || type == 4
                    H0_Rep_ANOVA_Gp_effect(:,:,:,:,B) = H0_Rep_ANOVA_Gp_effect_sub;
                    tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,:,:,B) = tmp_boot_H0_Rep_ANOVA_Interaction_with_gp_sub;
                end
            end
            
            % save
            for i=1:size(tmp_boot_H0_Rep_ANOVA,3)
                name = sprintf('H0_%s',Rep_filenames{i});
                H0_Rep_ANOVA = NaN(size(tmp_boot_H0_Rep_ANOVA,1), size(tmp_boot_H0_Rep_ANOVA, 2), size(tmp_boot_H0_Rep_ANOVA, 4), size(tmp_boot_H0_Rep_ANOVA, 5));
                H0_Rep_ANOVA(:,:,:,:) = squeeze(tmp_boot_H0_Rep_ANOVA(:,:,i,:,:)); % save each factor effect as F/p/LIMO.design.bootstrap values
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    H0_Rep_ANOVA = limo_tf_5d_reshape(H0_Rep_ANOVA);
                end
                save(['H0', filesep, name],'H0_Rep_ANOVA', '-v7.3');
            end
            
            if type == 3 || type ==4
                for i=1:size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp,3)
                    name = sprintf('H0_%s',IRep_filenames{i});
                    H0_Rep_ANOVA_Interaction_with_gp = NaN(size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp,1), size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp, 2), size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp, 4), size(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp, 5));
                    H0_Rep_ANOVA_Interaction_with_gp(:,:,:,:) = squeeze(tmp_boot_H0_Rep_ANOVA_Interaction_with_gp(:,:,i,:,:)); % save each interaction effect as F/p values
                    if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                        H0_Rep_ANOVA_Interaction_with_gp = limo_tf_5d_reshape(H0_Rep_ANOVA_Interaction_with_gp);
                    end
                    save(['H0', filesep, name],'H0_Rep_ANOVA_Interaction_with_gp', '-v7.3'); clear H0_Rep_ANOVA_Interaction_with_gp;
                end
                
                if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                    H0_Rep_ANOVA_Gp_effect = limo_tf_5d_reshape(H0_Rep_ANOVA_Gp_effect);
                end
                save(['H0', filesep, 'H0_Rep_ANOVA_Gp_effect'], 'H0_Rep_ANOVA_Gp_effect', '-v7.3');
            end
        end
        LIMO.design.bootstrap = LIMO.design.bootstrap;
        save(fullfile(LIMO.dir,'LIMO.mat'));
        
        % ------------------------- TFCE ---------------
        if LIMO.design.tfce ~= 0 % check if tfce is on and if more than one channel
            fprintf('Thresholding bootstrapped Rep ANOVA using TFCE \n');
            for i=1:nb_effects
                limo_tfce_handling(fullfile(LIMO.dir,Rep_filenames{i}));
            end
            LIMO.design.tfce = 1;
            save(fullfile(LIMO.dir,'LIMO.mat'));
        end
        disp('Repeated Measures ANOVA done')
end



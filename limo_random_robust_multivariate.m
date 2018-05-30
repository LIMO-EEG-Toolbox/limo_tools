function filepath = limo_random_robust_multivariate(varargin)

% This function makes the result files for the random effects of various tests
% after running mulitvarite analyses on first level. i.e., tests on
% timevectors
% This function also organizes and makes files for boostrap. It is interfaced with
% limo_random_effect which itself interfaces with the user to select and pass
% the data in the appropriate format. limo_random_robust_multivariate calls low level
% functions like limo_trimci to perform the actual computation.
%
% FORMAT limo_random_robust_multivariate(test,data,nboot,tfce)
%        filepath = limo_random_robust_multivariate(test,data,parameter, nboot,tfce);
%
% INPUTS --- note that a LIMO.mat is also expected in the current directory
%            see limo_random_select_multivariate
%
% limo_random_robust_multivariate(1,y,parameter number,nboot,tfce)
%                    1 = a one-sample t-test
%                    y = data (dim  time/freq x subjects)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%                    nboot = nb of resamples (0 for none)
%                    tfce = 0/1 to compute tcfe (only if nboot ~=0).
%
% limo_random_robust_multivariate(2,y1,y2,parameter number,nboot,tfce);
%                    2 = two samples t-test
%                    y1 = data (dim  time or freq, subjects)
%                    y2 = data (dim time or freq, subjects)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%                    nboot = nb of resamples (0 for none)
%                    tfce = 0/1 to compute tcfe (only if nboot ~=0).
%
% limo_random_robust_multivariate(3,y1,y2,parameter number,nboot,tfce);
%                    3 = paired t-test
%                    y1 = data (dim time or freq, subjects)
%                    y2 = data (dim time or freq, subjects)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%                    nboot = nb of resamples (0 for none)
%                    tfce = 0/1 to compute tcfe (only if nboot ~=0).
%
% limo_random_robust_multivariate(4,y,X,parameter number,nboot,tfce);
%                    4 = regression analysis
%                    y = data (time or freq, subjects)
%                    X = continuous regressor(s)
%                    parameter number = describe which parameter is currently analysed (e.g. 1 - use for maming only)
%                    nboot = nb of resamples (0 for none)
%                    tfce = 0/1 to compute tcfe (only if nboot ~=0).
%
% limo_random_robust_multivariate(5,y,cat,cont,LIMO,nboot,tfce)
%                    5 = N-way ANOVA/ANCOVA
%                    y = data (time or freq, subjects)
%                    cat = categorical variable(s)
%                    cont = continuous regressors (covariates)
%                    LIMO the basic structure with data, design and channel info
%                    nboot = nb of resamples (0 for none)
%                    tfce = 0/1 to compute tcfe (only if nboot ~=0).
%
% limo_random_robust_multivariate(6,y,gp,factor_levels,LIMO,nboot,tfce)
%                    6 = Repeated measures ANOVA/ANCOVA using multivariate approach
%                    y = data (dim time or freq, subjects, measures)
%                    gp = a vector defining gps
%                    factor_levels = a vector specifying the levels of each repeated measure factor
%                    LIMO the basic structure with data, design and channel info
%                    nboot = nb of resamples (0 for none)
%                    tfce = 0/1 to compute tcfe (only if nboot ~=0).
%
% OUPUT
% write on the disk matrices correponding to the test (Yr and LIMO.mat are generated in limo_random_select_multivariate,
% and for Regression, ANOVA, the LIMO.mat structure is updated)
%
% 1 one_sample_parameter_X (frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_one_sample_ttest_parameter_X (frames, [T values under H0, p values under H0], nboot)
%
% 2 two_samples_parameter_X (frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_two_samples_ttest_parameter_X (frames, [T values under H0, p values under H0], nboot)
%
% 3 paired_samples_parameter_X (frames [time, freq or freq-time], [mean value, se, df, t, p])
%   H0_paired_samples_ttest_parameter_X (frames, [T values under H0, p values under H0], nboot)
%
% 4 R2 (frames [time, freq or freq-time], [F p values])
%   H0_R2 (frames, [F p], nboot)
%   Covariate_effect_X (frames [time, freq or freq-time], [F p values])
%   H0_Covariate_effect_X (frames, [F p], nboot)
%
% 5 Condition_effect_X (frames [time, freq or freq-time], [F p values])
%   H0_Condition_effect_X (frames, [F p], nboot)
%   Covariate_effect_X (frames [time, freq or freq-time], [F p values])
%   H0_Covariate_effect_X (frames, [F p], nboot)
%
% 6 Rep_ANOVA_Factor_X (frames [time, freq or freq-time], [F p values])
%   Rep_ANOVA_Gp_effect (frames [time, freq or freq-time], [F p values])
%   Rep_ANOVA_Interaction_gp_Factor_X (frames [time, freq or freq-time], [F p values])
%   H0_XXXXX same as above, including nboot on the last dimension
%
% filepath - Path to the contrast result file. Mainly for EEGALB functionality to
%            allow loading test directly.
%
% See also LIMO_TRIMCI LIMO_YUEN_TTEST LIMO_YUEND_TTEST LIMO_ROBUST_1WAY_ANOVA
% LIMO_GLM1 LIMO_EEG(4) LIMO_EEG_TF(4) LIMO_REP_ANOVA LIMO_CREATE_BOOT_TABLE
% -----------------------------
% Copyright (C) LIMO Team 2015

% v1: Cyril Pernet and Guillaume Rousselet 24-08-2009
% v2: Cyril Pernet 12-07-2010
% v3: GAR 25-08-2010 - fixed bug in one-sample bootstrap + added new NaN check of boot indices
% 29-08-2010: Cyril Pernet and Guillaume Rousselet ANOVAs sorted out (boot_index per cell)
% 11/12-2012 Marianne Latinus added TFCE computation. Also added a check
% that if a bootstrap file with same characteristics (nboot, and nelec) on
% the same data exists the bootstrap step is skipped.
% v4: 02/06/2013 cleaned up + changed ANOVA to create matrices amd run
% limo_glm(4) - updated repeated measure - boot_tables are made via
% function -- thx to Benedikt Ehinger for spotting a few bugs here and
% there
% Novembre 2013 - Matt Craddock fixes applied on rep ANOVA for single
% electrode analyses
% v5: May 2014 - added time-frequency analyses + added check for matrices with
% NaNs everywhere (ie empty channel) + changed test of hypotheses to
% trimmed means when possible
% v6: May 2018 - added possibility to perform analyses on time x subjects
% vector (for example, coming from multivariate analyses). 


%% inputs checks
filepath = pwd;
type = varargin{1};
if type == 1
    tfce = varargin{5};
else
    tfce = varargin{6};
end
alpha = .05; % used for a basic computation like in t-test but data aren't thresholded here

% data check
data = varargin{2};
if ndims(data) ~=2
    errordlg('data should have 2 dimensions');
end
%% start

switch type
    %--------------------------------------------------------------------------
    % One Sample t-test // bootstrap-t method
    %--------------------------------------------------------------------------
    case {1}
        
        load LIMO
        parameter = varargin{3};
        nboot     = varargin{4};
        tfce      = varargin{5};
        clear varargin
        
        one_sample = NaN(size(data,1), 5);
        name = sprintf('one_sample_ttest_parameter_%g',parameter);
        Y = data;
        [one_sample(:,4),one_sample(:,1),trimci,one_sample(:,2),one_sample(:,5),tcrit,one_sample(:,3)]=limo_trimci(Y,5, 0.05, 1/LIMO.nb_conditions_fl);

        save ([name],'one_sample', '-v7.3')
        if nargout ~= 0, filepath = [fullfile(pwd,[name]),'.mat']; end            

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
                H0_one_sample = NaN(size(data,1),2,nboot); % stores T and p values for each boot under H0
                % subtract chance level prediction
                data = data - 1/LIMO.nb_conditions_fl;
                % create centered data to estimate H0
                centered_data = data - squeeze(repmat(limo_trimmed_mean(data),[1 1 size(data,2)]));
                % centered_data = data - repmat(nanmean(data,3),[1 1 size(data,3)]);
                % get boot table
                disp('making boot table ...')
                boot_table = randi(size(data,2), size(data,2), nboot);
                save(['H0', filesep, 'boot_table'], 'boot_table')

                % get results under H0
                Y = centered_data; 
                t = NaN(size(Y,1), nboot);
                p = NaN(size(Y,1), nboot);
%                 if exist('parfor','file') ~=0
%                     parfor b=1:nboot
%                         [t(:,b),~,~,~,p(:,b),~,~]=limo_trimci(Y(:,boot_table(:,b)), 5);
%                     end
% 
%                     for b=1:nboot
%                         H0_one_sample(:,1,b) = t(:,b);
%                         H0_one_sample(:,2,b) = p(:,b);
%                     end
%                 else
                    for b=1:nboot
                        [H0_one_sample(:,1,b),tmdata,trimci,se,H0_one_sample(:,2,b),tcrit,df]=limo_trimci(Y(:,boot_table(:,b)), 5);
                        fprintf('boot %d\n',b)
                        % [m,dfe,ci,sd,n,H0_one_sample(electrode,:,1,b),H0_one_sample(electrode,:,2,b)] = limo_ttest(1,Y(1,:,boot_table{electrode}(:,b)),0,5/100);
                    end
%                 end
                clear tmp Y
           

            save (['H0', filesep, boot_name],'H0_one_sample','-v7.3');
            end

            if tfce ~= 0
            end
        end % closes if nboot > 0
        disp('one sample t-test done')            
        
        
 
        %--------------------------------------------------------------------------
        % Two Samples t-test // percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {2}
        
        %--------------------------------------------------------------------------
        % Paired t-test // percentile bootstrap technique
        %--------------------------------------------------------------------------
    case {3}
        
       
        %------------------------------------------------------------------
        % Regression // percentile bootstrap under H0
        %------------------------------------------------------------------
    case {4}
         
        
        %----------------------------------------------------------------------------------------------
        % Repeated Measure ANOVA (multivariate approach) - bootstrap centering data
        %----------------------------------------------------------------------------------------------
    case {6}
        

end



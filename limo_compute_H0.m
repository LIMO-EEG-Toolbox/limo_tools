function varargout=limo_compute_H0(varargin)

% this function allows comnputing the null distribution of any of the
% LIMO EEG 2nd level statistical tests using bootstrap
%
% FORMAT [H0_data,boot_table] = limo_compute_H0(0,data,nboot,test);
%        status  = limo_compute_H0(type,data,label,nboot)
%
% -------------------------------------------------------------------------
% [H0_data,boot_table] = limo_compute_H0(0,data,nboot,test)
% the generic case to compute H0 using limo stat functions
%
% INPUTS: 0 stands for type = 0 and indicates the generic case
%        data this is a 4D or 3D matrix of data
%        nboot is the number of boostrap to use (recommended at least 800)
%        test is one of limo statistical test: 'limo_trimci' 'limo_yuen'
%                                              'limo_yuend' 'limo_glm'
%                                              'limo_'
%
% OUTPUTS: H0_data is the matrix of statistical values (t/F and p) needed
%          to do a correction for multiple comparisons
%          boot_table is the resampling table used
%
% -------------------------------------------------------------------------
% status = limo_compute_H0(type,data,label,boot_table)
% this is the format used by LIMO EEG, and in this case it creates data on
% the drive using dedicated names and format - see limo_random_robust for
% details of inputs and outputs
%
% -------------------------------------------------------------------------
% see also limo_create_boot_table
%
%

type = varargin(1);

% get the data and boot_table
if type == 0
    test = varargin(end);
    if strcmpi(test,'limo_yuen_ttest') || strcmpi(test,'limo_yuend_ttest')
        data1 = varagin(2);
        data2 = varagin(3);
        nboot = varargin(4);
    else
        data = varagin(2);
        nboot = varargin(3);
    end
elseif type == 1 || type > 3
    data = varagin(2);
    nboot = varargin(3);
else
    data1 = varagin(2);
    data2 = varagin(3);
    nboot = varargin(4);
end
clear varargin

% -------------------------------------------------------------------------

switch type
    case {1}
        % one sample t-test
        % -----------------
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
            centered_data = data - repmat(limo_trimmed_mean(data),[1 1 size(data,3)]);
            % centered_data = data - repmat(nanmean(data,3),[1 1 size(data,3)]);
            % get boot table
            disp('making boot table ...')
            boot_table = limo_create_boot_table(data,nboot);
            save(['H0', filesep, 'boot_table'], 'boot_table')
            
            % get results under H0
            for electrode = 1:size(data,1)
                fprintf('bootstrap: electrode %g parameter %g \n',electrode,parameter);
                tmp = centered_data(electrode,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                if exist('parfor','file') ~=0
                    parfor b=1:nboot
                        [t{b},~,~,~,p{b},~,~]=limo_trimci(Y(1,:,boot_table{electrode}(:,b)));
                    end
                    
                    for b=1:nboot
                        H0_one_sample(electrode,:,1,b) = t{b};
                        H0_one_sample(electrode,:,2,b) = p{b};
                    end
                else
                    for b=1:nboot
                        [H0_one_sample(electrode,:,1,b),tmdata,trimci,se,H0_one_sample(electrode,:,2,b),tcrit,df]=limo_trimci(Y(1,:,boot_table{electrode}(:,b)));
                        % [m,dfe,ci,sd,n,H0_one_sample(electrode,:,1,b),H0_one_sample(electrode,:,2,b)] = limo_ttest(1,Y(1,:,boot_table{electrode}(:,b)),0,5/100);
                    end
                end
                clear tmp Y
            end % closes for electrode
            
            if strcmp(LIMO.Analysis,'Time-Frequency') ||  strcmp(LIMO.Analysis,'ITC')
                H0_one_sample = limo_tf_5d_reshape(H0_one_sample);
            end
            save (['H0', filesep, boot_name],'H0_one_sample','-v7.3');
        end
        
        % ----------------------------------------------------------------------
        % ----------------------------------------------------------------------
        %                    THIS IS TYPE 0 - GENRERIC CASE
        % ----------------------------------------------------------------------
        % ----------------------------------------------------------------------
    case {0}
        % one sample t-test
        % -----------------
        if strcmpi(test,'limo_trimci')
            
            reshape = 0;
            if numel(size(data)) == 4
                reshape = 1; [elect,freq,time,obs]=size(data);
                data = limo_tf_4d_reshape(data,[elect,freq*time,obs]);
            end
            
            H0 = NaN(size(data,1), size(data,2),2,nboot); % stores T and p values for each boot under H0
            centered_data = data - repmat(limo_trimmed_mean(data),[1 1 size(data,3)]);
            disp('making boot table ...'); boot_table = limo_create_boot_table(data,nboot);
            
            for electrode = 1:size(data,1)
                fprintf('bootstrap: electrode %g parameter %g \n',electrode,parameter);
                tmp = centered_data(electrode,:,:); Y = tmp(1,:,find(~isnan(tmp(1,1,:))));
                parfor b=1:nboot
                    [t{b},~,~,~,p{b},~,~]=limo_trimci(Y(1,:,boot_table{electrode}(:,b)));
                end
                
                for b=1:nboot
                    H0(electrode,:,1,b) = t{b};
                    H0(electrode,:,2,b) = p{b};
                end
                clear tmp Y
            end % closes for electrode
            
            if reshape == 1
                H0 = limo_tf_5d_reshape(H0,[elect,freq,time,nboot]);
            end
        end
        clear data; varargout{2} = boot_table; cleaer boot_table
        varargout{1} = H0; clear H0;
end



function limo_itc(varargin)

% Function to batch process Inter-Trial Coherence (ITC) from multiple
% conditions in limo. Partly adapted from limo_batch.
%
%
% FORMAT limo_itc
%        limo_itc(option,model,contrast)
%
% INPUT if empty uses GUI
%       option should be 'model specification' 'contrast only' or 'both'
%       model is a structure that specifiy information to build a model
%             model.set_files: a cell array of EEG.set (full path) for the different subjects
%             model.cat_files: a cell array of categorial variable files
%             model.cont_files: a cell array of continuous variable files
%             model.defaults: specifiy the parameters to use for each subject
%             model.defaults.analysis 'Time' 'Frequency' or 'Time-Frequency'
%             model.defaults.fullfactorial     0/1
%             model.defaults.zscore            0/1
%             model.defaults.start             starting time in ms
%             model.defaults.end               ending time in ms
%             model.defaults.lowf              starting point in Hz
%             model.defaults.highf             ending point in Hz
%             model.defaults.bootstrap         0/1
%             model.defaults.tfce              0/1
%             model.defaults.channloc          common channel locations (necessary if bootstrap = 1)
%      contrast is a structure that specify which contrasts to run for which subject
%             contrast.LIMO_files: a list of LIMO.mat (full path) for the different subjects
%                                  this is optional if option 'both' is selected
%             contrast.mat: a matrix of contrasts to run (assumes the same for all subjects)
%
% OUTPUT none - generate a directory per subject with GLM results in it
%
% see also limo_eeg limo_batch limo_import_t limo_import_f limo_import_tf and psom in external folder
%
%
% Andrew Stewart, May 2014 - adapted from limo_batch
% -----------------------------
% Copyright (C) LIMO Team 2014

% psom stuff see mode for parallel computing
opt.mode = 'session'; % run one after the other in the current matlab session
opt.flag_pause = false;

global EEGLIMO
global LIMO

%% what to do




itc_prompt = sprintf('LIMO ITC allows loading of Inter-Trial Coherence (ITC) data and use of LIMO tools on that data. \n\n  It requires that that data has previously been generated and saved in EEG.etc.itc for each subject.');
%helpdlg(itc_prompt, 'Using LIMO ITC')

option='model specification';



% ITC gui - loading dataset list and trim info
[model.set_files,model.cat_files,model.cont_files,model.defaults] = limo_itc_gui;
model.info = 'Running ITC data loading';

for subject = 1:size(model.set_files,1)
    
    subpath = fileparts(model.set_files{subject});
    cd(subpath);
    
    limo_itc_import_data(model.set_files{subject},model.cat_files,model.cont_files,model.defaults);
    load('ITC_analysis/LIMO.mat')
    LIMO.Analysis = 'ITC';
    EEGLIMO=pop_loadset(LIMO.data.data);
    cd 'ITC_analysis'
    
    
    disp('loading ITC data...');
    
    % Let's treat 2-condition ITC data like TF data with 2 trials
    Yitc = EEGLIMO.etc.itc(:,LIMO.data.trim_low_f:LIMO.data.trim_high_f,LIMO.data.trim1:LIMO.data.trim2);
    
    
    if size(Yitc,1)/EEGLIMO.nbchan == 2  % If double electrode count in ITC data, check for 2 now
        disp('*** - Found double electrode count in ITC data - taking as two conditions')
        model.Ncond{subject} = 2;
        %model.info = strcat(model.info,' x2 electrode count in ITC data, taking as two conditions');
        
        if model.Ncond{subject} ~= model.Ncond{1}
            subject; error('Some subjects appear to have differing numbers of conditions')
        end
        

        Y = nan(EEGLIMO.nbchan,size(Yitc,2),size(Yitc,3),2);
        Y(:,:,:,1) = abs(Yitc(1:EEGLIMO.nbchan,:,:));
        Y(:,:,:,2) = abs(Yitc(EEGLIMO.nbchan+1:end,:,:));
              
        % Load a 4D Y with ITC data
        LIMO.data.size4D= size(Y);
        LIMO.data.size3D= [LIMO.data.size4D(1) LIMO.data.size4D(2)*LIMO.data.size4D(3) LIMO.data.size4D(4)];
        
        
    elseif size(Yitc,1)/EEGLIMO.nbchan == 3  % If double electrode count in ITC data, check for 3 now
        disp('*** - Found triple electrode count in ITC data - taking as three conditions')
        model.Ncond{subject} = 3;
        %model.info = strcat(model.info,' x3 electrode count in ITC data, taking as three conditions');
        
        if model.Ncond{subject} ~= model.Ncond{1}
            subject; error('Some subjects appear to have differing numbers of conditions')
        end
        

        Y = nan(EEGLIMO.nbchan,size(Yitc,2),size(Yitc,3),3);
        Y(:,:,:,1) = abs(Yitc(1:EEGLIMO.nbchan,:,:));
        Y(:,:,:,2) = abs(Yitc(EEGLIMO.nbchan+1:EEGLIMO.nbchan*2,:,:));
        Y(:,:,:,3) = abs(Yitc(EEGLIMO.nbchan*2+1:EEGLIMO.nbchan*3,:,:));
        
        
              
        % Load a 4D Y with ITC data
        LIMO.data.size4D= size(Y);
        LIMO.data.size3D= [LIMO.data.size4D(1) LIMO.data.size4D(2)*LIMO.data.size4D(3) LIMO.data.size4D(4)];
        
        
    elseif size(Yitc,1)/EEGLIMO.nbchan == 1
        
        
         model.Ncond{subject} = 1;
        %model.info = strcat(model.info,' x3 electrode count in ITC data, taking as three conditions');
        
        if model.Ncond{subject} ~= model.Ncond{1}
            subject; error('Some subjects appear to have differing numbers of conditions')
        end
        model.info{subject} = '1 condition';
        Y = nan(EEGLIMO.nbchan,size(Yitc,2),size(Yitc,3),1);
        Y(:,:,:,1) = abs(Yitc(1:EEGLIMO.nbchan,:,:));
        
        
    else
        error('Check ITC data length')
    end
    
    LIMO.model = model;
    
    save Y Y -v7.3;
    save LIMO LIMO -v7.3;
    model.itc_data{subject} = pwd;
    
    
    clear Yitc
    
    
    
end

disp('Finished loading')


% ITC gui 2 - get stat test selection
model.test_select = limo_itc_gui2;
disp(model.test_select)


current = pwd;
% mkdir('limo_batch_report')




% Have a check of set files and data being in agreement


%% Check test, run

original_LIMO_dir = model.itc_data{1};


if model.test_select{1} == 1
    %% 1-samp t
    
    % Build selected ITC data into correct format
    Nsub = length(model.set_files);
    Ncond = model.Ncond{1};
    total_elecs = length(model.defaults.chanlocs);
    
    
    
    
    disp(model.test_select{2})
    
    Ybig = nan(total_elecs,size(Y,2),size(Y,3),Nsub*Ncond);
    
    % Populate Ybig with ITC data already saved
   
    for sub = 1:Nsub
        
        cd(model.itc_data{sub})
        
        load LIMO
        load Y
        
        
        for elec = 1:size(Y,1)
            
            org_elec = LIMO.data.chanlocs(elec).urchan;  % Find original elec index
            
            j=0; % Additional cond count
            for cond = 1:Ncond
                Ybig(org_elec,:,:,sub+j) = Y(elec,:,:,cond);
                j=j+1;
            end
            
            % Leave the rest as nan
            
        end
    end
    
    % --- Check values
    if model.defaults.bootstrap == 1
        nboot = 1000;
    else
        nboot = [];
    end
    
    tfce = model.defaults.tfce;
    if tfce == 1 && isfield(LIMO.data,'neighbouring_matrix') == 0  % Check we have neighb matrix. If not, create it.
        EEGLIMO.chanlocs = model.defaults.chanlocs;
        neighbdis = inputdlg('What neighbourhood distance should be used for TFCE neighbourhood matrix? (Perhaps 0.37 for 128 electrode systems)','Enter neighb distance',1,{'0.37'});
        neighbdis = str2num(neighbdis{1});
        [tmpneighbs, LIMO.data.neighbouring_matrix] = limo_get_channeighbstructmat(EEGLIMO,neighbdis);
    end
    
    
    
    Analysis_type = 'ITC';
    parameters = 1;
    
    LIMO.data.cond_tested = 1;
    LIMO.data.chanlocs = model.defaults.chanlocs;
    save LIMO LIMO
    
    
    % Run stats on this Ybig
    cd(original_LIMO_dir);
    limo_random_robust(model.test_select{1},Ybig,parameters,nboot,tfce);
    
    
    % Plot
    plotnow = 0;
    if plotnow == 1
        load one_sample_ttest_parameter_1

        limo_display_results_tf(LIMO, one_sample(:,:,:,4),1,['ITC ', model.test_select{2}])
    end
    
    
    
    
    
    
    
elseif model.test_select{1} == 2
    %% 2-samp t-test
    
    
    
    disp(model.test_select{2})
    
    conds1 = 1;
    conds2 = 2;
    
    if model.Ncond{1} > 2
        cond_text{1} = sprintf('There are %d conditions. Which should be tested here? \n\nChoose first condition(s):', model.Ncond{1});
        cond_text{2} = 'Second condition(s) to test:';
        cond_tested = inputdlg(cond_text,'Which conditions should be tested?',1,{'1','2,3'});
        
        conds1 = str2num(cond_tested{1});
        conds2 = str2num(cond_tested{2});
    end
    
    
    
    % Build selected ITC data into correct format
    Nsub = length(model.set_files);
    Nconds = [length(conds1) length(conds2)];
    
    total_elecs = length(model.defaults.chanlocs);
    
    
    
    
 
    
    Y1 = nan(total_elecs,size(Y,2),size(Y,3),Nsub*Nconds(1));
    Y2 = nan(total_elecs,size(Y,2),size(Y,3),Nsub*Nconds(2));
    
   % Populate Y1 and Y2 with ITC data already saved
   

    for sub = 1:Nsub
        
        cd(model.itc_data{sub})
        
        load LIMO
        load Y
        
        
        for elec = 1:size(Y,1)
            
            org_elec = LIMO.data.chanlocs(elec).urchan;  % Find original elec index
            
            
            for cond = 1:Nconds(1)  % For each cond going into Y1
                Y1(org_elec,:,:,sub-1+cond) = Y(elec,:,:,conds1(cond));
            end
            
            
            for cond = 1:Nconds(2)  % For each cond going into Y2
                Y2(org_elec,:,:,sub-1+cond) = Y(elec,:,:,conds2(cond));
            end
            
            % Leave the rest as nan
            
        end
    end
    
    size(Y1)
    size(Y2)
    
    Y1nans = mean(isnan(Y1(:)))
    Y2nans = mean(isnan(Y2(:)))

    
    % --- Check values
    if model.defaults.bootstrap == 1
        nboot = 1000;
    else
        nboot = 0;
    end
    
    tfce = model.defaults.tfce;
    if tfce == 1 && isfield(LIMO.data,'neighbouring_matrix') == 0  % Check we have neighb matrix. If not, create it.
        EEGLIMO.chanlocs = model.defaults.chanlocs;
        neighbdis = inputdlg('What neighbourhood distance should be used for TFCE neighbourhood matrix? (Perhaps 0.37 for 128 electrode systems)','Enter neighb distance',1,{'0.37'});
        neighbdis = str2num(neighbdis{1});
        [tmpneighbs, LIMO.data.neighbouring_matrix] = limo_get_channeighbstructmat(EEGLIMO,neighbdis);
    end
    
    
    
    
    Analysis_type = 'ITC';
    parameters = 1;
    
    
    
    
    % Run stats on this Ybig
    cd(original_LIMO_dir); 
    tag = sprintf('%s_conds%s_%s',model.test_select{2},cond_tested{1},cond_tested{2});
    tag(~ismember(tag,['A':'Z' 'a':'z' '1':'9' '_']))= ''; % Clean up tag string
    mkdir(tag);cd(tag)
    LIMO.data.cond_tested = cond_tested;
    LIMO.data.chanlocs = model.defaults.chanlocs;
    save LIMO LIMO
    
    limo_random_robust(model.test_select{1},Y1,Y2,parameters,nboot,tfce);
    disp([model.test_select{2} ' done'])
    
    % Plot
    plotnow = 0;
    if plotnow == 1
        
        
        if nboot == 0
            load two_samples_ttest_parameter_1
            limo_display_results_tf(LIMO, two_samples(:,:,:,4),1,['ITC F ', model.test_select{2}])
        else
            load /H0/H0_two_samples_ttest_parameter_1
            limo_display_results_tf(LIMO, H0_two_samples(:,:,:,4),1,['ITC F ', model.test_select{2}])
        end
            
        
    end
    
    
elseif model.test_select{1} == 3  % Paired t
    
    disp(model.test_select{2})
    disp('Not yet implemented')
    
    
elseif model.test_select{1} == 4 % Reg
    
    
    cont_data = model.cont_files;
    
    disp(model.test_select{2})
    
    conds1 = 1;
    conds2 = 2;
    
    if model.Ncond{1} > 1
        cond_text{1} = sprintf('There are %d conditions. Which should be tested here? \n\nChoose first condition(s):', model.Ncond{1});
        %cond_text{2} = 'Second condition(s) to test:';
        cond_tested = inputdlg(cond_text,'Which conditions should be tested?',1,{'1,2'});
        
        conds1 = str2num(cond_tested{1});
    end
    
    
    
    % Build selected ITC data into correct format
    Nsub = length(model.set_files);
    Nconds = [length(conds1)];
    
    total_elecs = length(model.defaults.chanlocs);
    
    
    
    
    
    
    Y1 = nan(total_elecs,size(Y,2),size(Y,3),Nsub*Nconds(1));
    %Y2 = nan(total_elecs,size(Y,2),size(Y,3),Nsub*Nconds(2));
    
    % Check cont length
    if length(cont_data) ~= size(Y1,4)
        error('Continuous data is of different length to subjects*Conditions analysed')
    end
    
    
    % Populate Y1 and Y2 with ITC data already saved
    
    
    for sub = 1:Nsub
        
        cd(model.itc_data{sub})
        
        load LIMO
        load Y
        
        
        for elec = 1:size(Y,1)
            
            org_elec = LIMO.data.chanlocs(elec).urchan;  % Find original elec index
            
            
            for cond = 1:Nconds(1)  % For each cond going into Y1
                Y1(org_elec,:,:,sub-1+cond) = Y(elec,:,:,conds1(cond));
            end
            
            % Leave the rest as nan
            
        end
    end
    
    size(Y1)
    
    
    Y1nans = mean(isnan(Y1(:)))
    

    
    % --- Check values
    if model.defaults.bootstrap == 1
        nboot = 1000;
    else
        nboot = 0;
    end
    
    tfce = model.defaults.tfce;
    if tfce == 1 && isfield(LIMO.data,'neighbouring_matrix') == 0  % Check we have neighb matrix. If not, create it.
        EEGLIMO.chanlocs = model.defaults.chanlocs;
        neighbdis = inputdlg('What neighbourhood distance should be used for TFCE neighbourhood matrix? (Perhaps 0.37 for 128 electrode systems)','Enter neighb distance',1,{'0.37'});
        neighbdis = str2num(neighbdis{1});
        [tmpneighbs, LIMO.data.neighbouring_matrix] = limo_get_channeighbstructmat(EEGLIMO,neighbdis);
    end
    
    
    Analysis_type = 'ITC';
    parameters = 1;
    

    % Run stats on this Ybig
    cd(original_LIMO_dir); 
    tag = sprintf('%s_conds%s',model.test_select{2},cond_tested{1});
    tag(~ismember(tag,['A':'Z' 'a':'z' '1':'9' '_']))= ''; % Clean up tag string
    mkdir(tag);cd(tag)
    LIMO.data.cond_tested = cond_tested;
    LIMO.data.chanlocs = model.defaults.chanlocs;
    save LIMO LIMO
    
    limo_random_robust(model.test_select{1},Y1,cont_data,parameters,nboot,tfce);
    disp([model.test_select{2} ' done'])
    
    
       % Plot
    plotnow = 1;
    if plotnow == 1
        load R2        
        
        limo_display_results_tf(LIMO, R2(:,:,:,3),1,['ITC R2 ', model.test_select{2}])
    end
    
    
elseif model.test_select{1} == 5 % ANOVA
    
    disp(model.test_select{2})
    disp('Not yet implemented')
    
    
elseif model.test_select{1} == 6 % Central tendancy?
    
    disp(model.test_select{2})
    disp('Not yet implemented')
    
    
end








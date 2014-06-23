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

global EEG
global LIMO

%% what to do




itc_prompt = sprintf('LIMO ITC allows loading of Inter-Trial Coherence (ITC) data and use of LIMO tools on that data. \n\n  It requires that that data has previously been generated and saved in EEG.etc.itc for each subject.');
%helpdlg(itc_prompt, 'Using LIMO ITC')

option='model specification';



% ITC gui - loading dataset list and trim info
[model.set_files,model.cat_files,model.cont_files,model.defaults] = limo_itc_gui;

for subject = 1:size(model.set_files,1)
    
    subpath = fileparts(model.set_files{subject});
    cd(subpath);
    
    limo_itc_import_data(model.set_files{subject},model.cat_files,model.cont_files,model.defaults);
    load('ITC_analysis/LIMO.mat')
    LIMO.Analysis = 'ITC';
    EEG=pop_loadset(LIMO.data.data);
    cd 'ITC_analysis'
    
    
    disp('loading ITC data...');
    
    % Let's treat 2-condition ITC data like TF data with 2 trials
    Yitc = EEG.etc.itc(:,LIMO.data.trim_low_f:LIMO.data.trim_high_f,LIMO.data.trim1:LIMO.data.trim2);
    
    
    if size(Yitc,1)/EEG.nbchan == 2  % If double electrode count in ITC data, check for 2 now
        disp('*** - Found double electrode count in ITC data - taking as two conditions')
        
        Y = nan(EEG.nbchan,size(Yitc,2),size(Yitc,3),2);
        Y(:,:,:,1) = abs(Yitc(1:EEG.nbchan,:,:));
        Y(:,:,:,2) = abs(Yitc(EEG.nbchan+1:end,:,:));
        
        %Y1_3d = abs(Yitc(1:EEG.nbchan,:,:));
        %Y2_3d = abs(Yitc(EEG.nbchan+1:end,:,:));
        
        %save('Y_3d', 'Y1_3d', 'Y2_3d')
        
        % Load a 4D Y with ITC data
        LIMO.data.size4D= size(Y);
        LIMO.data.size3D= [LIMO.data.size4D(1) LIMO.data.size4D(2)*LIMO.data.size4D(3) LIMO.data.size4D(4)];
        
        
    elseif size(Yitc,1)/EEG.nbchan == 1
        
        Y = nan(EEG.nbchan,size(Yitc,2),size(Yitc,3),1);
        Y(:,:,:,1) = abs(Yitc(1:EEG.nbchan,:,:));
        
        
    else
        error('Check ITC data length')
    end
    
    LIMO.model = model;
    
    save Y Y -v7.3;
    save LIMO LIMO -v7.3;
    model.itc_data{subject} = pwd;
    
    
    clear EEG Yitc
    
    
    
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
    Ncond = 2;
    total_elecs = length(model.defaults.chanlocs);
    
    
    disp(model.test_select{2})
    
    Ybig = nan(total_elecs,size(Y1_3d,2),size(Y1_3d,3),Nsub*Ncond);
    
    % Populate Ybig with ITC data already saved
    for sub = 1:Nsub
        
        cd(model.itc_data{sub})
        
        load LIMO
        load Y
        
        
        for elec = 1:size(Y,1)
            
            org_elec = LIMO.data.chanlocs(elec).urchan;  % Find original elec index
            
            Ybig(org_elec,:,:,sub*2-1) = Y(elec,:,:,1);  % Write low phase
            Ybig(org_elec,:,:,sub*2) = Y(elec,:,:,2);    % Write high phase
            
            % Leave the rest as nan
            
        end
    end
    
    % --- Check values
    if model.defaults.bootstrap == 1
        nboot = 1000;
    end
    
    tfce = model.defaults.tfce; Analysis_type = 'ITC';
    parameters = 1;
    
    
    % Run stats on this Ybig
    cd(original_LIMO_dir);
    limo_random_robust(model.test_select{1},Ybig,parameters,nboot,tfce);
    
    
    % Plot
    plot = 0;
    if plot == 1
        load one_sample_ttest_parameter_1
        
        
        limo_display_results_tf(LIMO, one_sample(:,:,:,4),1,['ITC ', model.test_select{2}])
    end
    
    
    
    
    
    
    
elseif model.test_select{1} == 2
    %% 2-samp t-test
    
    disp(model.test_select{2})
    
    disp(model.test_select{2})
    disp('Not yet implemented')
    
    
    
    
    
elseif model.test_select{1} == 3
    
    disp(model.test_select{2})
    disp('Not yet implemented')
    
    
elseif model.test_select{1} == 4
    
    disp(model.test_select{2})
    disp('Not yet implemented')
    
    
elseif model.test_select{1} == 5
    
    disp(model.test_select{2})
    disp('Not yet implemented')
    
    
end








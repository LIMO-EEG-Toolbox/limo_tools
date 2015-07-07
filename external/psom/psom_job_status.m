function [curr_status,tab_refresh] = psom_job_status(path_logs,list_jobs,mode_pipe,tab_refresh)
% Get the current status of a list of jobs.
%
% SYNTAX :
% CURR_STATUS = PSOM_JOB_STATUS(PATH_LOGS,LIST_JOBS,MODE_PIPE)
%
% _________________________________________________________________________
% INPUTS:
%
% PATH_LOGS
%   (string) the folder where the logs of a pipeline are stored.
%
% LIST_JOBS
%   (cell of strings) a list of job names
%
% MODE_PIPE
%   (string) the execution mode of the pipeline.
%   Possible values : 'session', 'batch', 'background', 'qsub', 'msub'
%   'condor'.
%
% TAB_REFRESH
%   (array, optional) TAB_REFRESH(N,:) is the latest known activity time for job 
%   LIST_JOBS{N}. If jobs become active and stay idle more than 30 seconds, 
%   they are declared failed (and a kill signal is sent). 
%
% _________________________________________________________________________
% OUTPUTS:
%
% CURR_STATUS
%   (cell of string) CURR_STATUS{K} is the current status of
%   LIST_JOBS{K}. Status can be :
%
%      'running' : the job is currently being processed.
%
%      'failed' : the job was processed, but the execution somehow
%         failed. That may mean that the function produced an
%         error, or that one of the expected outputs was not
%         generated. See the log file of the job for more info
%         using PSOM_PIPELINE_VISU.
%
%      'finished' : the job was successfully processed.
%
%      'none' : no attempt has been made to process the job yet 
%         (neither 'failed', 'running' or 'finished').
%
%      'absent' : there is no tag file and no job file. It looks like
%         the job name does not exist in the pipeline.
%
% TAB_REFRESH 
%    (array) same as the input, with refreshed active times for each job.
%
% _________________________________________________________________________
% COMMENTS: 
%
% The conditions for achieving a 'finished' status vary from one execution
% mode to another : 
%   
%   1. In 'session' mode, the 'finished' status is achieved when a
%   'finished' tag file can be found in the log folder. This means that the
%   job has completed without producing any error, and all the outputs have
%   been created.
% 
%   2. In 'batch' mode, in addition of the 'finished' tag file, an 'exit'
%   tag file is also required, indicating that the batch script has
%   terminated.
%
%  3. In 'qsub' mode, in addition to the 'finished' and 'exit' tag files,
%  the presence of a 'oqsub' and 'eqsub' log files are required, indicating
%  that the qsub job was properly terminated.
%   
% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008.
% Maintainer : pbellec@bic.mni.mcgill.ca
% See licensing information in the code.
% Keywords : pipeline

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

%% SYNTAX
if ~exist('path_logs','var') || ~exist('list_jobs','var') || ~exist('mode_pipe','var')
    error('SYNTAX: CURR_STATUS = PSOM_JOB_STATUS(PATH_LOGS,LIST_JOBS,MODE). Type ''help psom_job_status'' for more info.')
end

%% Table of refresh time
if nargin < 4
    tab_refresh = [];
end

%% Loop over all job names, and check for the existence of tag files
nb_jobs = length(list_jobs);
curr_status = cell([nb_jobs 1]);

for num_j = 1:nb_jobs        
    
    name_job = list_jobs{num_j};
        
    file_running   = [path_logs name_job '.running'];
    file_failed    = [path_logs name_job '.failed'];
    file_finished  = [path_logs name_job '.finished'];
    file_exit      = [path_logs name_job '.exit'];
    file_oqsub     = [path_logs name_job '.oqsub'];
    file_heartbeat = [path_logs name_job '.heartbeat.mat'];
            
    flag_exit      = psom_exist(file_exit);
    flag_oqsub     = psom_exist(file_oqsub);   
    flag_failed    = psom_exist(file_failed);
    flag_finished  = psom_exist(file_finished);            
    flag_running   = psom_exist(file_running); 
    flag_heartbeat = psom_exist(file_heartbeat); 
    
    if (flag_finished+flag_failed)>1
        error('I am confused : job %s has multiple tags. Sorry dude, I must quit ...',name_job);
    end          
                     
   if flag_finished
        
        switch mode_pipe
            case {'qsub','msub','condor'}
                if flag_oqsub&&flag_exit
                    curr_status{num_j} = 'finished';
                else
                    curr_status{num_j} = 'running';
                end
            case 'batch'
                if flag_exit
                    curr_status{num_j} = 'finished';
                else
                    curr_status{num_j} = 'running';
                end
            otherwise
                curr_status{num_j} = 'finished';
        end

    elseif flag_failed

        switch mode_pipe
            case {'qsub','msub','condor'}
                if flag_oqsub&&flag_exit
                    curr_status{num_j} = 'failed';
                else
                    curr_status{num_j} = 'running';
                end
            case 'batch'
                if flag_exit
                    curr_status{num_j} = 'failed';
                else
                    curr_status{num_j} = 'running';
                end
            otherwise
                curr_status{num_j} = 'failed';
        end

    elseif flag_running

        %% Try to access refresh times
        if flag_heartbeat
            if isempty(tab_refresh)
                curr_status{num_j} = 'running';
            elseif any(tab_refresh(num_j,:,1)<0)
                % this is the first time an active time is collected
                % simply update tab_refresh
                tab_refresh(num_j,:,1) = 0;
                curr_status{num_j} = 'running';
            else
                try
                    refresh_time = load(file_heartbeat);
                    test_change = etime(refresh_time.curr_time,tab_refresh(num_j,:,1))>1;
                catch
                    % The heartbeat is unreadable
                    % Assume this is a race condition
                    % Consider no heartbeat was detected
                    test_change = false;
                end

                if test_change
                    % I heard a heartbeat!    
                    tab_refresh(num_j,:,1) = refresh_time.curr_time;
                    tab_refresh(num_j,:,2) = clock;
                    curr_status{num_j} = 'running';
                else 
                    % how long has it been without a heartbeat?
                    elapsed_time = etime(clock,tab_refresh(num_j,:,2));
                    if elapsed_time > 30
                        % huho 30 seconds without a heartbeat, he's dead Jim
                        curr_status{num_j} = 'failed';
                    else
                        curr_status{num_j} = 'running';
                    end
                end
            end
        else
            curr_status{num_j} = 'running';
        end

    else
        
        curr_status{num_j} = 'submitted';        
                
    end
end

% Demo script for the pipeline system for Octave and Matlab (PSOM).
%
% The blocks of code follow the following tutorial: 
% http://simexp.github.io/psom/how_to_use_psom.html
%
% The demo will create and populate a folder 'psom_demo' in the current direction. 
% To restart the demo from scratch, simply delete this folder.
%
% In matlab, you can run a specific block of code by selecting it and press F9, 
% or by putting the cursor anywhere in the block and press CTRL+ENTER.
% Otherwise just copy paste the code in the command window. 
% Please make sure PSOM is in the matlab/octave search path.
%
% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008-2010.
% Departement d'informatique et de recherche operationnelle
% Centre de recherche de l'institut de Geriatrie de Montreal
% Universite de Montreal, 2011-2014
% Maintainer : pierre.bellec@criugm.qc.ca
% See licensing information in the code.
% Keywords : pipeline, PSOM, demo

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% An example of (toy) pipeline %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% where to generate the outputs of the demo
clear, path_demo = [pwd filesep 'demo_psom' filesep]; 

% Job "sample" :    No input, generate a random vector a
command = 'a = randn([opt.nb_samps 1]); save(files_out,''a'')';
pipeline.sample.command      = command;
pipeline.sample.files_out    = [path_demo 'sample.mat'];
pipeline.sample.opt.nb_samps = 10;

% Job "quadratic" : Compute a.^2 and save the results
command = 'load(files_in); b = a.^2; save(files_out,''b'')';
pipeline.quadratic.command   = command;
pipeline.quadratic.files_in  = pipeline.sample.files_out;
pipeline.quadratic.files_out = [path_demo 'quadratic.mat']; 

% Adding a job "cubic" : Compute a.^3 and save the results
command = 'load(files_in); c = a.^3; save(files_out,''c'')';
pipeline.cubic.command       = command;
pipeline.cubic.files_in      = pipeline.sample.files_out;
pipeline.cubic.files_out     = [path_demo 'cubic.mat']; 

% Adding a job "sum" : Compute a.^2+a.^3 and save the results
command = 'load(files_in{1}); load(files_in{2}); d = b+c, save(files_out,''d'')';
pipeline.sum.command       = command;
pipeline.sum.files_in{1}   = pipeline.quadratic.files_out;
pipeline.sum.files_in{2}   = pipeline.cubic.files_out;
pipeline.sum.files_out     = [path_demo 'sum.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize the dependency graph %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%options
psom_visu_dependencies(pipeline)

%%%%%%%%%%%%%%%%%%%%%%
%% Run the pipeline %%
%%%%%%%%%%%%%%%%%%%%%%
opt.path_logs = [path_demo 'logs' filesep]; % The logs folder
opt.mode = 'background';                    % Execute jobs in the background
opt.max_queued = 2;                         % Use two threads
psom_run_pipeline(pipeline,opt);            % Now, actually run the pipeline

%%%%%%%%%%%%%%%%%%
%% Change a job %%
%%%%%%%%%%%%%%%%%%
pipeline.quadratic.command = 'BUG!'; % Change the job quadratic to introduce a bug
psom_run_pipeline(pipeline,opt)      % Restart the pipeline

psom_pipeline_visu(opt.path_logs,'log','quadratic'); % Visualize the log file of the failed job

pipeline.quadratic.command  = 'load(files_in); b = a.^2; save(files_out,''b'')'; % fix the bug 
psom_run_pipeline(pipeline,opt); % restart the pipeline

%%%%%%%%%%%%%%%
%% Add a job %%
%%%%%%%%%%%%%%%
pipeline.cleanup.command = 'delete(files_clean)'; % Clean up intermediate files
pipeline.cleanup.files_clean = pipeline.sample.files_out; % New job field for files deleted by the job: files_clean
psom_run_pipeline(pipeline,opt); % Restart the pipeline
psom_visu_dependencies(pipeline); % Visualize the new dependecy graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Manually restarting a job %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.restart = {'quadratic'};
psom_run_pipeline(pipeline,opt);

%%%%%%%%%%%%%%%%%%%%%%%%
%% Monitor a pipeline %%
%%%%%%%%%%%%%%%%%%%%%%%%

%% Display flowchart
psom_pipeline_visu(opt.path_logs,'flowchart');

%% List the finished jobs
psom_pipeline_visu(opt.path_logs,'finished');

%% Display log
psom_pipeline_visu(opt.path_logs,'log','quadratic');

%% Display Computation time
psom_pipeline_visu(opt.path_logs,'time','');

%% Monitor history
psom_pipeline_visu(opt.path_logs,'monitor');
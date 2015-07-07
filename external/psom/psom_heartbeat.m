function [] = psom_heartbeat(file_heart,file_kill,pid)
% This function is internal to PSOM and not meant to be used directly.
%
% PSOM_HEARTBEAT(FILE_HEART,FILE_KILL,PID)
%
% FILE_HEART (string) the name of a .mat file. The following variables will be updated every 5s
%    inside the .mat file:
%    CURR_TIME (vector) the output of clock
%    TELAPSED (scalar) the time (s) elapsed since the heartbeat was started.
% FILE_KILL (string) the name of a file. If this file is detected at any point in time 
%    the function will kill the process PID and exit.
% PID (scalar) a process ID.
%
% See licensing information in the code.
%  
%system(['octave --eval "cd /home/pbellec/, build_path_std, cd /home/pbellec/tmp/tmp, psom_heartbeat(''toto.mat'',''tata.mat'',' num2str(getpid) '); exit"'],false,'async')

% Copyright (c) Pierre Bellec, Montreal Neurological Institute, 2008-2010.
% Departement d'informatique et de recherche operationnelle
% Centre de recherche de l'institut de Geriatrie de Montreal
% Universite de Montreal, 2011
% Maintainer : pierre.bellec@criugm.qc.ca
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
tic;
flag_beat = true;
while flag_beat
    if exist('OCTAVE_VERSION','builtin')  
        [err,msg] = kill(pid,0); % use the kill octave command
    else
        [err,msg] = system(sprintf('kill %i',pid)); % kill is not available, try a system call
    end
    flag_beat = err==0;
    curr_time = clock;
    save(file_heart,'curr_time');
    if exist(file_kill,'file')
        if exist('OCTAVE_VERSION','builtin')  
            kill(pid,9)
        else
            system(sprintf('kill -9 %i',pid));
        end 
        exit
    end
    if exist('OCTAVE_VERSION','builtin')  
        [res,msg] = system('sleep 5');
    else
        sleep(5); 
    end
end
    
function [subname,valid] = limo_get_subname(stringin,mode)

% returns the subject name from the string
%
% FORMAT subname = limo_get_subname(stringin)
%        [subname,valid] = limo_get_subname(stringin,'silent')
%
% INPUT stringin can be a path that includes the BIDS sub- something
%                       a filename sub-sthg_xxx_eeg.set for instance
%                       a full filename that includes the path (in which
%                       case a check for consistency of the sub in path and
%                       in filename is also peformed)
%       mode if a full filename is passed, specifies if 'interruptive'
%       (default) meaning an error pops-up or 'silent' and users is
%       expected to deal with a possible error based on output
%
% OUTPUT subname is the sub- name (empty if not a BIDS path/filename)
%        valid is 0/1 testing if path and fiename have the same sub- name
%
% ------------------------------------------------------------------
% Cyril R. Pernet - Copyright (C) LIMO Team 2025

%% defaults
subname = [];
if nargin == 1
    mode = 'interruptive';
end

%% compute
if isfolder(stringin)
    if contains(stringin,'sub-')
        subname = extractAfter(stringin,'sub-');
        subname = ['sub-' subname(1:min(strfind(subname,filesep))-1)];   
    end
elseif isfile(stringin)
    [fpath,fname] = fileparts(stringin);
    if contains(fname,'sub-')
        subname = extractAfter(fname,'sub-');
        subname = ['sub-' subname(1:min(strfind(subname,'_'))-1)];   
    end
    if ~isempty(fpath)
        subname_from_path = limo_get_subname(fpath);
        if isempty(subname)
            invalid = 1;
        else
            invalid = strcmp(subname_from_path,subname);
        end
        if strcmpi(mode,'interruptive') && invalid == 0
            limo_errordlg('conflict in sub- name between path and file')
            error('%s',subname)
        end
    else
        
    end        
else
     limo_errordlg('the string passed in is not valid')
     error('%s',stringin)
end

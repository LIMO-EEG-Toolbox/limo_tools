function channeighbstructmat = limo_expected_chanlocs(varargin)

% This function loads an EEG dataset to create a file with the
% location of all expected channels and to create a neighbourhood
% distance matrix used to control for multiple comparisons.
%
% FORMAT
% limo_expected_chanlocs
% limo_expected_chanlocs(data set name, path)
% limo_expected_chanlocs(data set name, path,neighbour distance)
%
% INPUT
% data set name is the name of a eeglab.set
% path is the location of that file
% neighbour distance is the distance between channels to buid the
% neighbourhood matrix
%
% channeighbstructmat is the neighbourhood matrix
%
% See also LIMO_NEIGHBOURDIST LIMO_GET_CHANNEIGHBSTRUCMAT
%
% Guillaume Rousselet v1 11 June 2010
% Rewritten by Cyril Pernet so we don't have to know which subject has the
% largest channel description v2 16 July 2010 - update 18 July 2012 to get
% output channeighbstructmat so we can update subjects for tfce
% ----------------------------------
%  Copyright (C) LIMO Team 2010

 current_dir = pwd;

 % ask if data are from one subject or a set then get data
 % ---------------------------------------------------------
 if nargin == 0
     quest = questdlg('Make the Expected Chanlocs file from 1 subject or search throughout a set of subjects?','Selection','Set','One','Set');
     FileName = []; PathName = []; FilterIndex = [];
 elseif nargin >= 2
     quest = 'One';
     FileName = varargin{1}; PathName = varargin{2}; FilterIndex = 1;
 else
     error('wrong number of arguments')
 end
 
 neighbourdist = [];
 if nargin == 3
     neighbourdist = varargin{3};
 end
 
 % from 1 subject
 % -----------------------
 if strcmp(quest,'One')
     
     if isempty(FileName)
         [FileName,PathName,FilterIndex]=uigetfile('*.set','EEGLAB EEG dataset before electrode removal');
         if FilterIndex == 0
             return
         end
     end
     cd(PathName)
     
     try
         EEG=pop_loadset(FileName);
         expected_chanlocs = EEG.chanlocs;
         fprintf('Data set %s loaded',FileName); disp(' ')
         if isempty(neighbourdist)
             neighbourdist = eval(cell2mat(inputdlg('enter neighbourhood distance','neighbourhood distance'))); % 0.37 for biosemi 128;
         end
         [neighbours,channeighbstructmat] = limo_get_channeighbstructmat(EEG,neighbourdist);
         if sum(channeighbstructmat(:)) == 0
            msg = sprintf('the neighbouring matrix is empty, it''s likely a distance issue \n see imo_ft_neighbourselection.m');
            error(msg)
         end
         cd (current_dir); 
         if nargout == 0
             save expected_chanlocs expected_chanlocs channeighbstructmat % save all in one file
             fprintf('expected_chanlocs & channeighbstructmatfile saved\n');
         end
     catch ME
         errordlg('pop_loadset eeglab function not found','error');
     end
     
     
 else   % more likely from a set of subjects
     % -------------------------------------------
     
     % get data
     go = 1; index = 1;
     while go == 1
         [name,path] = uigetfile('LIMO.mat',['select LIMO file subject ',num2str(index),go]);
         if name == 0
             go = 0;
         else
             if ~strcmp(name,'LIMO.mat')
                 error(['you selected the file ' name ' but a LIMO.mat file is expected']);
             else
             Names{index} = name;
             Paths{index} = path;
             Files{index} = sprintf('%s\%s',path,name);
             cd(path); cd ..
             index = index + 1;
             end
         end
     end
     
     % retreive all chanlocs
     for i=1:length(Paths)
         cd(Paths{i})
         load LIMO
         chanlocs{i} = LIMO.data.chanlocs;
         size_chanlocs(i) = size(LIMO.data.chanlocs,2);
         clear LIMO
     end
     
%      % simply use the largest .set
     [v,loc]=max(size_chanlocs);
     cd(Paths{loc})
     local_files = dir;
     for i=1:length(local_files)
         name = local_files(i).name;
         if length(name)>4
             name = name(end-2:end);
             if strcmp(name,'set')
                 data_file_name = local_files(i).name;
             end
         end
     end
     
     % find unique labels across subjects
     % make a cap with those unique channels
     % LIMO handles NaN so that not all electrodes
     % have to be present in each subject :-)
%      [n,ref]=max(size_chanlocs);
%      caps = [1:size(chanlocs,2)];
%      caps(find(size_chanlocs - n));
%      for i=1:n
%          for c=caps
%              chanlocs{ref}(i).labels
%              chanlocs{c}.labels
%          end
%      end
     
      
     
     cd (current_dir); 
     % now we have 1 subject cap we can do as above
     limo_expected_chanlocs(data_file_name,Paths{loc})
     
 end




    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

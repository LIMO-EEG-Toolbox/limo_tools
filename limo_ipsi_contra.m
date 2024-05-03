function limo_ipsi_contra(varargin)

% routine taking 2 conditions (beta, con) and contrasting contralateral vs ipsilateral cases
% contralateral = (right hemisphere channels condition 1 + left hemisphere channels for condition 2) / 2
% ipsilateral   = (right hemisphere channels for condition 2 + left hemisphere channels for condition 1) / 2
%
% those type of analyses are common for lateralized experiments like e.g.
% N2pc https://en.wikipedia.org/wiki/N2pc
%
% FORMAT LIMO_LI(betafiles,'parameters',[A B],'chanlocs',chanstruct,'channelpairs',channelindices)
%        LIMO_LI(confiles1,confiles2,'chanlocs',chanstruct,'channelpairs',channelindices)
%
% INPUTS List of beta files + parameters (ie column of the design matrix) or two lists of con files
%       'channelpairs' key matches the channelindices values, a n*2 matrix for chanloc to pair
%                      this is can be obtained using LIMO_pair_channels.m
%       'chanlocs' key matches the EEGLAB/LIMO channel structure and neighbouring matrix
%
% OUTPUTS contralateral and ipsilateral files are created on the drive 
%         and a paired t-test is computed between those files
%
% see also LIMO_pair_channels, LIMO_random_robust, LIMO_yuend_ttest
%
% Cyril Pernet Decembre 2023

LIMO.dir  = pwd;
LIMO.Type = 'Channels'; % in theory we could do source as well

if nargin == 6 % con files
    [Names{1},Paths{1},LIMO.data.data{1}] = limo_get_files([],[],[],varargin{1});
    LIMO.data.data_dir{1}                 = Paths{1};
    [Names{2},Paths{2},LIMO.data.data{2}] = limo_get_files([],[],[],varargin{2});
    LIMO.data.data_dir{1}                 = Paths{2};
    if isempty(Names{1}) || isempty(Names{2})
        LIMO_warndlg('Could not parse files names - function aborded')
        return
    end

    for n=3:2:6
        if contains(varargin{n},'loc','IgnoreCase',true)
            LIMO.data.chanlocs            = varargin{n+1}.expected_chanlocs;
            LIMO.data.neighbouring_matrix = varargin{n+1}.channeighbstructmat;
        elseif contains(varargin{n},'pair','IgnoreCase',true)
            channelindices = varargin{n+1};
        else
            LIMO_errordlg('expecting text as 3rd and 5th input given the number of arguments in')
        end
    end
    LIMO.design.parameters = [1 2];

elseif nargin == 7 % beta files
    [~,Paths,Files] = limo_get_files([],[],[],varargin{1});
    for n=2:2:7
        if contains(varargin{n},'parameter','IgnoreCase',true)
            LIMO.design.parameters  = varargin{n+1};
        elseif contains(varargin{n},'loc','IgnoreCase',true)
            LIMO.data.chanlocs            = varargin{n+1}.expected_chanlocs;
            LIMO.data.neighbouring_matrix = varargin{n+1}.neighbouring_matrix;
        elseif contains(varargin{n},'pair','IgnoreCase',true)
            channelindices = varargin{n+1};
        else
            LIMO_errordlg('expecting text as 2nd, 4th and 6th nputs given the number of arguments in')
        end
    end
end

%% iterate through subjects to create files
% contralateral = (right hemisphere channels condition 1 + left hemisphere channels for condition 2) / 2
% ipsilateral   = (right hemisphere channels for condition 2 + left hemisphere channels for condition 1) / 2

if exist('Files','var') % betas
    LIMO.data.data     = Files;
    LIMO.data.data_dir = Paths;
    error('not implemented for Betas yet')

else % con
    [first_frame,last_frame,subj_chanlocs,LIMO] = limo_match_frames(Paths{1},LIMO);
    % LIMO.data.start,LIMO.data.end
    for subject = length(LIMO.data.data{1}):-1:1

         tmp = load(fullfile(Paths{1}{subject} ,'LIMO.mat'));
         if tmp.LIMO.Type ~= LIMO.Type
             limo_errordlg('Only Channel type of analysis supported')
         end
        tmp = load(LIMO.data.data{1}{subject});
        if size(tmp.con,2) == length(first_frame(subject):last_frame(subject))
            first = 1; last = size(tmp.con,2);
        else
            first = first_frame(subject); last = last_frame(subject);
        end
        data1(:,:,subject) = limo_match_elec(subj_chanlocs(subject).chanlocs,...
            LIMO.data.chanlocs,first,last,...
            squeeze(tmp.con(:,:,1)));
        tmp = load(LIMO.data.data{2}{subject});
        if size(tmp.con,2) == length(first_frame(subject):last_frame(subject))
            first = 1; last = size(tmp.con,2);
        else
            first = first_frame(subject); last = last_frame(subject);
        end
        data2(:,:,subject) = limo_match_elec(subj_chanlocs(subject).chanlocs,...
            LIMO.data.chanlocs,first,last,...
            squeeze(tmp.con(:,:,1)));
    end

    contralateral = (data1(channelindices(:,2),:,:)+data2(channelindices(:,1),:,:))./2;
    ipsilateral   = (data1(channelindices(:,1),:,:)+data2(channelindices(:,2),:,:))./2;
    save(fullfile(LIMO.dir,'contralateral.mat'),'contralateral')
    save(fullfile(LIMO.dir,'ipsilateral.mat'),'ipsilateral')
end

%% T-test
LIMO.Level                    = 2;
for n = size(channelindices,1):-1:1
    labels{n} = [LIMO.data.chanlocs(channelindices(n,1)).labels '/' ...
        LIMO.data.chanlocs(channelindices(n,2)).labels];
end
LIMO.data.chanlocs            = LIMO.data.chanlocs(channelindices(:,1)); % half now
for n = size(channelindices,1):-1:1
    LIMO.data.chanlocs(n).labels = labels{n};
end
LIMO.data.neighbouring_matrix = LIMO.data.neighbouring_matrix(channelindices(:,1),channelindices(:,1));
LIMO.design.name              = ['Paired t-test all ' LIMO.Type(1:end-1)];
LIMO.design.bootstrap         = 1000;
LIMO.design.tfce              = 1;
save(fullfile(LIMO.dir,'LIMO.mat'),'LIMO')
limo_random_robust(3,contralateral,ipsilateral,[1 2],LIMO);


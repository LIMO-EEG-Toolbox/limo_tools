function limo_create_single_trials(EEG,varargin)

% FORMAT limo_create_single_trials(EEG,options)
%
% example: 
%     options = {'format','matrix', 'datatype', 'channels', ...
%     'erp','off','spec','on','ersp','off','itc','on','rmicacomps','on', ...
%     'erpparams',[],'specparams',[],'erspparams',[], 'interp','off', ...
%     'scalp','off','recompute','on','savetrials','on'};
%     limo_create_single_trials('my_subject.set',options{:})
%
% INPUT EEG is a EEG.set of a single subject from EEGLAB
%       option is a structure with various fiedls depending on what data
%       type needs to be computed -- the idea here is that the data in
%       EEG.data are the 'raw' data and we here apply IC removals and then
%       create output -- see std_precomp
%
%       option.format: 'matrix' (default) or 'cell'
%       option.datatype: 'channels' or 'components'
%
%  Optional fields
%    'erp'      - ['on'|'off'] pre-compute ERPs for each dataset (default 'on')
%    'spec'     - ['on'|'off'] pre-compute spectrum for each dataset (default 'on')
%                 Use 'specparams' to set spectrum parameters.
%    'ersp'     - ['on'|'off'] pre-compute ERSP for each dataset (default 'on')
%                 Use 'erspparams' to set time/frequency parameters.
%    'itc'      - ['on'|'off'] pre-compute ITC for each dataset (default 'on')
%                 Use 'erspparams' to set time/frequency parameters.
%    'erpparams'   - [cell array] Parameters for the std_erp function. See
%                    std_erp for more information.
%    'specparams'  - [cell array] Parameters for the std_spec function. See
%                    std_spec for more information.
%    'erspparams'  - [cell array] Optional arguments for the std_ersp function.
%    'rmicacomps'  - ['on'|'off'|'processica'] remove ICA components pre-selected in
%                    each dataset (EEGLAB menu item, "Tools > Reject data using ICA
%                    > Reject components by map). This option is ignored when
%                    precomputing measures for ICA clusters. Default is 'off'.
%                    'processica' forces to process ICA components instead of
%                    removing them (default 'on')
%
% OUTPUT file_names the names of the files creates
%        + single trials data files stored on the drive
%
% see also std_precomp
%
% Author: Cyril Pernet (LIMO Team), The university of Edinburgh, 2014
%         Arnaud Delorme (EEGLAB), SCCN, 2014
%
% ------------------------------------------
% Copyright (C) LIMO Team 2014

%% deal with options
% set defaults
try
    options = varargin;
    if ~isempty( varargin ),
        for i = 1:2:numel(options)
            opt.(options{i}) = options{i+1};
        end
    else opt= []; end;
catch
    error('limo_create_single_trials() error: calling convention {''key'', value, ... } error'); return;
end;

try, opt.format;           catch, opt.format        = 'matrix';     end;
try, opt.datatype;         catch, opt.datatype      = 'channels';   end;
try, opt.erp;              catch, opt.erp           = 'on';        end;
try, opt.spec;             catch, opt.spec          = 'on';        end;
try, opt.ersp;             catch, opt.ersp          = 'on';        end;
try, opt.itc;              catch, opt.itc           = 'on';        end;
try, opt.rmicacomps;       catch, opt.rmicacomps    = 'on';         end;
try, opt.interp;           catch, opt.interp        = 'off';        end;

try, opt.erpparams;        catch, opt.erpparams     = [];          end;
try, opt.specparams;       catch, opt.specparams    = [];          end;
try, opt.erspparams;       catch, opt.erspparams    = [];          end;

if isfield(options,'itc') && isfield(options,'ersp')
    opt.itc = options.ersp; % for consitancy use ersp
end

if isfield(options,'ersp') && ~isfield(options,'itc')
    opt.itc = 'off';
end

opt.savetrials     = 'on';
opt.scalp          = 'off';
opt.recompute      = 'on';

fields = fieldnames(opt);
c = 1;
for i = 1: length(structfun(@numel,opt))
    if ~any ([isempty(eval(['opt.' fields{i}])),strcmp(fields(i),'format'), strcmp(fields(i),'datatype')])
    in_options{c}   = fields{i};
    in_options{c+1} = eval(['opt.' fields{i}]);
    c = c+2;
    end
end


%% create a STUDY just for that subject
if ~isstruct(EEG)
    ALLEEG = pop_loadset(EEG);
else
    ALLEEG = EEG; clear EEG;
end

STUDY = struct;
[trash,name]=fileparts(ALLEEG.filepath); clear trash
[STUDY ALLEEG] = std_editset(STUDY, ALLEEG, 'commands',{{'index' 1 'subject' name}},'updatedat','off' );
STUDY.name = ALLEEG.filename; [STUDY ALLEEG] = std_checkset(STUDY, ALLEEG);
name = [ALLEEG.filepath filesep ALLEEG.filename '_single_trials'];
STUDY.design.cell.filebase = name;

%% compute
[STUDY ALLEEG] = std_precomp(STUDY,ALLEEG,opt.datatype,in_options{:});
clear STUDY

%% save-delete data / update EEG.set file
if strcmp(opt.datatype,'channels')
     
    if strcmp(opt.erp,'on')
        data = load('-mat',[name '.daterp']);
        ALLEEG.etc.timeerp = data.times;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_daterp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.daterp = [name '_daterp.mat'];
            delete([name '.daterp']);
        else
            ALLEEG.etc.datafiles.daterp = [name '.daterp'];
        end
    end
    
    if strcmp(opt.spec,'on')
        if strcmp(opt.format,'matrix')
            [data,trash,ALLEEG.etc.freqspec] = limo_struct2mat([name '.datspec']);
            save([name '_datspec.mat'],'data'); clear data trash
            ALLEEG.etc.datafiles.datspec = [name '_datspec.mat'];
            delete([name '.datspec']);
        else
            data = load('-mat',[name '.datspec']);
            ALLEEG.etc.freqspec = data.freqs;
            ALLEEG.etc.datafiles.datspec = [name '.datspec'];
        end
    end
        
    
    if strcmp(opt.ersp,'on')
        disp('reading single trials ersp, be patient ..')
        if strcmp(opt.format,'matrix')
            [data,ALLEEG.etc.timeersp,ALLEEG.etc.freqersp] = limo_struct2mat([name '.dattimef']);
            save([name '_datersp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datersp = [name '_datersp.mat'];
            delete([name '.dattimef']);
        else
            data = load('-mat',[name '.dattimef']);
            ALLEEG.etc.timeersp = data.times;
            ALLEEG.etc.freqersp = data.freqs;
            ALLEEG.etc.datafiles.datersp = [name '.dattimef'];
        end
    end
    
    if strcmp(opt.itc,'on')
        data = load('-mat',[name '.datitc']);
        ALLEEG.etc.timeitc = data.times;
        ALLEEG.etc.freqitc = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_datitc.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datitc = [name '_datitc.mat'];
            delete([name '.datitc']);
        else
            ALLEEG.etc.datafiles.datitc = [name '.datitc'];
        end
    end
    
end


if strcmp(opt.datatype,'ica')

    if strcmp(opt.erp,'on')
        data = load('-mat',[name '.icaerp']);
        ALLEEG.etc.timeerp = data.times;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaerp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.daterp = [name '_icaerp.mat'];
            delete([name '.icaerp']);
        else
            ALLEEG.etc.datafiles.icaerp = [name '.icaerp'];
        end
    end
    
    if strcmp(opt.spec,'on')
        data = load('-mat',[name '.icaspec']);
        ALLEEG.etc.freqspec = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaspec.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datspec = [name '_icaspec.mat'];
            delete([name '.icaspec']);
        else
            ALLEEG.etc.datafiles.icaspec = [name '.icaspec'];
        end
    end
        
    
    if strcmp(opt.ersp,'on')
        data = load('-mat',[name '.icatimef']);
        ALLEEG.etc.timeersp = data.times;
        ALLEEG.etc.freqersp = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaersp.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datersp = [name '_icaersp.mat'];
            delete([name '.icatimef']);
        else
            ALLEEG.etc.datafiles.icaersp = [name '.icatimef'];
        end
    end
    
    if strcmp(opt.itc,'on')
        data = load('-mat',[name '.icaitc']);
        ALLEEG.etc.timeitc = data.times;
        ALLEEG.etc.freqitc = data.freqs;
        if strcmp(opt.format,'matrix')
            data = limo_struct2mat(data);
            save([name '_icaitc.mat'],'data'); clear data
            ALLEEG.etc.datafiles.datitc = [name '_icaitc.mat'];
            delete([name '.icaitc']);
        else
            ALLEEG.etc.datafiles.icaitc = [name '.icaitc'];
        end
    end
end

EEG = ALLEEG; clear ALLEEG;
pop_saveset(EEG, 'filename', EEG.filename, 'filepath',EEG.filepath,'savemode' ,'twofiles');


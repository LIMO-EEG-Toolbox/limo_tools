function limo_create_files(varargin)

% creates EEG files from design matrix
%
% INPUTS
% limo_creates_files(File,parameters,type)
%                    File name (with path) of the LIMO file to use
%                    parameters list the regressors to use
%                    type should be 'raw' or 'modeled'
%
% Cyril Pernet v1 26-May-2010
% ---------------------------
% Copyright (C) LIMO Team


%% simple checking

if isempty(varargin)
    [file,dir] = uigetfile('LIMO.mat','select a LIMO.mat file');
    cd (dir); load LIMO.mat; X = LIMO.design.X;
    parameters = eval(cell2mat(inputdlg('which parameters to test e.g [1:3]','parameters option')));
    if isempty(parameters)
        return
    elseif max(parameters) >= size(X,2)
        errordlg('incorrect parameter(s)');return
    end
    type = questdlg('File option','type of data to write?','raw','modeled','raw');
    if isempty(type)
        return;
    end

else
    load(varargin{1}); % LIMO file
    parameters = varargin{2};
    type = varargin{3};
    if ~strcmp(type,'raw') || ~strcmp(type,'modeled')
        error('type of output not recognized') ;
    end
end

%% create files

if strcmp(type,'raw')
    load Yr
    for i=parameters
        if i <=  LIMO.design.nb_conditions
            name = sprintf('raw_data_regressor_%g',i);
            fprintf('making file %s ... ',name); disp(' ');
            Reg = squeeze(X(:,i));
            for j=1:size(Yr,1)
                for k=1:size(Yr,2)
                    R(j,k,:) = Reg(find(Reg)).* squeeze(Yr(j,k,find(Reg)));
                end
            end
            save ([name],'R'); clear Reg R
        else
            errordlg(['no raw data are saved for continuous regressors - parameter ',num2str(i)]);
        end
    end
else % strcmp(type,'modeled')
    load Betas
    for i=parameters
        if i <size(X,2)
            name = sprintf('modeled_data_regressor_%g',i);
            fprintf('making file %s ... ',name); disp(' ');
            Reg = squeeze(LIMO.design.X(:,[i end]));
            for j=1:size(Betas,1)
                for k=1:size(Betas,2)
                    M(j,k,:) = Reg(:,1).*squeeze(Betas(j,k,i)) + Reg(:,2).*squeeze(Betas(j,k,end)); % model of param + model of cst
                end
            end
            save ([name],'M'); clear Reg M
        end
    end
end

disp('file creation done')


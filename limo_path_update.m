function limo_path_update

% silly routine to update all LIMO path in LIMO.mat files
% By default it updates the LIMO.mat to it's new path and look for the 
% associated data set either in the same directory or above - in cases
% where it is not found, a question dialogue pops out to ask if you want to 
% select the relevant data set
%
% Cyril Pernet v2 13-Sept-2010
% -----------------------------
%  Copyright (C) LIMO Team 2010


source_dir = pwd;

% 1st get the LIMO.mat files
% --------------------------

go = 1; index = 1;
while go == 1
    [name,path] = uigetfile('LIMO.mat',['select LIMO file subject ',num2str(index),go]);
    if name == 0
        go = 0;
    else
        Names{index} = name;
        Paths{index} = path;
        Files{index} = sprintf('%s\%s',path,name);
        cd(path); cd ..
        index = index + 1;
    end
end

% 2nd update the path and look for .set
% -------------------------------------
for i=1:size(Paths,2)
    cd(Paths{i})
    load LIMO
    
    % same directory for LIMO.dir and data
    if strcmp(LIMO.dir,LIMO.data.data_dir) || strcmp(LIMO.dir,LIMO.data.data_dir(1:end-1)) % happens because from EEGLAB there is an extra /
        LIMO.dir = pwd;
        LIMO.data.data_dir = pwd;
        name = dir(LIMO.data.data);
        if ~isempty(name)
            file_found(i) = 1;
            fprintf('LIMO file subject %g successfully updated\n',i);
        else
            file_found(i) = 0;
        end
        save LIMO LIMO
        
    % different directory
    elseif length(LIMO.dir) > length(LIMO.data.data_dir)
        test = LIMO.dir(1:length(LIMO.data.data_dir)) == LIMO.data.data_dir;
        if sum(test) == length(LIMO.data.data_dir)  % data are located in a directory above LIMO
            current_dir = pwd; LIMO.dir = current_dir;
            [LIMO.data.data_dir,name,ext]=fileparts(current_dir); 
            cd(LIMO.data.data_dir)
            name = dir(LIMO.data.data);
            if ~isempty(name)
                file_found(i) = 1;
                fprintf('LIMO file subject %g successfully updated\n',i);
            else
                file_found(i) = 0;
            end
        end
        cd(LIMO.dir); save LIMO LIMO
                
    else
        file_found(i) = 0;
    end
end


% 3rd update missing .set
% ------------------------
if sum(file_found) ~= length(file_found)
    q = questdlg('.set missing - do you want to select data or simply update the LIMO file','association between file missing','Only update LIMO','Let me select .set','Only update LIMO');
    if strcmp(q,'Let me select .set')
        for i=1:length(file_found)
            if file_found(i) == 0
                cd(Paths{i}); load LIMO
                [name,path] = uigetfile('*.set','data set missing - please select data' );
                % check
                if ~strcmp(name,LIMO.data.data)
                    errordlg('a different set was previously assicoated to this LIMO.mat file, please check')
                else
                    cd(path); LIMO.data.data_dir = pwd;
                    cd(Paths{i}); save LIMO LIMO
                end
            end
        end
    end
end

cd(source_dir)


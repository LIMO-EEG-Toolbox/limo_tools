function [Names,Paths,Files] = limo_get_files(gp)

% routine to get multifiles from different directories
% Name, path, Files are returned as cells
%
% Cyril Pernet 25-08-09
% Ncolas Chauveau 07-04-11 - allows txt files
% --------------------------------------------
%  Copyright (C) LIMO Team 2010

go = 1; index = 1;
while go == 1
    try
        [name,path] = uigetfile('*.mat; *.txt',['select beta or con file subject ',num2str(index),' ',gp,' or list file']);
    catch ME
        [name,path] = uigetfile('*.mat; *.txt',['select beta or con file subject or list file']);
    end
    
    if name == 0
        go = 0; % exit
    elseif name(end-2:end)=='mat'  % select mat files
        Names{index} = name;
        Paths{index} = path;
        % Files{index} = sprintf('%s\%s',path,name);
        Files{index} = [path name];
        cd(path); cd ..
        index = index + 1;
    else
        groupe_files=textread([path name],'%s'); % ,'whitespace','');  % select a txt file listing all mat files
        for f=1:size(groupe_files,1)
            Files{f} = groupe_files{f};
            [Paths{f},NAME,EXT] = fileparts(groupe_files{f});
            Names{f}=[NAME EXT];
        end
        index = f;
        go = 0;
    end
end

if index == 1
    Names = []; Paths = []; Files = [];
end
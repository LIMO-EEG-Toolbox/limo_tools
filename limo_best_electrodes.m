function electrode_vector = limo_best_electrodes

% This function finds the electrode with the maximum F value in each subject.
% The function works on files R2.mat, Condition_effect.mat, Continuous.mat. 
% The function returns a map of frequency showing how often an electrode is
% selected across subjects. Values are the number associated to the channel
% in LIMO.data.chanlocs.urchan
%
% Cyril Pernet 15 July 2010
% Code stollen from GAR :-)
% GAR, 15 August 2010: updated description, uigetfile prompt, frequency
% figure
% -----------------------------
%  Copyright (C) LIMO Team 2010

%% get the 'best' electrodes

current_dir = pwd;

[name,path] = uigetfile('*.mat',['Select an expected chanloc file']);
cd(path); load(name); cd(current_dir)

go = 1; index = 1;
while go == 1
    [name,path] = uigetfile('*.mat',['Select a R2, Condition_effect, or Continuous file for subject ',num2str(index),go]);
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

Ns = length(Names);
electrode_vector = NaN(Ns,1);
abs_electrode_vector = NaN(Ns,1);

for i=1:Ns
    cd(Paths{i})
    load(Names{i})
    name = Names{i}(1:end-4);
    load LIMO
    try
        % check type of map: R2, condition_effect or Continuous
        if strcmp(Names{i},'R2.mat') || strcmp(Names{i}(1:end-6),'Condition_effect')
            tmp = eval(name);
            data{i} = squeeze(tmp(:,:,end-1));
        elseif strcmp(Names{i},'Continuous.mat')
            if size(Continuous,3) > 1
                tmp = eval(name);
                data{i} = squeeze(tmp(:,:,end-1));
            else
                if i == 1
                    n = eval(cell2mat(inputdlg('which regressor to use?','several regressors found')));
                end
                tmp = eval(name);
                data{i} = squeeze(tmp(:,:,n,end-1));
            end
        else
            error('limo_best_electrodes: input file not supported')
        end
        data_size(:,:,i) = size(data{i});
        [v,f] = max(data{i},[],2); % max over electrodes
        electrode_vector(i)=find(v == max(v));
        % frames(i) = f(electrode_vector(i));
        try
            electrode_vector(i) = LIMO.data.chanlocs(electrode_vector(i)).urchan;
        catch 
            fprintf('can''t read data.chanlocs.urchan subject %g \n',i)
        end
        fprintf('subject %g analysed \n',i); 
    catch ME
        message = sprintf('file error, the map of subject %s%s is not recognized',Paths{i},Names{i});
        error([message])
    end
end

% usually no output just write it down on the disk for latter use
if nargout == 0
    cd(current_dir)
    name = cell2mat(inputdlg('Save electrode vector as','Name'));
    if isempty(name)
        return
    else
        save ([name],'electrode_vector')
        assignin('base',[name],electrode_vector)
    end
end


%% do the map

if sum(isnan(electrode_vector)) == 0
    
    origin = which('limo_eeg'); origin = origin(1:end-10); cd(origin);
    data = zeros(1,length(expected_chanlocs));
    
    for S=1:Ns
        data(electrode_vector(S)) = data(electrode_vector(S))+1;
    end
    
    % create the frequency map
    [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
        topoplot( data,expected_chanlocs,'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on','conv','on');
    freqmap = grid_or_val(end:-1:1,:); % reverse row order to get back of the head at the bottom of the image
    freqmap(34,67)=NaN;freqmap(67,34)=NaN;freqmap(34,1)=NaN;freqmap(1,34)=NaN;
    freqmap(isnan(freqmap))=max(freqmap(:))+1;
    
    % make the figure
    figure('Color','w','NumberTitle','off','Name','limo_tools: best electrode frequency map')
    imagesc(freqmap,[0 max(freqmap(:))])
    axis tight;axis square;axis off
    cc=colormap(jet);cc(1,:)=[.9 .9 .9];cc(end,:)=[1 1 1];colormap(cc);
    cd(current_dir)
end

% % get the subjects maps
% maxmaps = zeros(67,67,Ns); % 67 = default topoplot map size
% for S=1:Ns 
%     % data has format E x F x subjects with subjects in cell
%     [h grid_or_val plotrad_or_grid, xmesh, ymesh]= ...
%         topoplot( data{S}(:,frames(S)),chanlocs{S},'style','both','electrodes','off','hcolor','none','numcontour',0,'whitebk','on','noplot','on');
%     maxmaps(:,:,S) = grid_or_val(end:-1:1,:); % reverse row order to get back of the head at the bottom of the image
% end
% 
% 
% % create the frequency map
% topomask = zeros(size(maxmaps,1),size(maxmaps,2));
% topomask = topomask.*mean(maxmaps,3); % using mean we keep the NaNs
% freqmap  = zeros(size(topomask,1),size(topomask,2),Ns); % default size is 67 in EEGlab
% 
% % take the max and neighbours
% for S=1:Ns
%     tmp = squeeze(maxmaps(:,:,S));
%     [R,C] = find(tmp==max(tmp(:)));
%     freqR2maps(R,C,S) = 1;
% end
%  
% % make the figure
% figure('Color','w','NumberTitle','off','Name','limo_tools: best electrode frequency map')
% cst=4; % constant to add so that the background has one colour and the lowest frequency value another one
% toplot = squeeze(sum(freqmap,3));
% imagesc(toplot+topomask.*cst,[0 max(toplot(:))])
% axis tight;axis square;axis off
% colormap(jet)
% h=colormap;
% colormap(h(end:-1:1,:))
% cd(current_dir)















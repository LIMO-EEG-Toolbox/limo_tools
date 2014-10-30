function out = limo_struct2mat(in)

% routine to read EEGLAB structures files .daterp .datersp .icaerp etc ..
% and convert the cell arraw of channels or components into a matrix
%
% FORMAT Matrix = limo_struct2mat(Structure)
%
% INPUT Structure is a strcture containing cells of channels or components
%       such files are usually called e.g. xxxx.icaerp xxxx.datersp
% OUTPUT Matrix is a matrix of data
%
% Cyril Pernet Novembre 2014
% ----------------------------------
%  Copyright (C) LIMO Team 2014

F = fieldnames(in);
% how many components or channels
n = 0;
for f=1:size(F,1)
    if strncmp(cell2mat(F(f)),'comp',4) || strncmp(cell2mat(F(f)),'chan',4)
        A = strfind(cell2mat(F(f)),'base');
        B = strfind(cell2mat(F(f)),'boot');
        C = strfind(cell2mat(F(f)),'label');
        if isempty(A) && isempty(B) && isempty(C)
            if n==0
                tmp_data = getfield(in,cell2mat(F(f)));
            end
            n = n +1;
        end
    end
end

% other dim
if numel(size(tmp_data)) == 2
    [frames,trials]=size(tmp_data);
elseif numel(size(tmp_data)) == 3;
    [freq,time,trials]=size(tmp_data);
end

% build the matrix
if numel(size(tmp_data)) == 2 % erp spec itc
    
    out=NaN(n,frames,trials);
    out(1,:,:) = tmp_data; index = 2;
    for f=2:size(F,1)
        if strncmp(cell2mat(F(f)),'comp',4) || strncmp(cell2mat(F(f)),'chan',4)
            A = strfind(cell2mat(F(f)),'base');
            B = strfind(cell2mat(F(f)),'boot');
            C = strfind(cell2mat(F(f)),'label');
            if isempty(A) && isempty(B) && isempty(C)
                
                out(index,:,:) = getfield(in,cell2mat(F(f)));
                index = index +1;
            end
        end
    end
    
elseif numel(size(tmp_data)) == 3; % ersp
    
    out=NaN(n,freq,time,trials);
    out(1,:,:,:) = tmp_data; index = 2;
    for f=2:size(F,1)
        if strncmp(cell2mat(F(f)),'comp',4) || strncmp(cell2mat(F(f)),'chan',4)
            A = strfind(cell2mat(F(f)),'base');
            B = strfind(cell2mat(F(f)),'boot');
            C = strfind(cell2mat(F(f)),'label');
            if isempty(A) && isempty(B) && isempty(C)
                out(index,:,:,:) = getfield(in,cell2mat(F(f)));
                index = index +1;
            end
        end
    end
end






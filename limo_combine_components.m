function combined_ica = limo_combine_components(varargin)

% compute a weighted average of components
%
% FORMAT combined_ica = limo_combine_components(data,weights)
%
% INPUT data a 3D matrix representing set of components 
%            dim compoments*time frames*epochs
%       weights are the rows from EEG.icawinv
%       method 'sum'
%
% OUTPUT combined_ica is the weighted average of comnponents
%
% Cyril Pernet & Ramon Martinez-Cancino 23-10-2014 updates for components (ICA)
%% ----------------------------------
%  Copyright (C) LIMO Team 2010

%% check data orientation
data = varargin{1};
weights = varargin{2};
if nargin == 2
    method = 'sum';
elseif nargin == 3
    method = varargin{3};
end
clear varargin
nb_comp = size(data,1);
[n,p]=size(weights);
if p==nb_comp
    weights = weights';
    [n,p]=size(weights);
end

%% down to business
switch method
    
    %------
    case {'sum'}

        W = sqrt(sum(weights.^2,2));
end

combined_ica = mean(data.*repmat(W,[1,size(data,2),size(data,3)]),1);

% figure;plot(mean(data,3)','--','LineWidth',1.5)
% hold on; plot(mean(combined_ica,3),'k','LineWidth',3); axis tight

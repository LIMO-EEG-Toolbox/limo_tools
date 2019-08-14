function combined_ica = limo_combine_components(varargin)

% compute a weighted average of components
%
% FORMAT combined_ica = limo_combine_components(data,weights,invweights,whichic)
%
% INPUT data a 3D matrix representing set of components 
%            dim compoments*time frames*epochs
%       weights are the columns EEG.weights
%       invweights are the columns EEG.icawinv
%       whichic are the com ponents to 'combine'
%       method 'sum' 'max' 'maxvar' (default)
%
% OUTPUT combined_ica is the weighted average of comnponents
%
% Cyril Pernet & Ramon Martinez-Cancino 23-10-2014 updates for components (ICA)
%% ----------------------------------
%  Copyright (C) LIMO Team 2010

%% check data orientation
data       = varargin{1};
weights    = varargin{2};
invweights = varargin{3};
whichic    = varargin{4};
if nargin == 4
    method = 'maxvar';
elseif nargin == 5
    method = varargin{5};
end
clear varargin

%% down to business
switch method
    
    %------
    case {'sum'} 
        % Average power of ica weights on eletrodes
        W = sqrt(sum((invweights(:,whichic)).^2,1)); 
        combined_ica = mean(data(whichic,:,:).*repmat(W',[1,size(data,2),size(data,3)]),1);

    case ('max') 
        % Takes the ica with max activation variance
        meanicaact = mean(data,3);
        [tmp, ind] = max(var(meanicaact(whichic,:), [], 2));      
        combined_ica = data(whichic(ind),:,:);
        
    case ('maxvar') 
        % Takes the ica which has the max mean variance accounted for
        meanicaact = mean(data,3);
        for iComp = 1:length(whichic)
            comp = whichic(iComp);
            dataMinusIca = invweights*meanicaact - invweights(:,comp)*meanicaact(comp,:);
            varica(iComp) = mean(var(dataMinusIca, [], 2));
        end
        [tmp, ind] = max(varica);
        combined_ica = data(whichic(ind),:,:);   
        
end
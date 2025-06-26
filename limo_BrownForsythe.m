function dataout = limo_BrownForsythe(Y1r,Y2r)

% computes the Brown-Forsythe test for variance equality
% the analysis is restricted to two EEG data sets running
% a GLM on z transormed data with z = abs(xij-median(xj))
%
% FORMAT dataout = limo_BrownForsythe(data1,data2,alphavalue)
% INPUT data1 and data 2 are EEG matrices [channels * freq/time frames * trials]
% OUTPUT dataout is a matrix of results [channels * frames * F/p]
%
% Cyril Pernet 2020
% ------------------------------
%  Copyright (C) LIMO Team 2020

%% check data
if ischar(Y1r)
    Y1r = load(Y1r); Y1r = Y1r.(cell2mat(fieldnames(Y1r)));
end

if ischar(Y2r)
    Y2r = load(Y2r); Y2r = Y2r.(cell2mat(fieldnames(Y2r)));
end

istf = 0;
if ndims(Y1r)== 4 && ndims(Y2r) == 4
    Y1r = limo_tf_4d_reshape(Y1r);
    Y2r = limo_tf_4d_reshape(Y2r);
    istf = 1;
end

if ndims(Y1r)~=3 && ndims(Y2r) ~=3
    error('data must be 4D/3D matrices')
end

if numel(Y1r) ~= numel(Y2r)
   error('data in do not have the same dimension - use NaN for missing values') 
end

dataout = NaN(size(Y1r,1),size(Y1r,2),2);

%% compute
array = intersect(find(~isnan(Y1r(:,1,1))),find(~isnan(Y2r(:,1,1))));  % get common channels
for channel = length(array):-1:1
    % z transform using the median
    Z1                   = squeeze(Y1r(channel,:,:));
    Z1                   = abs(Z1 - repmat(median(Z1,2),[1 size(Z1,2)]));
    Z1(Z1==0)            = min(Z1(:)); % if exactly 0 use min
    Z2                   = squeeze(Y2r(channel,:,:));
    Z2                   = abs(Z2 - repmat(median(Z2,2),[1 size(Z2,2)]));
    Z2(Z2==0)            = min(Z2(:)); 
    % GLM
    Y                    = [Z1 Z2]';
    X                    = [[ones(size(Z1,2),1);zeros(size(Z2,2),1)], [zeros(size(Z1,2),1);ones(size(Z2,2),1)],ones(size(Y,1),1)];
    T                    = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total 
    R                    = eye(size(Y,1)) - X*pinv(X);                                        
    E                    = Y'*R*Y;                                                            % SS Error
    df                   = rank(X)-1;
    dfe                  = size(Y,1)-rank(X);
    Betas                = pinv(X)*Y;
    C                    = eye(size(X,2));
    C(:,size(X,2))       = 0;                              
    C0                   = eye(size(X,2)) - C*pinv(C);     
    X0                   = X*C0;                          
    R0                   = eye(size(Y,1)) - (X0*pinv(X0)); 
    M                    = R0 - R;                         
    H                    = (Betas'*X'*M*X*Betas);                                             % SS Effects 
    dataout(channel,:,1) = (diag(H)./df) ./ (diag(E)/dfe);
    dataout(channel,:,2) = 1 - fcdf(squeeze(dataout(channel,:,1)), df, dfe);
end

% matlab version - gives almost identical results (precision and case value=0)
% for channel = length(array):-1:1
%     for frame = size(dataout,2):-1:1
%         Y = [squeeze(Y1r(channel,frame,:)) squeeze(Y2r(channel,frame,:))];
%         [testout(channel,frame,2),stats] = vartestn(Y,'TestType','BrownForsythe','Display','off');
%         testout(channel,frame,1) = stats.fstat;
%     end
% end



%% output
if istf == 1
    dataout = limo_tf_4d_reshape(dataout);
end

if nargout == 0
   save('BrownForsythe_test.mat','dataout','-v7.3') 
end

hfig = figure;
mask = (squeeze(dataout(:,:,2)) <= 0.05);
toplot = squeeze(dataout(:,:,1));
imagesc(toplot.*mask);
xlabel('frames'); ylabel('channels')
toplot(mask==0)=NaN;
colormap(limo_color_images(toplot))
title('Variance Homogeneity test')
drawnow; saveas(hfig, 'BrownForsythe','fig'); close(hfig)





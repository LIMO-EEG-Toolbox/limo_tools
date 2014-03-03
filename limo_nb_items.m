function nb_items = limo_nb_items(LIMO)


% FORMAT limo_nb_items(LIMO)
% INPUT  LIMO the LIMO structure
% OUTPUT nb_items the number of trials per condition
%
% silly routine to return the number of trials used in each column of X
% Cyril Pernet 01/02/2012  updated for all regressors Jan 2014
% ------------------------------------------------------------
%  Copyright (C) LIMO Team 2014

X = LIMO.design.X;
P = sum(LIMO.design.nb_conditions)+sum(LIMO.design.nb_interactions);
nb_items = [sum(X(:,1:P)) repmat(size(X,1),1,LIMO.design.nb_continuous)];


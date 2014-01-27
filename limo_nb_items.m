function [nb_items] = limo_nb_items(LIMO)

% silly routine to return the number of trials used in each column
% Cyril Pernet 01/02/2012
% -----------------------------------------------------
%  Copyright (C) LIMO Team 2010

X = LIMO.design.X;
nb_continuous = LIMO.design.nb_continuous;
P = size(X,2)-1-nb_continuous;
nb_items = sum(X(:,1:P));

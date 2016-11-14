function [est,HDI,Outliers,KDE]=data_plot(Data,varargin)

% plots the data split by groups showing each data points with the
% distribution and a summary statistics estimator with 95&% Bayes boot HDI
%
% FORMAT: [est,HDI]= rst_data_plot(Data,options)
%         [est,HDI]= rst_data_plot(Data,'between',0.25,'within',0.025,'pointsize',50,'estimator','median','coverage',0.95)
%
% INPUT: Data is a matrix, data are plotted colmun-wise
%        options are
%                'between' for the distance between distributions
%                'within' for the distance/width between points in a group
%                'point_size' the size of data points
%                'estimator' can be 'median', 'mean', 'trimmed mean'
%                'coverage' indicate the probability coverage of the HDI
% OUTPUT: est is the summary statistics
%         HDI the 95% High Density Interval (Bayes bootstrap)
%         Ouliers is a binary matrix indicating S-outliers
%         KDE is the hernel density estimated (returns bins and frequencies)
%
% example: G = gamrnd(3,2,100,1); G = zscore(G);
%          B = betarnd(0.5,0.5,100,1).*3; B = zscore(B);
%          Data = [zscore(randn(100,1)) G B]; 
%          Data(randi(100,10,1)) = NaN; % cool it takes missing data as well
%          [est,HDI,Outliers,KDE]=data_plot(Data,'between',0.25,'within',0.025,'pointsize',50,'estimator','mean')
%
% -------------------------------------------------------------------------
%     Copyright (C) <2016>  <Cyril Pernet>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
   
disp(' ')
disp('Copyright (C) 2016  Cyril Pernet')
disp('Data_plot comes with ABSOLUTELY NO WARRANTY.')
disp('This is free software, and you are welcome to redistribute it under')
disp('certain conditions; see http://www.gnu.org/licenses/ for details')
disp(' ')
   
%% Defaults

% hard coded
Nb = 1000;      % number of bootstrap samples, no need more really
trimming = 20;  % 20% trimming is the standard approach
decile = .5;    % Median is estimated using the 5th decile of Harell Davis

% soft coded, see options
between_gp_dispersion = 0.25;
within_gp_dispersion = 0.025;
point_size = 50;
prob_coverage = 95/100; % prob coverage of the HDI

% check inputs
if ~exist('Data','var')
    [name,place,sts]=uigetfile({'*.csv;*.tsv;*.txt'},'Select data file');
    if sts == 0
        return
    else
        if strcmp(name(end-2:end),'csv')
            Data = csvread([place name]); % expect a comma between variables
        elseif strcmp(name(end-2:end),'tsv')
            Data = dlmread([place name]); % expect a tab between variables
        elseif strcmp(name(end-2:end),'txt')
            Data = load([place name]); % assumes just a space between variables
        end
    end
end

if nargin <= 1
    estimator = questdlg('Which summary statistics to plot?','Stat question','Mean',' 20% Trimmed mean','Median','Median');
    if strcmp(estimator,' 20% Trimmed mean')
        estimator = 'Trimmed mean';
    end
else
    for n=1:(nargin-1)
        if strcmpi(varargin(n),'between')
            between_gp_dispersion = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'within')
            within_gp_dispersion = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'point_size')
            point_size = cell2mat(varargin(n+1));
        elseif strcmpi(varargin(n),'estimator')
            if ~strcmpi(varargin(n+1),'median') && ...
                    ~strcmpi(varargin(n+1),'mean') && ...
                    ~strcmpi(varargin(n+1),'trimmed mean')
                error(['estimator ' cell2mat(varargin(n+1)) ' is not recognized'])
            else
                estimator = cell2mat(varargin(n+1));
            end
        elseif strcmpi(varargin(n),'coverage')
            prob_coverage = cell2mat(varargin(n+1));
        end
    end
end

%% how many groups
grouping = size(Data,2);
gp_index = 1:(1+between_gp_dispersion):(grouping*(1+between_gp_dispersion));
KDE = cell(1,grouping);

%% Summary statistics
if strcmpi(estimator,'Mean')
    est = nanmean(Data,1);
elseif strcmpi(estimator,'Trimmed mean')
    est = trimmean(Data,trimming);
elseif strcmpi(estimator,'Median')
    est = harelldavis(Data,decile); % Median estimated using the 5th decile of Harell Davis
end

%% Compute High Density Intervals (Bayesian Bootstrap)
% sample with replcaement from Dirichlet
% sampling = number of observations, e.g. participants
n=size(Data,1);
bb = zeros(Nb,grouping);
for boot=1:Nb % bootstrap loop
    theta = exprnd(1,[n,1]);
    weigths = theta ./ repmat(sum(theta),n,1);
    resample = (datasample(Data,n,'Replace',true,'Weights',weigths));

    % compute the estimator
    if strcmpi(estimator,'Mean')
        bb(boot,:) = nanmean(resample,1);
    elseif strcmpi(estimator,'Trimmed Mean')
        bb(boot,:) = trimmean(resample,20);
    elseif strcmpi(estimator,'Median')
        bb(boot,:) = harelldavis(resample,.5);
    end
end

sorted_data = sort(bb,1); % sort bootstrap estimates
upper_centile = floor(prob_coverage*size(sorted_data,1)); % upper bound
nCIs = size(sorted_data,1) - upper_centile;
HDI = zeros(2,grouping);

%% outliers
Outliers = S_outliers(Data);

%% start
figure; hold on

% select color scheme
color_scheme = select_colours(grouping);

% iterate per group/condition
for u=1:grouping
    % remove NaN
    [tmp,sort_index] = sort(Data(~isnan(Data(:,u)),u));

    % find outliers removing NaN
    outliers = Outliers(~isnan(Data(:,u)),u);
    outliers = find(outliers(sort_index));

    %% Scatter plot
    % creater a matrix with spread = 2
    change = find(diff(tmp) < 0.1);
    if isempty(change)
        Y = tmp;
    else
        Y = NaN(length(tmp),2);
        Y(1:change(1),1) = tmp(1:change(1));
        c_index = 2;
        for c=2:length(change)
            Y((change(c-1)+1):change(c),c_index) = tmp((change(c-1)+1):change(c));
            if mod(c,2) == 0
                c_index = 1;
            else
                c_index = 2;
            end
        end
        Y((change(c)+1):length(tmp),c_index) = tmp((change(c)+1):length(tmp));
    end
    X = repmat([gp_index(u)-(within_gp_dispersion/2) gp_index(u)+(within_gp_dispersion/2)],[length(tmp),1]);
    X(isnan(Y)) = NaN;
    for p=1:size(Y,2)
        scatter(X(:,p),Y(:,p),point_size,color_scheme(u,:));
        scatter(X(outliers,p),Y(outliers,p),point_size,color_scheme(u,:),'filled');
    end

    %% Get the kernel density estimate
    [bc,K]=RASH(tmp,100);
    KDE{u} = [bc',K'];
    % remove 0s
    bc(K==0)=[]; K(K==0)= [];
    % create symmetric values
    K = (K - min(K)) ./ max(K); % normalize to [0 1] interval
    high=(K/2); low=(-high);

    % plot contours
    y1 = plot(high+gp_index(u),bc); set(y1,'Color',color_scheme(u,:)); hold on
    y2 = plot(low+gp_index(u),bc); set(y2,'Color',color_scheme(u,:));
    if isnumeric(y1)
        y1 = get(y1); y2 = get(y2); % old fashion matlab
    end
    % fill
    xpoints=[y2.XData',y1.XData']; filled=[y2.YData',y1.YData'];
    % check that we have continuous data, otherwise 0 pad
    if diff(xpoints(1,:)) > 0.001*(range(xpoints(:,1)))
        xtmp = NaN(size(xpoints,1)+1,size(xpoints,2));
        xtmp(1,:) = [gp_index(u) gp_index(u)];
        xtmp(2:end,:) = xpoints;
        xpoints = xtmp; clear xtmp
        ytmp = NaN(size(filled,1)+1,size(filled,2));
        ytmp(1,:) = filled(1,:);
        ytmp(2:end,:) = filled;
        filled = ytmp; clear ytmp
    end

    if diff(xpoints(end,:)) > 0.001*(range(xpoints(:,1)))
        xpoints(end+1,:) = [gp_index(u) gp_index(u)];
        filled(end+1,:)  = filled(end,:);
    end
    hold on; fillhandle=fill(xpoints,filled,color_scheme(u,:));
    set(fillhandle,'LineWidth',2,'EdgeColor',color_scheme(u,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color

    %% add IQR - using again Harell-Davis Q
    ql = harelldavis(tmp,0.25);
    [~,position] = min(abs(filled(:,1) - ql));
    plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
    if  strcmpi(estimator,'median')
        qm = est(u);
    else
        qm = harelldavis(tmp,0.5);
    end
    [~,position] = min(abs(filled(:,1) - qm));
    plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
    qu = harelldavis(tmp,0.75);
    [~,position] = min(abs(filled(:,1) - qu));
    plot(xpoints(position,:),filled(position,:),'Color',color_scheme(u,:));
    clear xpoints filled

    %% Add the High Density Intervals
    tmp = sorted_data(:,u);
    ci = 1:nCIs; ciWidth = tmp(ci+upper_centile) - tmp(ci); % all centile distances
    [~,index]=find(ciWidth == min(ciWidth)); % densest centile
    if length(index) > 1; index = index(1); end % many similar values
    HDI(1,u) = tmp(index);
    HDI(2,u) = tmp(index+upper_centile);

    % plot this with a rectangle
    X = (gp_index(u)-0.3):0.1:(gp_index(u)+0.3);
    plot(X(3:5),repmat(est(u),[1,length(X(3:5))]),'LineWidth',6,'Color',color_scheme(u,:));
    rectangle('Position',[X(3),HDI(1,u),X(5)-X(3),HDI(2,u)-HDI(1,u)],'Curvature',[0.4 0.4],'LineWidth',3,'EdgeColor',color_scheme(u,:))
end

%% finish off the figure
cst = max(abs(diff(Data(:)))) * 0.1;
axis([0.3 gp_index(grouping)+0.7 min(Data(:))-cst max(Data(:))+cst])

if size(Data,1) == 1
    title(sprintf('Data distribution with %s and 95%% Highest Density Interval',estimator),'FontSize',16);
else
    title(sprintf('Data distributions with %ss and 95%% Highest Density Intervals',estimator),'FontSize',16);
end
grid on ; box on; drawnow

% add an output if not specified during the call
if nargout == 0
    S = questdlg('Save computed data?','Save option','Yes','No','No');
    if strcmp(S,'Yes')
        if exist(place,'var')
            place = pwd;
        end
        cd(place); tmp = [est; HDI];
        if strcmpi(estimator,'Mean')
            save('Mean_and_HDI.txt','tmp','-ascii');
        elseif strcmpi(estimator,'Trimmed Mean')
            save('Trimmed-Mean_and_HDI.txt','tmp','-ascii');
        elseif strcmpi(estimator,'Median')
            save('Median_and_HDI.txt','tmp','-ascii');
        end
        save('S-outliers.txt','Outliers','-ascii');
        fprintf('Data saved in %s\n',place)
    end
end

end

% ---------------
%% Subfunctions
% --------------
function TM = trimmean(varargin)

% the trimmed mean is the mean of the data excluding the highest and lowest
% K data values (K=round(N*percent)) and N is the number of values
% in data) - see e.g. Wilcox, Rand R. "Introduction to Robust Estimation
% and Hypothesis Testing." New York: Academic Press. 2005.

data = varargin{1};
percent = 20/100; % default

if nargin>3
    error('too many arguments')
elseif nargin == 2
    percent = varargin{2};
    if ~isnumeric(percent)
        error('percent of trimming must be a numeric')
    elseif percent > 1
        percent = percent / 100;
    end
end

[n,p]=size(data);
if n== 1 && p>2
    data = data'; % transpose row vector in column vector if needed
    [n,p]=size(data);
end
TM = NaN(1,p);

% compute
for c=1:p
    tmp = data(:,c);
    tmp(isnan(tmp)) = [];
    N = size(tmp,1);
    K=round(N*percent);
    if K == 0
        error('not enough data points to trim')
    end
    X = sort(tmp,1);
    TM(c) = mean(X((K+1):(N-K),:),1);
end

end % closes the trimmean sub-function

function HD = harelldavis(varargin)

% FRANK E. HARRELL and C. E. DAVIS (1982).
% A new distribution-free quantile estimator
% Biometrika 69 (3): 635-640. doi: 10.1093/biomet/69.3.635

% defaults
X = varargin{1};
[p,N]=size(X); % number of estimates to compute
if p==1 % if X row vector, transpose to column
    X=X';
    [p,N]=size(X);
end

q=.5; % median
if nargin > 2
    disp('too many arguments in, only the 2 first ones will be used')
elseif nargin == 2
    q = varargin{2};
end
HD = NaN(1,N);

% compute
for i=1:N
    x = X(~isnan(X(:,i)),i);
    n=length(x);
    m1=(n+1).*q;
    m2=(n+1).*(1-q);
    vec=1:n;
    w=betacdf(vec./n,m1,m2)-betacdf((vec-1)./n,m1,m2);
    y=sort(x);
    HD(i)=sum(w(:).*y(:));
    clear x n m1 m2 vec w y
end

end % closes the harelldavis subfunction

function colours = select_colours(n)

% this is a 124 units semi-continous color scale generated with cubehelix
% Green, D. A., 2011, `A colour scheme for the display of
% astronomical intensity images', Bulletin of the Astronomical
% Society of India, 39, 289.
% The map is available here https://figshare.com/articles/Brain_Colour_Scales/3574173/1

R = [0.15945	0.16168	0.16302	0.16349	0.16311	0.1619	0.15989	0.1571	0.15358	0.14937	0.1445	0.13904	0.13303	0.12652	0.11958	0.11227	0.10464	0.09677	0.088724	0.08057	0.072379	0.064223	0.056172	0.0483	0.040678	0.033377	0.026469	0.020021	0.014102	0.0087784	0.0041128	0.00016682	0	0	0	0	0	0	0	0.00079018	0.0055218	0.011388	0.018407	0.02659	0.035945	0.046473	0.058173	0.071036	0.085048	0.10019	0.11644	0.13378	0.15216	0.17155	0.19191	0.21319	0.23535	0.25832	0.28205	0.30648	0.33155	0.35719	0.38333	0.4099	0.43682	0.46403	0.49144	0.51897	0.54656	0.57413	0.60159	0.62887	0.6559	0.6826	0.70889	0.73472	0.76	0.78468	0.80868	0.83196	0.85445	0.87611	0.89687	0.91671	0.93557	0.95342	0.97023	0.98597	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	0.99531	0.99088];
G = [0	0	0	0	0.0071101	0.02196	0.037099	0.052499	0.068131	0.083964	0.099964	0.1161	0.13233	0.14863	0.16495	0.18126	0.19752	0.2137	0.22976	0.24565	0.26135	0.27681	0.292	0.30689	0.32144	0.33562	0.3494	0.36274	0.37563	0.38804	0.39994	0.41131	0.42213	0.43239	0.44206	0.45115	0.45963	0.4675	0.47477	0.48141	0.48745	0.49288	0.4977	0.50193	0.50558	0.50867	0.5112	0.51321	0.51471	0.51572	0.51628	0.5164	0.51613	0.5155	0.51453	0.51326	0.51173	0.50998	0.50805	0.50597	0.50378	0.50154	0.49927	0.49702	0.49483	0.49275	0.4908	0.48903	0.48749	0.4862	0.48521	0.48455	0.48425	0.48435	0.48487	0.48586	0.48732	0.48929	0.49178	0.49483	0.49844	0.50263	0.50741	0.51279	0.51878	0.52539	0.53261	0.54044	0.54888	0.55792	0.56755	0.57776	0.58854	0.59987	0.61173	0.6241	0.63695	0.65027	0.66401	0.67816	0.69268	0.70754	0.72271	0.73815	0.75382	0.7697	0.78574	0.8019	0.81816	0.83446	0.85078	0.86707	0.8833	0.89944	0.91544	0.93128	0.94692	0.96233	0.97748	0.99235	1	1	1	1];
B = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0.016251	0.038138	0.060952	0.084641	0.10915	0.13442	0.16038	0.18698	0.21413	0.24177	0.26982	0.2982	0.32683	0.35563	0.38452	0.41342	0.44223	0.47088	0.49929	0.52736	0.55503	0.5822	0.6088	0.63475	0.65999	0.68443	0.70802	0.73068	0.75236	0.77301	0.79257	0.81099	0.82824	0.84427	0.85906	0.87257	0.88479	0.8957	0.9053	0.91357	0.92052	0.92616	0.9305	0.93355	0.93535	0.93592	0.9353	0.93352	0.93063	0.92668	0.92172	0.91582	0.90902	0.90141	0.89303	0.88398	0.87432	0.86412	0.85347	0.84245	0.83113	0.81961	0.80795	0.79625	0.78458	0.77303	0.76168	0.7506	0.73988	0.72958	0.71978	0.71055	0.70195	0.69405	0.68691	0.68057	0.6751	0.67053	0.66691	0.66427	0.66265	0.66207	0.66256	0.66411	0.66676	0.6705	0.67533	0.68124	0.68823	0.69627	0.70534	0.71542	0.72648	0.73848	0.75137	0.76511];

sampling = round(1:(124/n-1):124);
colours = [R(sampling)' G(sampling)' B(sampling)'];

end

function class = S_outliers(data)

% This methods allows to detect outliers in the data using a S estimator
% Like the MAD, the S estimator is robust and relies on the median of
% absolute distances. However, MAD relies on the distance to the median
% whereas S relies on the median of absolute distances. Because it doesn't
% rely on an estimator of central tendency (like the MAD) it works well for
% non symmetric distributons
%
% Ref. Rousseeuw, P.J. and Croux C. (1993). Alternatives to the the median
% absolute deviation. Journal of the American Statistical Association, 88
% (424) p 1273-1263

k = 2.2414; % = sqrt(chi2inv(0.975,1))
[n,p]=size(data);
class = NaN(n,p);
distance = NaN(n,p);

for p=1:size(data,2)
    tmp = data(:,p);
    points = find(~isnan(tmp));
    tmp(isnan(tmp)) = [];

    % count all distances
    n = length(tmp);
    for i=1:n
        j = points(i);
        indices = 1:n; indices(i) = [];
        distance(j,p) = median(abs(tmp(i) - tmp(indices)));
    end

    % get the S estimator
    % consistency factor c = 1.1926;
    Sn = 1.1926*median(distance(points,p));

    % get the outliers in a normal distribution
    class(:,p) = (distance(:,p) ./ Sn) > k; % no scaling needed as S estimates already std(data)
    class(:,p) = class(:,p)+isnan(data(:,p));
end

end % closes the S_outliers subfunction

function [bc,K]= RASH(varargin)

% Computes the Random Average Shifted Histogram
% the algorithm is coded based on Bourel et al. (2014)
% Computational Statistics and Data Analysis 79

data = varargin{1};
if size(data,1) == 1
    data = data';
end
n = length(data);
m = 100; % this is the parameter setting the number of hist in ASH and RASH

if nargin == 2
    m = varargin{2};
end
clear varargin

h = 2.15*sqrt(var(data))*n^(-1/5);
delta = h/m;

% 1st make a mesh with size delta
t0 = min(data)-min(diff(data))/2;
tf = max(data)+min(diff(data))/2;
nbin = ceil((tf-t0)/delta);
binedge = t0:delta:(t0+delta*nbin);
out = find(binedge>tf);
if out == 1
    binedge(out) = tf;
else
    binedge(out(1)) = tf;
    binedge(out(2:end)) = [];
end

% 2nd Get the weight vector.
kern = inline('(15/16)*(1-x.^2).^2');
ind = (1-m):(m-1);
den = sum(kern(ind/m));% Get the denominator.
wm = m*(kern(ind/m))/den;% Create the weight vector.

% 3rd compute bin with shifted edges
RH=zeros(1,nbin);
RSH=zeros(m,nbin);
for e=1:m
    v = binedge + (delta*randn(1,1)); % e is taken from N(0,h);
    v(v<t0) = t0; % lower bound
    v(v>tf) = tf; % upper bound
    nu = histc(data,v);
    nu = [zeros(1,m-1) nu' zeros(1,m-1)];
    for k=1:nbin
        ind=k:(2*m+k-2);
        RH(k)=sum(wm.*nu(ind));
    end
    RSH(e,:) = RH/(n*h);
end
K = mean(RSH,1);
bc = t0+((1:nbin)-0.5)*delta;

end % closes the RASH subfunction

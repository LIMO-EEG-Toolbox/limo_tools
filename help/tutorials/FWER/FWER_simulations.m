
% famillywise error has to do with the number of tests performed and
% assuming tests are independents from each other FWER = 1 ? (1 ? ?)^n 
% so for ?=5/100 if we do 2 tests we should get about 1-(1-5/100)^2 ~ 0.09
% Bonferroni is conservative if tests are correlated

%% illustration of the effect of number of tests
clear
N= 100;
nb_var = 5;
for MC = 1:1000
    r = randn(N,nb_var);
    x1(:,MC) = r(:,1);
    x2(:,MC) = r(:,2);
    x3(:,MC) = r(:,3);
    x4(:,MC) = r(:,4);
    x5(:,MC) = r(:,5);
end

H(1) = mean(ttest(x1));
H(2) = mean(ttest(x2));
H(3) = mean(ttest(x3));
H(4) = mean(ttest(x4));
H(5) = mean(ttest(x5));
FH(1)= mean(ttest(x1)+ttest(x2));
FH(2)= mean(ttest(x1)+ttest(x2)+ttest(x3));
FH(3)= mean(ttest(x1)+ttest(x2)+ttest(x3)+ttest(x4));
FH(4)= mean(ttest(x1)+ttest(x2)+ttest(x3)+ttest(x4)+ttest(x5));

figure('Name','type 1 error under H0: 1000 MC');
subplot(1,2,1); bar(H); set(gca,'XTickLabel',{'1','2','3','4','5'});
hold on; plot([0:6],repmat(0.05,7,1),'--r','LineWidth',2); axis([0 6 0 0.07])
grid on; title('type one error rate per variable','Fontsize',14)

subplot(1,2,2)
bar(FH); hold on; set(gca,'XTickLabel',{'12','123','1234','1234','12345'});
for i=1:4; FWER = 1 - (1 - 5/100)^(i+1);  plot([0:i+1],repmat(FWER,i+2,1),'--r','LineWidth',2); end
grid on; title('Famillywise error rate','Fontsize',14); axis([0 5 0 0.25])
    
%% 

clear all
sample_sizes = [10 50 100 500];
shape = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];

% compare type 1 error for 2 variables from a multivariate normal
% change there corr value
for MC = 1:1000
    for ss =1:4
        for s = 1:11
            SIGMA = [1 shape(s) ; shape(s) 1];
            r = mvnrnd([0 0],SIGMA,sample_sizes(ss));
            X = r(:,1); Y = r(:,2);
            clear r ;
            
            [H(s,ss,1,MC),p(s,ss,1,MC),ci,stats] = ttest(X);
            T(s,ss,1,MC) = stats.tstat;
            [H(s,ss,2,MC),p(s,ss,2,MC),ci,stats] = ttest(Y);
            T(s,ss,2,MC) = stats.tstat;
            FH(s,ss,MC) = (H(s,ss,1,MC)+H(s,ss,2,MC));
        end
    end
end
    
% correlation between t-tests
% shows that the more variables are correlated, the more the tests are correlated
for ss=1:4
    for s=1:11
        r(s,ss) = corr(squeeze(T(s,ss,1,:)),squeeze(T(s,ss,2,:)));
    end
end
figure('Name','correlation between t-tests for variable X and Y')
plot(shape,r,'LineWidth',3); grid on; axis([0 1 0 1]);
xlabel('correlation between variables','Fontsize',14); ylabel('correlation between tests','Fontsize',14)
for f=1:4; R(f) = corr(shape',r(:,f)); end
legend(['N=10 corr=' num2str(R(1))],['N=50 corr=' num2str(R(2))],['N=100 corr=' num2str(R(3))],['N=500 corr=' num2str(R(4))]);
title('correlation between t-tests for variable X and Y as a function of the corr(X,Y)','Fontsize',16)

% illustrate for N=10 and corr between variables X and Y 0 0.5 1
ss = 1; 
for s=1:5:11
   Pearson(squeeze(T(s,ss,1,:)),squeeze(T(s,ss,2,:))); % need to robust correlation toolbox
end

figure('Name','type 1 error under H0: 1000 MC');
index = 1; for f=1:4
tmp = squeeze(H(:,f,:,:)); subplot(2,4,index); plot(shape,mean(tmp,3),'LineWidth',3); grid on
axis([0 1 0.03 0.07]); legend('X','Y'); xlabel('correlation between variables')
hold on; plot([0:5],repmat(.05,1,6),'--r','LineWidth',3);
ylabel('type 1 error rate'); title('error rate per variable','Fontsize',16);
subplot(2,4,index+1); tmp = squeeze(FH(:,f,:)); plot(shape,mean(tmp,2),'LineWidth',3); grid on
axis tight; xlabel('correlation between variables'); ylabel('famillywise type 1 error rate'); 
hold on; plot([0:3],repmat(.09,1,4),'--r','LineWidth',3);
mytitle = sprintf('family wise error rate for N=%g',sample_sizes(f)); title(mytitle,'Fontsize',16); 
index = index+2; end

   
%% to like eeg data but random
clear

% say we have completely random data
data= randn(100,100,10);
figure; imagesc(squeeze(data(:,:,1))); 
colormap('gray'); title('100 electrodes, 100 time frames, 10 trials','Fontsize',14)

% do some tests, if alpha 5% per tests, for 10000 tests we should get
% 1-(1-5/100)^1000
% for the whole space - meaning 
for i=1:100
    H(:,i)=(ttest(squeeze(data(:,i,:))'))';
end
figure; subplot(1,2,1);imagesc(H); colormap('gray'); title( [num2str(mean(H(:))) '% significant cells']);

% Bonferroni correction
N = 10000; alphav = 5/100/N;
for i=1:100
    H2(:,i)=(ttest(squeeze(data(:,i,:))',0,'alpha',alphav))';
end
subplot(1,2,2);imagesc(H2); colormap('gray'); title( [num2str(mean(H(:))) '% significant cells']);


figure; subplot(1,2,1);imagesc(H); colormap('gray'); title( [num2str(mean(H(:))) '% significant cells']);
[clusters,num]=bwlabel(H);
clear cluster_size
for i=1:num
    cluster_size(i) = length(find(clusters == i));
    if cluster_size(i)<5
        H(find(clusters == i)) = 0;
    end
end
subplot(1,2,2);hist(cluster_size); title('cluster size distribution')

% smooth data
G = gauss2d(10,4);
for i=1:10
    newdata(:,:,i)= conv2(squeeze(data(:,:,i)),G,'same');
end
figure; imagesc(squeeze(newdata(:,:,1)));

clear H
for i=1:100
    H(:,i)=(ttest(squeeze(newdata(:,i,:))'))';
end
figure; subplot(1,2,1);imagesc(H); colormap('gray'); title( [num2str(mean(H(:))) '% significant cells']);
[clusters,num]=bwlabel(H);
clear cluster_size
for i=1:num
    cluster_size(i) = length(find(clusters == i));
    if cluster_size(i)<10
        H(find(clusters == i)) = 0;
    end
end
subplot(1,2,2);hist(cluster_size); title('cluster size distribution')
figure; imagesc(H); colormap('gray'); title( [num2str(mean(H(:))) '% significant cells k=10']);


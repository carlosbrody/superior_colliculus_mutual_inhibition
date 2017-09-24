clear
clc

aaa=dir('farm_AG_runs*.mat');
for i=1:length(aaa)
    
    n=aaa(i).name;
    load(n);
    a=(sign(dif1)+1)/2;
    b=(sign(dif2)+1)/2;
    vals(i,1)=mean(a(:));
    vals(i,2)=mean(b(:));    
    
end

vals
% load ba
figure;
hold on
% plot(vals(ba==1,1),vals(ba==1,2),'.b')
% plot(vals(ba==2,1),vals(ba==2,2),'.r')
tol=0.08;
min1=0.9-tol;
max1=0.9+tol;
min2=0.7-tol;
max2=0.7+tol;
plot(vals(:,1),vals(:,2),'.b')
plot(0.9,0.7,'.k','MarkerSize',30)
plot(xlim,[min2 min2],'k--')
plot(xlim,[max2 max2],'k--')
plot([min1 min1],ylim,'k--')
plot([max1 max1],ylim,'k--')

goodind=find(vals(:,1)>min1 & vals(:,1)<max1 & vals(:,2)>min2 & vals(:,2)<max2);
save goodind2 goodind

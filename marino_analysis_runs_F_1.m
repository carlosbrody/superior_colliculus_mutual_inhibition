clear
clc

aaa=dir('farm_F_good_runs*.mat');
for i=1:length(aaa)
    
    n=aaa(i).name;
    load(n);
    a=(sign(dif1)+1)/2;
    b=(sign(dif2)+1)/2;
    vals(i,1)=mean(a(:));
    vals(i,2)=mean(b(:));    
    
end

vals
load ba
figure;
hold on
plot(vals(ba==1,1),vals(ba==1,2),'.b')
plot(vals(ba==2,1),vals(ba==2,2),'.r')
plot(0.9,0.7,'.k','MarkerSize',30)
plot(xlim,[.75 .75],'k--')
plot(xlim,[.65 .65],'k--')
plot([.85 .85],ylim,'k--')
plot([.95 .95],ylim,'k--')

goodind=find(vals(:,1)>.85 & vals(:,1)<.95 & vals(:,2)>.65 & vals(:,2)<.75);
save goodind goodind

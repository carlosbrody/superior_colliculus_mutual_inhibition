% function [] = cluster_farms()
clear
clc


%%%%%%%%% BAYES INFORMATION CRITERION TO DETERMINE OPTIMAL NUMBER OF
%%%%%%%%% CLUSTERS IN THE REDUCED RESPONSE SPACE%%%%%%%%%%%%


nclusters=2:4;
mus=[1 1;-1 1;1 -1;-1 -1];
sigmas=repmat(eye(2),[1,1,4])/2.5;
ps=[1 1 1 1];

figure


for iii=1:length(nclusters)
    iii
    truek=nclusters(iii);
    
    mu=mus(1:truek,:);
    sigma=sigmas(:,:,1:truek);
    p=ps(1:truek);
    
    obj = gmdistribution(mu,sigma,p);
    [y,idx] = random(obj,2000);
    
    subplot(2,2,iii)
    hold on
    for i=1:truek
        plot(y(idx==i,1),y(idx==i,2),'.')
    end
    
    
    for jjj=1:length(nclusters)
        
        options = statset('MaxIter',10000,'Display','off');
        clear bicval
        z=1;
        for zzz=1:20
            try
                GM = fitgmdist(y,nclusters(jjj),'Options',options);
                bicval(z)=GM.BIC;
                z=z+1;
            catch
            end
        end
        bic(iii,jjj)=median(bicval);
    end
end

figure
subplot(2,2,1)
plot(nclusters,bic(1,:),'.-','MarkerSize',30)
title('2')
subplot(2,2,2)
plot(nclusters,bic(2,:),'.-','MarkerSize',30)
title('3')
subplot(2,2,3)
plot(nclusters,bic(3,:),'.-','MarkerSize',30)
title('4')



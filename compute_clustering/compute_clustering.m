clear
clc
close all


%%% to begin, run the julia script savematlab.jl to export the file 
%%% MiniOptimizedC17_SVD_response_matrix.jld into the matlab file
%%% MiniOptimizedC17_SVD_response_matrix.mat

load MiniOptimizedC17_SVD_response_matrix

%%% load the names of the args 
load args



para=results(1).params;
costo=[results(1).cost{:}];


ind=find(costo<-0.0002);

response=response(ind,:);
para=para(ind,:);

para_orig=para;



%normalize the parameters by imposing 0 mean and 1 st.dev.
for i=1:12
    p=para(:,i);
    p=p-mean(p);
    p=p/std(p);
    para(:,i)=p;
end



%%%% RUN PCA IN RESPONSE SPACE %%%
[c,s,l] = pca(response);






%%%%%%%%% BAYES INFORMATION CRITERION TO DETERMINE OPTIMAL NUMBER OF
%%%%%%%%% CLUSTERS IN THE REDUCED RESPONSE SPACE%%%%%%%%%%%%

score=s(:,1:2);

nclusters=2:6;
for i=1:length(nclusters)
    
    options = statset('MaxIter',10000,'Display','off');
    
    z=1;
    clear goodgm
    for j=1:20
        try
            GM = fitgmdist(score,nclusters(i),'Options',options);
            goodgm(z)=GM.BIC;
            z=z+1;
        catch
        end
        
        bic(i)=mean(goodgm);
    end
end


figure
plot(nclusters,bic,'b.-');
xlabel('Number of clusters')
ylabel('Bayes Information Criterion')








%%%% RUN K-MEANS WITH 3 GROUPS %%%
ngroups=3;
indici=kmeans(response,ngroups);
% save indici indici
load indici




%%%% PLOT PCA IN RESPONSE SPACE, COLORED BY GROUP IDENTITY %%%
colo='rbkmcy';
figure
hold on
for i=1:ngroups
    ind=find(indici==i);
    for j=1:length(ind)
        plot(s(ind(j),1),s(ind(j),2),['.' colo(i)],'MarkerSize',15)
    end
end
xlabel('First PCA component')
ylabel('Second PCA component')
title('PCA of dynamics space');





%%%% RUN PCA IN PARAMETER SPACE %%%
[c,s,l] = pca(para);



%%%% PLOT PCA IN PARAMETER SPACE, COLORED BY GROUP IDENTITY %%%
colo='rbkmcy';
figure
hold on
for i=1:ngroups
    ind=find(indici==i);
    for j=1:length(ind)
        plot(s(ind(j),1),s(ind(j),2),['.' colo(i)],'MarkerSize',15)
    end
end
xlabel('First PCA component')
ylabel('Second PCA component')
title('PCA of parameter space');



%%%% COMPUTE AND PLOT LDA IN RESPONSE SPACE, COLORED BY GROUP IDENTITY %%%
X=para;
Y=indici';
r=2;
[~,coef] = FDA(X',Y',r);
newcoord=(para*coef)';
figure
hold on
for i=1:3
    plot(newcoord(1,indici==i),newcoord(2,indici==i),['.' colo(i)],'MarkerSize',15);
end
xlabel('First LDA component')
ylabel('Second LDA component')
title('Linear Discriminant projection in parameter space');




%%% HISTOGRAM OF THE PARAMETERS
figure
for i=1:12
    subplot(3,4,i);
    hold on
    [~,x]=hist(para_orig(:,i),10);
    maxu=0;
    for j=1:ngroups
        ind=find(indici==j);
        y=hist(para_orig(ind,i),x);
        y=y/sum(y);
        plot(x,y,colo(j),'LineWidth',2);
        xlim([min(x),max(x)])
        maxu=max(maxu,max(y));
    end
    ylim([0 maxu])
    
    %underscores act weird in matlab labels!
    str=args{i};
    indi=strfind(str,'_');
    str(indi)=32;
    title(str)
    
end



%%% MEANS OF THE PARAMETERS

for i=1:3
    ind=find(indici==i);
    means(:,i)=mean(para_orig(ind,:));
    stderr(:,i)=std(para_orig(ind,:))/sqrt(size(para_orig,1));
end
figure
for i=1:12
    subplot(3,4,i);
    
    hold on
    
    for j=1:3
        bar(j,means(i,j),colo(j))
    end
    errorbar(1:3,means(i,:),stderr(i,:),'.k','LineWidth',2)
    
    %underscores act weird in matlab labels!
    str=args{i};
    indi=strfind(str,'_');
    str(indi)=32;
    title(str)
    
end






%%% WEIGHT OF THE COEFFICIENTS FOR THE PARAMETERS
figure
for i=1:12
    subplot(3,4,i);
    bar(1:2,coef(i,:))
    ylim([-max(abs(coef(:))) max(abs(coef(:)))])
    %underscores act weird in matlab labels!
    str=args{i};
    indi=strfind(str,'_');
    str(indi)=32;
    title(str)
    
end



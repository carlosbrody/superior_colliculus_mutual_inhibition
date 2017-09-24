clear
clc
close all

first_run=0;

if(first_run)
    %%% add .mat extension to all farm files (only run once!)
    a=dir('farm_*_*');
    for i=1:length(a)
        n=a(i).name;
        str=['mv ' n ' ' n '.mat'];
        system(str);
    end


    
    %%% extract and save parameters and cost (only run once!)
    a=dir('farm_AB*.mat');
    
    b=nan(length(a),9);
    costo=nan(length(a),1);
    for i=1:length(a)
        load(a(i).name,'pars','scost','traj');
        b(i,:)=pars;
        costo(i)=scost;
        
        traj=traj(3:end,:);
        tra{i}=traj;        
        
    end
    tra=tra';
    save saved_data b costo tra %parameters are saved in b

end







load saved_data



% figure
% hist(costo,30)
% return

le=cellfun(@length,tra)';
% ind=find(le<1000);

ind=find(costo<0.005);




b=b(ind,:);
le=le(ind);
costo=costo(ind);



b_orig=b;
costo_orig=costo;

%normalize the parameters by imposing 0 mean and 1 st.dev.
for i=1:9
    p=b(:,i);
    p=p-mean(p);
    p=p/std(p);
    b(:,i)=p;
end


%normalize cost between 0 and 1 (and exclude top 5% outliers)
costo=costo-min(costo);
costo=costo/prctile(costo,95);
costo(costo>1)=1;


%run pca
[c,s,l] = pca(b);

% 
% %plot % of variance explained by each component
% figure;
% plot(1:9,[100*l/sum(l)],'.-')
% xlabel('component')
% ylabel('% of variance')




pcadim1=1;
pcadim2=2;



%%% scatter plot of first 2 pca components, colored by cost
%%% where red=low cost ; green=high cost
figure
hold on
for i=1:size(b,1)
    c=costo(i);
    plot(s(i,pcadim1),s(i,pcadim2),'.','Color',[1-c c 0])
end
xlabel('First PCA component')
ylabel('Second PCA component')
title('PCA Scatter Plot with Cost(red=good)');



%%% scatter plot of first 2 pca components, colored by cost
%%% where red=low cost ; green=high cost
figure
hold on
for i=1:size(b,1)
%     c=costo(i);
    c=le(i)/1000;
    plot(s(i,pcadim1),s(i,pcadim2),'.','Color',[1-c c 0])
end
xlabel('First PCA component')
ylabel('Second PCA component')
title('PCA Scatter Plot with length(green=long)');





%%% run kmeans using ngroups (I already saved a "good run" with 5 groups)
ngroups=2;
% indici=kmeans(b,ngroups);save indici indici
load indici


%%% scatter plot of first 2 pca components, colored by group identity
colo='kbcmry';
figure
hold on
for i=1:ngroups
    ind=find(indici==i);
    for j=1:length(ind)
        plot(s(ind(j),pcadim1),s(ind(j),pcadim2),['.' colo(i)])
    end
end
xlabel('First PCA component')
ylabel('Second PCA component')
title('PCA Scatter Plot with K-means Grouping');




%%% histogram of cost, superimposed for the different groups
figure
hold on
[~,x]=hist(costo_orig,30);
for j=1:ngroups
    ind=find(indici==j);
    y=hist(costo_orig(ind),x);
    y=y/sum(y);
    plot(x,y,colo(j),'LineWidth',2);
    xlim([min(x),max(x)])
end
xlabel('cost')
ylabel('fraction')
title('Cost for the K-means groups');




%%% load names of arguments
load farm_AB_0001 args


%%% histogram of the parameters, superimposed for the different groups
figure
for i=1:9
    subplot(3,3,i);
    hold on
    [~,x]=hist(b_orig(:,i),10);
    maxu=0;
    for j=1:ngroups
        ind=find(indici==j);
        y=hist(b_orig(ind,i),x);
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

    

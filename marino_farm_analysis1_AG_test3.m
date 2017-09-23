clear
clc

first_run=1;

if(first_run)
    
    %%% extract and save parameters and cost (only run once!)
    a=dir('farm_AG_0*');
    
    b=nan(length(a),9);
    costo=nan(length(a),1);
    for i=1:length(a)
        load(a(i).name,'-mat','pars','scost');
        b(i,:)=pars;
        costo(i)=scost;
    end
%     
%     
%     return
%     save saved_data b costo %parameters are saved in b

end


% load saved_data

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
pcadim2=4;



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
title('PCA Scatter Plot with Cost');



%%% run kmeans using ngroups (I already saved a "good run" with 5 groups)
ngroups=2;
indici=kmeans(b,ngroups);
% indici=kmeans(b,ngroups);save indici indici
% load indici

inde{1}=find(indici==1);
inde{2}=find(indici==2);


load goodind2
% goodind=1:241;

inde{1}=intersect(inde{1},goodind);
inde{2}=intersect(inde{2},goodind);



% 
% for i=1:length(good)
%     if(~isempty(find(inde{1}==good(i))))
%         ba(i)=1;
%     else
%         ba(i)=2;
%     end
% end
%     
% 
% return
% 
% for i=1:length(good)    
%     in=['farm_F_' sprintf('%04d',good(i))];        
%     load(in,'-mat');
%     
%     
%     out=['farm_F_good_' sprintf('%04d',i)];
%     save(out);
%     
% end
% 
% 
% return


%%% scatter plot of first 2 pca components, colored by group identity
colo='rb';
figure
hold on
for i=1:ngroups
    ind=inde{i};
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
    ind=inde{j};
    y=hist(costo_orig(ind),x);
    y=y/sum(y);
    plot(x,y,colo(j),'LineWidth',2);
    xlim([min(x),max(x)])
end
xlabel('cost')
ylabel('fraction')
title('Cost for the K-means groups');



%%% load names of arguments
load('FarmFields/farm_F_0001','-mat','args')


%%% histogram of the parameters, superimposed for the different groups
figure
for i=1:9
    subplot(3,3,i);
    hold on
    [~,x]=hist(b_orig(:,i),20);
    maxu=0;
    for j=1:ngroups
        ind=inde{j};
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

    

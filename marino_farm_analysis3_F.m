clear
clc

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
    a=dir('farm_F*.mat');
    
    b=nan(length(a),9);
    costo=nan(length(a),1);
    for i=1:length(a)
        load(a(i).name,'pars','scost');
        b(i,:)=pars;
        costo(i)=scost;
    end
    save saved_data b costo %parameters are saved in b

    

    
    %%% extract and save parameters and cost (only run once!)
    a=dir('farm_F*.mat');
    
    for i=1:length(a)
        load(a(i).name,'traj');
        traj=traj(3:end,:);
        tra{i}=traj;
    end
    save saved_traj tra %parameters are saved in b
    
    
end




load saved_data
load saved_traj

return

b_orig=b;
costo_orig=costo;

%normalize the parameters by imposing 0 mean and 1 st.dev.
for i=1:9
    p=b(:,i);
    meanp(i)=mean(p);
    p=p-meanp(i);
    stdp(i)=std(p);
    p=p/stdp(i);
    b(:,i)=p;
end


%normalize cost between 0 and 1 (and exclude top 5% outliers)
costo=costo-min(costo);
costo=costo/prctile(costo,95);
costo(costo>1)=1;


%run pca
[c,s,l] = pca(b);


% s=b*c;

pcadim1=1;
pcadim2=2;



%%% run kmeans using ngroups (I already saved a "good run" with 5 groups)
ngroups=5;
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
xlim1=[-4,4];
ylim1=[-4,3];
    xlim(xlim1);
    ylim(ylim1);

figure
% hold on
for i=1:5
%     subplot(3,3,i)
    ind=find(indici==i);
%     for j=1:length(ind)
    for j=1:20
        o=tra{ind(j)}';
        for k=1:9
            vec=o(:,k);
            vec=vec-meanp(k);
            vec=vec/stdp(k);
            o(:,k)=vec;
        end
        du=o*c;
        hold on
%         plot(du(:,1),du(:,2),colo(i))        
%         plot(du(end,1),du(end,2),['.' 'g'])
%         plot(du(1,pcadim1),du(1,pcadim2),['.' 'g'])
        plot(du(1,pcadim1),du(1,pcadim2),['.' colo(i)])
    end
    
    xlim(xlim1);
    ylim(ylim1);
%     return
       
    
end




le=cellfun(@length,tra)';

figure
hold on
for i=[2 5]
%     if(i~=2)
%         continue
%     end
%     subplot(3,3,i)
    ind=find(indici==i & le<400);
    length(ind)
%     for j=1:length(ind)
    for j=1:min(10,length(ind))
        o=tra{ind(j)}';
        for k=1:9
            vec=o(:,k);
            vec=vec-meanp(k);
            vec=vec/stdp(k);
            o(:,k)=vec;
        end
        du=o*c;
        hold on
        plot(du(:,pcadim1),du(:,pcadim2),colo(i))        
        plot(du(end,pcadim1),du(end,pcadim2),['.' 'k'])
        plot(du(1,pcadim1),du(1,pcadim2),['.' 'g'])
%         plot(du(1,1),du(1,2),['.' colo(i)],'MarkerSize',20)
    end
    
%     xlim(xlim1);
%     ylim(ylim1);
%     return
       
    
end


return

figure
hold on
for i=1:ngroups
    if(i~=2)
        continue
    end
%     subplot(3,3,i)
    ind=find(indici==i);
%     for j=1:length(ind)
    for j=1:10
        o=tra{j}';
        for k=1:9
            vec=o(:,k);
            vec=vec-meanp(k);
            vec=vec/stdp(k);
            o(:,k)=vec;
        end
        du=o*c;
        hold on
        plot(du(:,1),du(:,2),colo(i))        
        plot(du(end,1),du(end,2),['.' 'k'])
        plot(du(1,1),du(1,2),['.' 'g'])
%         plot(du(1,1),du(1,2),['.' colo(i)],'MarkerSize',20)
    end
    
    xlim(xlim1);
    ylim(ylim1);
%     return
       
    
end


return


return

o=tra{1}';
du=o*c;
plot(du(:,1),du(:,2),'k')



return

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
title('Cost for the 5 K-means groups');



%%% load names of arguments
load farm_F_0001 args


%%% histogram of the parameters, superimposed for the different groups
figure
for i=1:9
    subplot(3,3,i);
    hold on
    [~,x]=hist(b_orig(:,i),20);
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

    

clear
clc
acc=dir;
z=1;
for i=1:length(acc)
    n=acc(i).name;
    if(length(n)>5)
        load(n)
        pa(z,:)=pars;
        z=z+1;
    end
end

b=pa;
    
%normalize the parameters by imposing 0 mean and 1 st.dev.
for i=1:10
    p=b(:,i);
    p=p-mean(p);
    p=p/std(p);
    b(:,i)=p;
end

% 
% %normalize cost between 0 and 1 (and exclude top 5% outliers)
% costo=costo-min(costo);
% costo=costo/prctile(costo,95);
% costo(costo>1)=1;


%run pca
[c,s,l] = pca(b);


%plot % of variance explained by each component
figure;
plot(1:10,[100*l/sum(l)],'.-')
xlabel('component')
ylabel('% of variance')



pcadim1=1;
pcadim2=4;



%%% scatter plot of first 2 pca components, colored by cost
%%% where red=low cost ; green=high cost
figure
hold on
for i=1:size(b,1)
    plot(s(i,pcadim1),s(i,pcadim2),'.b')
end
xlabel('First PCA component')
ylabel('Second PCA component')
title('PCA Scatter Plot with Cost');

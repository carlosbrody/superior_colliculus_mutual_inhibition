clear
clc



first_run=0;

if(first_run)
%     %%% add .mat extension to all farm files (only run once!)
%     a=dir('FarmFields/farm_F_*');
%     for i=1:length(a)
%         n=a(i).name;
%         str=['mv ' n ' ' n '.mat'];
%         system(str);
%     end
% 

    
    %%% extract and save parameters and cost (only run once!)
    a=dir('FarmFields/farm_F*');
    
    b=nan(length(a),9);
    costo=nan(length(a),1);
    for i=1:length(a)
        load(['FarmFields/' a(i).name],'-mat','pars','scost','cb');
        b(i,:)=pars;
        costo(i)=scost;
        cbs(i)=cb;
    end
    save saved_data b costo cbs %parameters are saved in b

end

return


% 0.01

load saved_data

% figure
% hist(costo,50)
% return
z=1;
aaa=dir('giggio*.mat');
for i=1:length(aaa)
    if(costo(i)>-.005)
        continue
    end
    n=aaa(i).name;
    load(n);
    a=(sign(dif1)+1)/2;
    b=(sign(dif2)+1)/2;
    vals(z,1)=mean(a(:));
    vals(z,2)=mean(b(:));  
    vals(z,3)=cbs(i);
    z=z+1;
end

vals;

figure
hold on
plot(0.9,0.7,'ko','MarkerSize',30)
plot(0.9,0.7,'k.','MarkerSize',30)

ind1=find(vals(:,3)==0.02);
plot(vals(ind1,1),vals(ind1,2),'.b')

norm(vals(ind1,1)-0.9)+norm(vals(ind1,2)-0.7)


ind2=find(vals(:,3)==0.04);
plot(vals(ind2,1),vals(ind2,2),'.r')
xlim([0.5 1])
ylim([0.5 1])

norm(vals(ind2,1)-0.9)+norm(vals(ind2,2)-0.7)


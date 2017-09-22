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
    a=dir('/Users/mpagan/animals/farm_AB*');
    
    b=nan(length(a),9);
    costo=nan(length(a),1);
    for i=1:length(a)
        load(['/Users/mpagan/animals/' a(i).name],'-mat','pars','scost');
        b(i,:)=pars;
        costo(i)=scost;
    end
    save saved_data_AB b costo %parameters are saved in b

end

% 0.01

load saved_data_AB

% figure
% hist(costo,50)
% return
z=1;
aaa=dir('tuorlo*.mat');
for i=1:length(aaa)
    if(costo(i)>.005)
        continue
    end
    n=aaa(i).name;
    load(n);
    a=(sign(dif1)+1)/2;
    b=(sign(dif2)+1)/2;
    vals(z,1)=mean(a(:));
    vals(z,2)=mean(b(:));    
    z=z+1;
end

vals


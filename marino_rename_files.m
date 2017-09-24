clear
clc
close all



load saved_data

% figure

ind=find(costo<0.005);

for i=1:length(ind)    
    in=['farm_AB_' sprintf('%04d',ind(i))];        
    load(in,'-mat','pars');
    out=['farm_AC_' sprintf('%04d',i)];
    save(out,'pars');
    
end

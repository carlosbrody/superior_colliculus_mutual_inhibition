% clear
% clc

% load ../superior_colliculus_mutual_inhibition/grads1000

load gradsnew
grads

figure;
plot(1:size(grads,2),grads)


figure
hold on
for i=1:size(grads,1)
plot([0,grads(i,1)],[0,grads(i,2)],'.')
end


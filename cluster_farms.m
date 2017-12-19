clear

load compute_clustering/data_for_matlab_temp
nclusters = double(nclusters)';  % Came in as int, turn it into regular Matlab doubles vector

%%% select responses with cost lower than -0.0002
cost=[results(1).(cost_choice){:}];
ind=find(cost<threshold);

response=response(ind,:);
files = results(1).files(ind);



%%%% run PCA in response space
[c,s,l] = pca(response);



%%% find # of PCA dimensions that can recover at least 60% of the variance
ind=find(cumsum(l)/sum(l)>.6);
ndim=ind(1);



%%% reduce dimensionality of response space to ndim
score=s(:,1:ndim);

%%% use BIC to determine optimal number of clusters
if ~exist('nclusters', 'var')
    nclusters=2:7;
end

for i=1:length(nclusters)
    
    options = statset('MaxIter',10000,'Display','off');
    
    z=1;
    clear bicval
    %%% repeat the process 20 times; try-catch loop because it can fail
    for j=1:20
        try
            GM = fitgmdist(score,nclusters(i),'Options',options);
            bicval(z)=GM.BIC;
            z=z+1;
        catch
        end
    end
    %%% select median BIC value across the 20 repeats
    bic(i)=median(bicval);
end


%%% choose the number of groups with minimum BIC value
[~,i]=min(bic);
ngroups=nclusters(i);

%%% run k-means with chosen number of groups
idx=kmeans(response,ngroups);


%%% save:
%%% ndim = dimensionality of response space
%%% bic = BIC value for # of clusters from 2 to 5
%%% ngroups = number of groups with minimum BIC value
%%% idx = indices resulting from k-means clustering


save compute_clustering/data_for_julia_temp ndim bic ngroups idx files



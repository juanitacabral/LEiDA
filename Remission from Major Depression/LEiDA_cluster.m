function LEiDA_cluster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function computes k-means clustering of leading eigenvectors
%
% - Reads the pre-processed LEiDA_data
% - Saves the Cluster properties for all solutions (up to k=20)
% - Estimates the best number of clusters using Dunn's score
% 
% Saves the outputs to LEiDA_k_results.mat
%
% Joana Cabral May 2016
% Modified July 2016 with Kristina Rapuano
% joana.cabral@psych.ox.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create variable X with all vectors
% Rows are time points (observations)
% Collumns are brain areas (variables)

load LEiDA_data Leading_Eig
X=[];
for s=1:size(Leading_Eig,1)
	%X=cat(1,X,Leading_Eig{s});
    for task=1:2,
        X=cat(1,X,Leading_Eig{s,task});
    end
end
clear Leading_Eig
        
% Set maximum/minimum number of clusters
maxk=20;
mink=2;

%Kmeans clustering
opt= statset('UseParallel',1); %,'UseSubstreams',1);
Kmeans_results={};

for k=mink:maxk  
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','Replicates',20,'Display','final','Options',opt);   
    Kmeans_results{k}.IDX=IDX; % Cluster indices - numeric collumn vectos
    Kmeans_results{k}.C=C; % Cluster centroid locations
    Kmeans_results{k}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
    Kmeans_results{k}.D=D; % Distance from each point to every centroid

end

%Evaluate Clustering performance
distM_fcd=squareform(pdist(X,'cityblock'));
dunn_score=zeros(maxk,1);
for j=2:maxk
    dunn_score(j)=dunns(j,distM_fcd,Kmeans_results{j}.IDX);
    disp(['Performance for ' num2str(j) 'clusters'])
end;

[~,ind_max]=max(dunn_score);

disp(['Best clustering solution: ' num2str(ind_max) ' clusters']);

save('LEiDA_k_results_to_check_dunn','Kmeans_results_dunn','dunn_score_check')


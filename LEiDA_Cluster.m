function LEiDA_Cluster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function computes a k-means clustering algorith to the eigenvectors
%
% - Loads all the eigenvectors saved in LEiDA_data
% - Computes k-means clustering
% - Calculates the Dunn's score
% - Selects the optimal solution
%
% - Saves the optimal solution into LEiDA_Clusters.mat
%   Clusters =
%      IDX: Cluster index for each observation
%        C: Cluster centroids
%        D: Distance from each observation to every centroid
%     SUMD: Within-cluster sums of point-to-centroid distances
%
% Written by 
% Joana Cabral joana.cabral@psych.ox.ac.uk
% Paulo Marques paulo.c.g.marques@gmail.com
% Last edited May 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load LEiDA_data Leading_Eig

N_sub=size(Leading_Eig,1);

% Generate vector X concatenating all eigenvectors from all subjects and
% time points, where rows are observations and collumns are variables.

X=[];
for s=1:N_sub
    X=cat(1,X,squeeze(Leading_Eig(s,:,:))); 
end
clear Leading_Eig

%Kmeans clustering
maxk=20;
% opt= statset('UseParallel',1); %,'UseSubstreams',1);
% The options may vary according to the Matlab version
Kmeans_results=cell(1,20);

parfor k=2:maxk  
    disp(['Calculating for ' num2str(k) 'clusters'])
    [IDX, C, SUMD, D]=kmeans(X,k,'Distance','cityblock','Replicates',20,'Display','final'); %,'Options',opt);   
    Kmeans_results{k}.IDX=IDX;
    Kmeans_results{k}.C=C; 
    Kmeans_results{k}.SUMD=SUMD; 
    Kmeans_results{k}.D=D; 
end

%Evaluate Clustering performance
distM_fcd=squareform(pdist(X,'cityblock'));
dunn_score=zeros(maxk,1);
for j=2:maxk
    dunn_score(j)=dunns(j,distM_fcd,Kmeans_results{j}.IDX);
    disp(['Performance for ' num2str(j) 'clusters'])
end
[~,ind_max]=max(dunn_score);
disp(['Best clustering solution: ' num2str(ind_max) ' clusters']);

Clusters= Kmeans_results{ind_max};

save('LEiDA_Clusters','Clusters')
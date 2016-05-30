function LEiDA_analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function plots the cluster solutions on the cortex, and is matrix 
% format, in order to visualize the FC patterns.
% It alsso plots the probability of each cluster and the weighted sum of
% FC patterns in order to compare with the mean BOLD FC matrix
%
% Written by 
% Joana Cabral joana.cabral@psych.ox.ac.uk
% Last edited May 2016 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load LEiDA_Clusters Clusters
[N_Cl, N_ba]=size(Clusters.C);
h=hist(Clusters.IDX,N_Cl);
[y, ind]=sort(h,'descend');
V=Clusters.C(ind,:);
VVT_mean=zeros(N_ba);

% To reorder matrix plots
Order=[1:2:N_ba N_ba:-2:2];

figure
colormap(jet)    
% Pannel A
% Plot the cluster centroids over the cortex 
% Pannel C 
% Plot the centroids outerproduct

for c=1:N_Cl
    subplot(2,N_Cl+1,c)
    plot_nodes_in_cortex(V(c,:))
    title(['#' num2str(c)])
    axis off
    subplot(2,N_Cl+1,c+N_Cl+1)
    VVT=V(c,:)'*V(c,:);   
    imagesc(VVT(Order,Order))   
    axis square
    title('Outer product') 
    ylabel('Brain area #')
    xlabel('Brain area #') 
    VVT_mean=VVT_mean+squeeze(VVT)*y(c);
end

VVT_mean=VVT_mean/sum(y);

% Pannel B
% Plot the probabilities of each cluster
% The Colormap is adjusted for 5 clusters 
mymap=[1 0 0; 1 1 0; 0 1 0; 0 1 1; 0 0 1];

x=y/sum(y);
subplot(2,N_Cl+1,N_Cl+1)
hold on
for c=1:N_Cl
    bar(c,x(c),'Facecolor',mymap(c,:))
end
set(gca,'xtick',1:N_Cl)
box off
xlabel('State #')
ylabel('Probability')
xlim([0 6])

% Pannel D plot the weighted sum of outer products

subplot(3,N_Cl+1,(N_Cl+1)*2)
imagesc(VVT_mean(Order,Order))
title({'Weighted sum of VV^T'})
ylabel('Brain area')
xlabel('Brain area')
colorbar
axis square

load Static_FC
I_sup_diag=find(triu(ones(N_ba),1));
[cc p]=corrcoef(Static_FC_mean(I_sup_diag),VVT_mean(I_sup_diag))

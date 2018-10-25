function Bars_nodes_States

% Creates figure with bar plots of the 5 eigenvectors

load Kmeans_clusters.mat dunn_score Kmeans_results
[max_dunn,ind_max]=max(dunn_score);
Clusters= Kmeans_results{ind_max};
clear Kmeans_results dunn_score

[y, ind]=sort(Clusters.SUMD,'descend');
N=90;
Order=[1:2:N N:-2:2];
    
load AAL_labels.mat
figure

for c=1:ind_max
    cl=ind(c);
    subplot(1,5,c)
    barh(Clusters.C(cl,Order),'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
    set(gca,'YTick',1:90,'Fontsize',8)
    if c==1
        set(gca,'YTickLabel',label90(Order,:))
    else
        set(gca,'YTickLabel','')
    end
    title('Mean BOLD','Fontsize',10)
    title(['State #' num2str(c)])
    ylim([0 91])
    xlim([-.2 .2])
    grid on
end
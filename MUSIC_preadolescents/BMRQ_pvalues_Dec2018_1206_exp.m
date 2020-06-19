%function BMRQ_pvalues

load LEiDA_results_v6.mat P
load BMRQ_scores BMRQforanalysisOctexp

rangeK=3:15;

figure
for sc=1;
    
    scores=BMRQforanalysisOctexp(:,sc);
    
    for k=1:length(rangeK)
        
        for c=1:rangeK(k)
            
            Table=squeeze(P(:,:,k,c))';
            
            % P is a matrix with size 3x17x13x15
            
%             LEARNminREST=Table(:,2)-Table(:,1);
%             nLEARNminREST=Table(:,3)-Table(:,1);
            
            MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);         
            
%             [cc, p]=corrcoef(scores,LEARNminREST);
%             pval_LEARNminREST(sc,k,c)=p(2);
%             [cc, p]=corrcoef(scores,nLEARNminREST);
%             pval_nLEARNminREST(sc,k,c)=p(2);
            
            [cc, p]=corrcoef(scores,MUSICminREST);
            pval_MUSICminREST(sc,k,c)=p(2);
            
            
        end
    end
    
%     subplot(3,2,sc)
    hold on
    for k=1:length(rangeK) 
        semilogy(rangeK(k),squeeze(pval_MUSICminREST(sc,k,pval_MUSICminREST(sc,k,:)>0)),'*k'); 
%         semilogy(rangeK(k),squeeze(pval_LEARNminREST(sc,k,pval_LEARNminREST(sc,k,:)>0)),'*k'); 
%         semilogy(rangeK(k),squeeze(pval_nLEARNminREST(sc,k,pval_nLEARNminREST(sc,k,:)>0)),'ob');
    end
    semilogy(rangeK,0.05*ones(1,length(rangeK)),'r--','LineWidth',1)
    semilogy(rangeK,0.05./rangeK,'g--','LineWidth',1)
    
end

pval_score1=squeeze(pval_MUSICminREST(1,:,:));

% ylabel('Corr BMRQ Diff Prob')
% xlabel('Number of clusters K')
% set(gca,'XTick',3:15)
% ylim([1e-6 1])
% xlim([2 16])
% box off

load AAL_labels.mat label90

load LEiDA_results_v6.mat Kmeans_results rangeK
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

% COMPARE BMRQ 6 with k=7 c=2
figure
k=7;
c=2;
sc=1;

subplot(2,2,1)
scores=BMRQforanalysisOctexp(:,sc);
Table=squeeze(P(:,:,rangeK==k,c))';
MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);
plot(MUSICminREST,scores,'*')
xlabel(['Diff Prob of PC k=' num2str(k) ' c=' num2str(c)])
ylabel(['BMRQ score ' num2str(sc)])
[cc, p]=corrcoef(scores,MUSICminREST);
title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])
Vc=Kmeans_results{rangeK==k}.C;

subplot(2,2,[2 4])
V=Vc(c,Order);
hold on
barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
ylim([0 91])
xlim([-.15 .15])
set(gca,'YTick',1:N_areas,'Fontsize',8)
set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)

subplot(2,2,3)
plot_nodes_in_cortex(Vc(c,:))


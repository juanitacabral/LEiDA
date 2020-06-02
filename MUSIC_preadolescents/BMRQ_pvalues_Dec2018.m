%function BMRQ_pvalues

load LEiDA_results_v6.mat P
load BMRQ_scores BMRQforanalysisOct

rangeK=3:15;

figure
for sc=1:size(BMRQforanalysisOct,2)
    
    scores=BMRQforanalysisOct(:,sc);
    
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
    
    subplot(3,2,sc)
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
pval_score2=squeeze(pval_MUSICminREST(2,:,:));
pval_score3=squeeze(pval_MUSICminREST(3,:,:));
pval_score4=squeeze(pval_MUSICminREST(4,:,:));
pval_score5=squeeze(pval_MUSICminREST(5,:,:));
pval_score6=squeeze(pval_MUSICminREST(6,:,:));

ylabel('Corr BMRQ Diff Prob')
xlabel('Number of clusters K')
set(gca,'XTick',3:15)
ylim([1e-6 1])
xlim([2 16])
box off

load AAL_labels.mat label90

load LEiDA_results_v6.mat Kmeans_results rangeK
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

% COMPARE BMRQ 6 with k=7 c=2
figure
k=7;
c=2;
sc=6;

subplot(2,2,1)
scores=BMRQforanalysisOct(:,sc);
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


% COMPARE BMRQ 1 with k=14 c=10 
% score=1;
% k=14;
% c=10;
% 
% figure   
% subplot(2,2,1)
% scores=BMRQforanalysisOct(:,score);
% 
% Table=squeeze(P(:,:,rangeK==k,c))';
% MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);
% 
% plot(MUSICminREST,scores,'*')
% ylabel('Diff Prob of FC state 3 k 14')
% xlabel('BMRQ score 1')
% [cc, p]=corrcoef(scores,MUSICminREST);
% title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])
% 
% load LEiDA_results_v6.mat Kmeans_results rangeK
% load AAL_labels.mat label90
% Vc=Kmeans_results{rangeK==k}.C;
% N_areas=90;
% Order=[1:2:N_areas N_areas:-2:2];
% label90=label90(Order,:);
% 
% subplot(2,2,[2 4])
% V=Vc(c,Order);
% hold on
% barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
% barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
% ylim([0 91])
% xlim([-.15 .15])
% set(gca,'YTick',1:N_areas,'Fontsize',8)
% set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
% 
% % COMPARE BMRQ 1 with c=11 k=14 
% figure  
% 
% subplot(2,2,1)
% scores=BMRQforanalysisOct(:,1);
% k=14;
% c=11;
% Table=squeeze(P(:,:,rangeK==k,c))';
% MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);
% 
% plot(MUSICminREST,scores,'*')
% ylabel('Diff Prob of FC state 4 k 13')
% xlabel('BMRQ score 1')
% [cc, p]=corrcoef(scores,MUSICminREST);
% title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])
% 
% load LEiDA_results_v6.mat Kmeans_results rangeK
% load AAL_labels.mat label90
% Vc=Kmeans_results{rangeK==k}.C;
% N_areas=90;
% Order=[1:2:N_areas N_areas:-2:2];
% label90=label90(Order,:);
% 
% subplot(2,2,[2 4])
% V=Vc(c,Order);
% hold on
% barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
% barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
% ylim([0 91])
% xlim([-.15 .15])
% set(gca,'YTick',1:N_areas,'Fontsize',8)
% set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
% 
% % COMPARE BMRQ 1 with c=12 k=13 
% figure  
% 
% subplot(2,2,1)
% scores=BMRQforanalysisOct(:,1);
% k=14;
% c=12;
% Table=squeeze(P(:,:,rangeK==k,c))';
% MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);
% 
% plot(MUSICminREST,scores,'*')
% ylabel('Diff Prob of FC state 5 k 13')
% xlabel('BMRQ score 1')
% [cc, p]=corrcoef(scores,MUSICminREST);
% title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])
% 
% load LEiDA_results_v6.mat Kmeans_results rangeK
% load AAL_labels.mat label90
% Vc=Kmeans_results{rangeK==k}.C;
% N_areas=90;
% Order=[1:2:N_areas N_areas:-2:2];
% label90=label90(Order,:);
% 
% subplot(2,2,[2 4])
% V=Vc(c,Order);
% hold on
% barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
% barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
% ylim([0 91])
% xlim([-.15 .15])
% set(gca,'YTick',1:N_areas,'Fontsize',8)
% set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
% 
% % COMPARE BMRQ 1 with c=12 k=14 
% figure  
% 
% subplot(2,2,1)
% scores=BMRQforanalysisOct(:,1);
% k=14;
% c=12;
% Table=squeeze(P(:,:,rangeK==k,c))';
% MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);
% 
% plot(MUSICminREST,scores,'*')
% ylabel('Diff Prob of FC state 6 k 14')
% xlabel('BMRQ score 1')
% [cc, p]=corrcoef(scores,MUSICminREST);
% title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])
% 
% load LEiDA_results_v6.mat Kmeans_results rangeK
% load AAL_labels.mat label90
% Vc=Kmeans_results{rangeK==k}.C;
% N_areas=90;
% Order=[1:2:N_areas N_areas:-2:2];
% label90=label90(Order,:);
% 
% subplot(2,2,[2 4])
% V=Vc(c,Order);
% hold on
% barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
% barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
% ylim([0 91])
% xlim([-.15 .15])
% set(gca,'YTick',1:N_areas,'Fontsize',8)
% set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
% 
% % COMPARE BMRQ 1 with c=13 k=14 
% figure  
% 
% subplot(2,2,1)
% scores=BMRQforanalysisOct(:,1);
% k=14;
% c=13;
% Table=squeeze(P(:,:,rangeK==k,c))';
% MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);
% 
% plot(MUSICminREST,scores,'*')
% ylabel('Diff Prob of FC state 7 k 13')
% xlabel('BMRQ score 1')
% [cc, p]=corrcoef(scores,MUSICminREST);
% title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])
% 
% load LEiDA_results_v6.mat Kmeans_results rangeK
% load AAL_labels.mat label90
% Vc=Kmeans_results{rangeK==k}.C;
% N_areas=90;
% Order=[1:2:N_areas N_areas:-2:2];
% label90=label90(Order,:);
% 
% subplot(2,2,[2 4])
% V=Vc(c,Order);
% hold on
% barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
% barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
% ylim([0 91])
% xlim([-.15 .15])
% set(gca,'YTick',1:N_areas,'Fontsize',8)
% set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
% 

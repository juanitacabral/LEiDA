%%% SES socio-economic background (edu + income)

SES=[3	2.5	4.5	4.5	4	4.5	2	3.5	4	4.5	5 1.5 4 4.5	4.5	4.5	4.5];

load LEiDA_results_v6.mat P rangeK

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
        
        [cc, p]=corrcoef(SES,MUSICminREST);
        pval_MUSICminRESTvsYME(k,c)=p(2);
    end
end
    
figure
hold on
for k=1:length(rangeK) 
    semilogy(rangeK(k),squeeze(pval_MUSICminRESTvsYME(k,pval_MUSICminRESTvsYME(k,:)>0)),'*k'); 
%         semilogy(rangeK(k),squeeze(pval_LEARNminREST(sc,k,pval_LEARNminREST(sc,k,:)>0)),'*k'); 
%         semilogy(rangeK(k),squeeze(pval_nLEARNminREST(sc,k,pval_nLEARNminREST(sc,k,:)>0)),'ob');
end
semilogy(rangeK,0.05*ones(1,length(rangeK)),'r--')
semilogy(rangeK,0.05./rangeK,'g--')
% 



figure

subplot(2,2,1)
k=6;
c=1;
Table=squeeze(P(:,:,rangeK==k,c))';
MUSICminREST=(Table(:,2)+Table(:,3))./2-Table(:,1);

plot(MUSICminREST,SES,'*')
xlabel(['Diff Prob of CC k=' num2str(k) ' c=' num2str(c)])
ylabel('Months Music Education')
[cc, p]=corrcoef(SES,MUSICminREST);
title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])

load LEiDA_results_v6.mat Kmeans_results rangeK
load AAL_labels.mat label90
Vc=Kmeans_results{rangeK==k}.C;
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

subplot(2,2,[2 4])
V=Vc(c,Order);
hold on
barh(find(V<0),V(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
barh(find(V>=0),V(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
ylim([0 91])
xlim([-.15 .15])
set(gca,'YTick',1:N_areas,'Fontsize',8)
set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
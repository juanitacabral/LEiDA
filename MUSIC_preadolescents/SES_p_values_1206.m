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



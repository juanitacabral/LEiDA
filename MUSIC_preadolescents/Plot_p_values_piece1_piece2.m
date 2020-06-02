function Plot_p_values_MC_2

load LEiDA_results_v6.mat P_pval_LEARNvsNONL rangeK P_pval_RESTvsMUSIC Kmeans_results

rangeK=3:15;


P_pval=P_pval_RESTvsMUSIC;

sigC=zeros(1,length(rangeK));

for k=1:length(rangeK)   
    [Min_p_value(k), sigC(k)]=min(P_pval(k,P_pval(k,:)>0));    
end



figure
semilogy(rangeK,0.05*ones(1,length(rangeK)),'r--','LineWidth',1)
hold on
semilogy(rangeK,0.05./rangeK,'g--','LineWidth',1)
semilogy(rangeK,0.05./sum(rangeK)*ones(1,length(rangeK)),'b--','LineWidth',1)

for k=1:length(rangeK) 
    for c=1:rangeK(k)
        if P_pval(k,c)<0.05 && P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),P_pval(k,c),'*r');
        end
        if P_pval(k,c)>0.05
            semilogy(rangeK(k),P_pval(k,c),'*k');
        end
        if P_pval(k,c)<(0.05/rangeK(k)) && P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*g');
        end
        if P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*b');
        end
    end
   
    
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)<0.05/sum(rangeK)),'*b');  
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)>0.05/sum(rangeK)),'*g');  
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)>0.05/rangeK(k)),'*r');    
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)>0.05),'*k'); 
end
%semilogy(rangeK,Min_p_value,'*b')
%semilogy(rangeK,Min_p_value.*rangeK,'g*-')

ylabel('Prob Rest vs Music (p-value)')
xlabel('Number of clusters K')
set(gca,'XTick',3:15)
ylim([1e-6 1])
xlim([2 16])
box off

figure % FIGURE SHOWING MOST DIFFERENT NETWORK IN EACH K
load AAL_labels.mat label90
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

V_all=zeros(N_areas,length(rangeK));

for k=1:length(rangeK)
    c=sigC(k);
    V=Kmeans_results{k}.C(c,:);
    V=V(Order);
    V_all(:,k)=V;
    
    subplot(1,length(rangeK),k)
    hold on
    barh(V.*(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(V.*(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 91])
    xlim([-.15 .15])
    
    set(gca,'Ydir','reverse')
    
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',8)    
    if k==1
        set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
    else
        set(gca,'YTickLabel',[])
    end

    if Min_p_value(k)<0.05 && Min_p_value(k)>(0.05/k)
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','r')
    elseif Min_p_value(k)>0.05
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','k')
    elseif Min_p_value(k)<(0.05/k) && Min_p_value(k)>(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','g')
    elseif Min_p_value(k)<(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','b')        
    end

end


% LEARNED VS NON-LEARNED
P_pval=P_pval_LEARNvsNONL;

sigC=zeros(1,length(rangeK));

for k=1:length(rangeK)   
    [Min_p_value(k), sigC(k)]=min(P_pval(k,P_pval(k,:)>0));    
end

figure
semilogy(rangeK,0.05*ones(1,length(rangeK)),'r--','LineWidth',1)
hold on
semilogy(rangeK,0.05./rangeK,'g--','LineWidth',1)
semilogy(rangeK,0.05./sum(rangeK)*ones(1,length(rangeK)),'b--','LineWidth',1)


for k=1:length(rangeK) 
    for c=1:rangeK(k)
        if P_pval(k,c)<0.05 && P_pval(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),P_pval(k,c),'*r');
        end
        if P_pval(k,c)>0.05
            semilogy(rangeK(k),P_pval(k,c),'*k');
        end
        if P_pval(k,c)<(0.05/rangeK(k)) && P_pval(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*g');
        end
        if P_pval(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),P_pval(k,c),'*b');
        end
    end
   
    
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)<0.05/sum(rangeK)),'*b');  
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)>0.05/sum(rangeK)),'*g');  
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)>0.05/rangeK(k)),'*r');    
%     semilogy(rangeK(k),P_pval(k,P_pval(k,:)>0.05),'*k'); 
end

%semilogy(rangeK,Min_p_value,'*b')
%semilogy(rangeK,Min_p_value.*rangeK,'g*-')

ylabel('Prob Music 1 vs Music 2 (p-value)')
xlabel('Number of clusters K')
set(gca,'XTick',3:15)
ylim([1e-6 1])
xlim([2 16])
box off

figure % FIGURE SHOWING MOST DIFFERENT NETWORK IN EACH K
load AAL_labels.mat label90
N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

V_all=zeros(N_areas,length(rangeK));

for k=1:length(rangeK)
    c=sigC(k);
    V=Kmeans_results{k}.C(c,:);
    V=V(Order);
    V_all(:,k)=V;
    
    subplot(1,length(rangeK),k)
    hold on
    barh(V.*(V<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(V.*(V>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 91])
    xlim([-.15 .15])
    set(gca,'Ydir','reverse')
    
    grid on
    set(gca,'YTick',1:N_areas,'Fontsize',8)    
    if k==1
        set(gca,'YTickLabel',label90(1:end,:),'Fontsize',6)
    else
        set(gca,'YTickLabel',[])
    end
    
    if Min_p_value(k)<0.05 && Min_p_value(k)>(0.05/rangeK(k))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','r')
    elseif Min_p_value(k)>0.05
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','k')
    elseif Min_p_value(k)<(0.05/rangeK(k)) && Min_p_value(k)>(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','g')
    elseif Min_p_value(k)<(0.05/sum(rangeK))
        title({['K=' num2str(rangeK(k))],['p=' num2str(Min_p_value(k))]},'Fontsize',10,'color','b')        
    end
end


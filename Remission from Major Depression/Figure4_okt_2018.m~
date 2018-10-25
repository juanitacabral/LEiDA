
%for figure 
%Patients vs Controls, plot most significant networks

%For probability and lifetime k=10, c=8

load LEiDA_k_results Kmeans_results
load Between_results_sept_2018 P pvalue_PvsC LT LT_pval_PvsC

load AAL_labels.mat label90
Order=[1:2:89 90:-2:2];

%Figure 4. k=10 c=8 is significant after correction  for multiple
%comparisons for prob and lifetime. % to make the plot of brain areas in
%the network for figure 4B:



for t=1:2
    if t==1
        task='Neutral';
    else
        task='Sad';
    end
    for k=10
        for c=8
  
            if pvalue_PvsC(t,k-1,c)<(0.5) && pvalue_PvsC(t,k-1,c)>0
                V=Kmeans_results{k}.C(c,:);
                figure('Name',['Patients vs Controls under ' task ' mood K=' num2str(k) ' c=' num2str(c) ])
%                 subplot(4,2,1)
%                 plot_nodes_in_cortex(V);
                subplot(4,2,[2 4 6 8])
                %V=V(Order);
                i_red=find(V>0);
                i_blue=find(V<0);
                barh(i_red,V(i_red),'FaceColor',[0.98 0.85 0.37],'EdgeColor','none')
                hold on
                barh(i_blue,V(i_blue),'FaceColor',[0.3010, 0.7450, 0.9330],'EdgeColor','none')
                set(gca,'YTick',1:90,'Fontsize',6)
                set(gca,'YTickLabel',label90)
                %ylim([0 91])
                xlim([-.15 .15])
                grid off
                
            end
        end
    end
end

%For figure 4A

k=10;
[~, ind_sort]=sort(Kmeans_results{k}.SUMD,'descend');

%Plot the significant network for k=10, c=8 on the brain
V=Kmeans_results{k}.C(c,:);
plot_nodes_in_cortex_new(V);

%Make probability and lifetime for the significant state in one figure, first neutral, than sad.
%for figure 4D
%MK made figure 4C

%rrMDDneutral
    state1np=P(k-1,1:51,1,ind_sort(1))';
    state2np=P(k-1,1:51,1,ind_sort(2))';
    state3np=P(k-1,1:51,1,ind_sort(3))';
    state4np=P(k-1,1:51,1,ind_sort(4))';
    state5np=P(k-1,1:51,1,ind_sort(5))';
    state6np=P(k-1,1:51,1,ind_sort(6))';
    state7np=P(k-1,1:51,1,ind_sort(7))';
    state8np=P(k-1,1:51,1,ind_sort(8))';
    state9np=P(k-1,1:51,1,ind_sort(9))';
    state10np=P(k-1,1:51,1,ind_sort(10))';
    
    %controlsneutral
    state1cp=P(k-1,52:86,1,ind_sort(1))';
    state2cp=P(k-1,52:86,1,ind_sort(2))';
    state3cp=P(k-1,52:86,1,ind_sort(3))';
    state4cp=P(k-1,52:86,1,ind_sort(4))';
    state5cp=P(k-1,52:86,1,ind_sort(5))';
    state6cp=P(k-1,52:86,1,ind_sort(6))';
    state7cp=P(k-1,52:86,1,ind_sort(7))';
    state8cp=P(k-1,52:86,1,ind_sort(8))';
    state9cp=P(k-1,52:86,1,ind_sort(9))';
    state10cp=P(k-1,52:86,1,ind_sort(10))';
    
    %rrMDDsad
    state1nps=P(k-1,1:51,2,ind_sort(1))';
    state2nps=P(k-1,1:51,2,ind_sort(2))';
    state3nps=P(k-1,1:51,2,ind_sort(3))';
    state4nps=P(k-1,1:51,2,ind_sort(4))';
    state5nps=P(k-1,1:51,2,ind_sort(5))';
    state6nps=P(k-1,1:51,2,ind_sort(6))';
    state7nps=P(k-1,1:51,2,ind_sort(7))';
    state8nps=P(k-1,1:51,2,ind_sort(8))';
    state9nps=P(k-1,1:51,2,ind_sort(9))';
    state10nps=P(k-1,1:51,2,ind_sort(10))';
    
    %controlssad
    state1cps=P(k-1,52:86,2,ind_sort(1))';
    state2cps=P(k-1,52:86,2,ind_sort(2))';
    state3cps=P(k-1,52:86,2,ind_sort(3))';
    state4cps=P(k-1,52:86,2,ind_sort(4))';
    state5cps=P(k-1,52:86,2,ind_sort(5))';
    state6cps=P(k-1,52:86,2,ind_sort(6))';
    state7cps=P(k-1,52:86,2,ind_sort(7))';
    state8cps=P(k-1,52:86,2,ind_sort(8))';
    state9cps=P(k-1,52:86,2,ind_sort(9))';
    state10cps=P(k-1,52:86,2,ind_sort(10))';
    
    %for lifetimes
       %rrMDDneutral
    state1npLT=LT(k-1,1:51,1,ind_sort(1))';
    state2npLT=LT(k-1,1:51,1,ind_sort(2))';
    state3npLT=LT(k-1,1:51,1,ind_sort(3))';
    state4npLT=LT(k-1,1:51,1,ind_sort(4))';
    state5npLT=LT(k-1,1:51,1,ind_sort(5))';
    state6npLT=LT(k-1,1:51,1,ind_sort(6))';
    state7npLT=LT(k-1,1:51,1,ind_sort(7))';
    state8npLT=LT(k-1,1:51,1,ind_sort(8))';
    state9npLT=LT(k-1,1:51,1,ind_sort(9))';
    state10npLT=LT(k-1,1:51,1,ind_sort(10))';
    
    %controlsneutral
    state1ncLT=LT(k-1,52:86,1,ind_sort(1))';
    state2ncLT=LT(k-1,52:86,1,ind_sort(2))';
    state3ncLT=LT(k-1,52:86,1,ind_sort(3))';
    state4ncLT=LT(k-1,52:86,1,ind_sort(4))';
    state5ncLT=LT(k-1,52:86,1,ind_sort(5))';
    state6ncLT=LT(k-1,52:86,1,ind_sort(6))';
    state7ncLT=LT(k-1,52:86,1,ind_sort(7))';
    state8ncLT=LT(k-1,52:86,1,ind_sort(8))';
    state9ncLT=LT(k-1,52:86,1,ind_sort(9))';
    state10ncLT=LT(k-1,52:86,1,ind_sort(10))';
    
    %rrMDDsad
    state1spLT=LT(k-1,1:51,2,ind_sort(1))';
    state2spLT=LT(k-1,1:51,2,ind_sort(2))';
    state3spLT=LT(k-1,1:51,2,ind_sort(3))';
    state4spLT=LT(k-1,1:51,2,ind_sort(4))';
    state5spLT=LT(k-1,1:51,2,ind_sort(5))';
    state6spLT=LT(k-1,1:51,2,ind_sort(6))';
    state7spLT=LT(k-1,1:51,2,ind_sort(7))';
    state8spLT=LT(k-1,1:51,2,ind_sort(8))';
    state9spLT=LT(k-1,1:51,2,ind_sort(9))';
    state10spLT=LT(k-1,1:51,2,ind_sort(10))';
    
    %controlssad
    state1scLT=LT(k-1,52:86,2,ind_sort(1))';
    state2scLT=LT(k-1,52:86,2,ind_sort(2))';
    state3scLT=LT(k-1,52:86,2,ind_sort(3))';
    state4scLT=LT(k-1,52:86,2,ind_sort(4))';
    state5scLT=LT(k-1,52:86,2,ind_sort(5))';
    state6scLT=LT(k-1,52:86,2,ind_sort(6))';
    state7scLT=LT(k-1,52:86,2,ind_sort(7))';
    state8scLT=LT(k-1,52:86,2,ind_sort(8))';
    state9scLT=LT(k-1,52:86,2,ind_sort(9))';
    state10scLT=LT(k-1,52:86,2,ind_sort(10))';

subplot(4,2,3)
    ydata = [mean(state4np) mean(state4cp); mean(state4nps) mean(state4cps)];
    error=[std(state4np)/sqrt(numel(state4np)) std(state4cp)/sqrt(numel(state4cp)); std(state4nps)/sqrt(numel(state4nps)) std(state4cps)/sqrt(numel(state4cps))];
    h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
    set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})
    %set(gca,'XTickLabel',{'neutral', 'sad'})
    legend({'rrMDD','controls'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
    ylabel('Probability', 'Fontsize',10) % y-axis label
    set(h(1),'FaceColor',[0 0.8 0]);
    set(h(2),'FaceColor',[0 0.4 1]);
    title('Probability of FC-state in neutral and sad mood')
    text(2,max([mean(state4nps) mean(state4cps)])+0.005,'*','FontSize',25)
    text(1,max([mean(state4np) mean(state4cp)])+.01,'**','FontSize',25)
    box off
    
    subplot(4,2,5)
    ydata = [mean(state4npLT) mean(state4ncLT); mean(state4spLT) mean(state4scLT)];
    error=[std(state4npLT)/sqrt(numel(state4npLT)) std(state4ncLT)/sqrt(numel(state4ncLT)); std(state4spLT)/sqrt(numel(state4spLT)) std(state4scLT)/sqrt(numel(state4scLT))]; 
    h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
    %set(gca,'XTickLabel',{'neutral', 'sad'})
    set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})
    legend({'rrMDD','Control'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
    ylabel('Lifetime', 'Fontsize',10) % y-axis label
    set(h(1),'FaceColor',[1 0.8 0]);
    set(h(2),'FaceColor',[0.6 0 1]);
    title('Duration of FC-state in neutral and sad mood')
    text(1,max([mean(state4npLT) mean(state4ncLT)])+0.18,'**','FontSize',25)
    %text(1,max([mean(state4np) mean(state4cp)])+.01,'**','FontSize',25)
    box off 
    
   %to determine the AAL areas that are in the network, this is used for giving the networks a name:
   V=Kmeans_results{10}.C(ind_sort(3),:);
    find(V>0)
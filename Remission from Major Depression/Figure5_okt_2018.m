%figure 5

load LEiDA_k_results Kmeans_results
load Between_results_sept_2018 P LT pvalue_PvsC LT_pval_PvsC
load Within_results_okt_2018 pvalue_NvsS LT_pval_NvsS


k=10;
[~, ind_sort]=sort(Kmeans_results{k}.SUMD,'descend');


%for every state make plot for rrMDD and controls and for neutral and sad
%mood


%Make probability in one figure, first neutral, than sad.

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
   


%to check what areas the FC patterns consist of use:
 V=Kmeans_results{10}.C(ind_sort(3),:);
 load AAL_labels.mat label90
 label90(V>0,:)

 %Use this for figure!
    
%whole thing together Probability
% same for sad mood
    %for PvsC in neutral mood,only 4 is significant after MC
    %for PvsC in sad mood: only 4 is significant before MC
    %for sad vs neutral mood within rrMDD: state 3 and 7 differ after MC.
    %before MC: 1, 4 and 8 and 10 also differ
    %for sad vs neutral mood within controls: no differences before and after MC.
subplot(2,20,[1:20]) 
ydata = [mean(state1np) mean(state1cp); mean(state1nps) mean(state1cps); mean(state2np) mean(state2cp); mean(state2nps) mean(state2cps);  mean(state3np) mean(state3cp); mean(state3nps) mean(state3cps);  mean(state4np) mean(state4cp); mean(state4nps) mean(state4cps);  mean(state5np) mean(state5cp); mean(state5nps) mean(state5cps); mean(state6np) mean(state6cp); mean(state6nps) mean(state6cps); mean(state7np) mean(state7cp); mean(state7nps) mean(state7cps); mean(state8np) mean(state8cps); mean(state8nps) mean(state8cps); mean(state9np) mean(state9cp); mean(state9nps) mean(state9cps); mean(state10np) mean(state10cp); mean(state10nps) mean(state10cps);];
error=[std(state1np)/sqrt(numel(state1np)) std(state1cp)/sqrt(numel(state1cp)); std(state1nps)/sqrt(numel(state1nps)) std(state1cps)/sqrt(numel(state1cps)); std(state2np)/sqrt(numel(state2np)) std(state2cp)/sqrt(numel(state2cp));
std(state2nps)/sqrt(numel(state2nps)) std(state2cps)/sqrt(numel(state2cps));  std(state3np)/sqrt(numel(state3np)) std(state3cp)/sqrt(numel(state3cp));
std(state3nps)/sqrt(numel(state3nps)) std(state3cps)/sqrt(numel(state3cps));  std(state4np)/sqrt(numel(state4np)) std(state4cp)/sqrt(numel(state4cp)); 
std(state4nps)/sqrt(numel(state4nps)) std(state4cps)/sqrt(numel(state4cps)); std(state5np)/sqrt(numel(state5np)) std(state5cp)/sqrt(numel(state5cp)); 
std(state5nps)/sqrt(numel(state5nps)) std(state5cps)/sqrt(numel(state5cps)); std(state6np)/sqrt(numel(state6np)) std(state6cp)/sqrt(numel(state6cp)); std(state6np)/sqrt(numel(state6np)) std(state6cps)/sqrt(numel(state6cps)); std(state7np)/sqrt(numel(state7np)) std(state7cp)/sqrt(numel(state7cp)); std(state7nps)/sqrt(numel(state7nps)) std(state7cps)/sqrt(numel(state7cps)); std(state8np)/sqrt(numel(state8np)) std(state8cp)/sqrt(numel(state8cp));
std(state8nps)/sqrt(numel(state8nps)) std(state8cps)/sqrt(numel(state8cps));  std(state9np)/sqrt(numel(state9np)) std(state9cp)/sqrt(numel(state9cp));
std(state9nps)/sqrt(numel(state9nps)) std(state9cps)/sqrt(numel(state9cps));  std(state10np)/sqrt(numel(state10np)) std(state10cp)/sqrt(numel(state10cp));
std(state10nps)/sqrt(numel(state10nps)) std(state10cps)/sqrt(numel(state10cps))];

h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})    
%set(gca,'XTickLabel',{'Neutral', 'Sad', 'Neutral',  'Sad', 'Neutral',  'Sad', 'Neutral',  'Sad','Neutral', 'Sad', 'Neutral', 'Sad', 'Neutral', 'Sad', 'Neutral', 'Sad','Neutral', 'Sad', 'Neutral','Sad','Neutral', 'Sad', 'Neutral','Sad','Neutral', 'Sad','Neutral', 'Sad', 'Neutral', 'Sad', 'Neutral',  'Sad', 'Neutral', 'Sad', 'Neutral', 'Sad'})
%I added the label neutral and sad myself after making the figure    
legend({'rrMDD','Control'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
   % ylabel('Probability', 'Fontsize',10) % y-axis label
    set(h(1),'FaceColor',[0 0.8 0]);
    set(h(2),'FaceColor',[0 0.4 1]);
    title('Neutral mood probability')
    %text(3.8,0.09,'**','FontSize',28)
    box off
    
%Sad and neutral together LT
%for LT PvsC, in neutral mood state 4 is sign different after multiple
    %comparisons. before MC, no other states are different.
%for LT P vs C in sad mood no difference after MC. Before MC: 10 and 7
    %differ
    %for LT neutral vs sad, for rrMDD state 4 is different after MC.
    %before: 1,3,4,7,8,10. No differences for controls before and after MC.
    %I added this myself after making the figure
subplot(2,20,[21:40])
ydata=[mean(state1npLT) mean(state1ncLT); mean(state1spLT) mean(state1scLT); mean(state2npLT) mean(state2ncLT); mean(state2spLT) mean(state2scLT); mean(state3npLT) mean(state3ncLT); mean(state3spLT) mean(state3scLT); mean(state4npLT) mean(state4ncLT); mean(state4spLT) mean(state4scLT); mean(state5npLT) mean(state5ncLT); mean(state5spLT) mean(state5scLT); mean(state6npLT) mean(state6ncLT); mean(state6spLT) mean(state6scLT); mean(state7npLT) mean(state7ncLT); mean(state7spLT) mean(state7scLT); mean(state8npLT) mean(state8ncLT); mean(state8spLT) mean(state8scLT); mean(state9npLT) mean(state9ncLT); mean(state9spLT) mean(state9scLT); mean(state10npLT) mean(state10ncLT); mean(state10spLT) mean(state10scLT);];
error=[std(state1npLT)/sqrt(numel(state1npLT)) std(state1ncLT)/sqrt(numel(state1scLT)); std(state1spLT)/sqrt(numel(state1spLT)) std(state1scLT)/sqrt(numel(state1scLT));
std(state2npLT)/sqrt(numel(state2npLT)) std(state2ncLT)/sqrt(numel(state2ncLT));  std(state2spLT)/sqrt(numel(state2spLT)) std(state2scLT)/sqrt(numel(state2scLT));
std(state3npLT)/sqrt(numel(state3npLT)) std(state3ncLT)/sqrt(numel(state3ncLT));  std(state3spLT)/sqrt(numel(state3spLT)) std(state3scLT)/sqrt(numel(state3scLT)); 
std(state4npLT)/sqrt(numel(state4npLT)) std(state4ncLT)/sqrt(numel(state4ncLT)); std(state4spLT)/sqrt(numel(state4spLT)) std(state4scLT)/sqrt(numel(state4scLT));
std(state5npLT)/sqrt(numel(state5npLT)) std(state5ncLT)/sqrt(numel(state5ncLT)); std(state5spLT)/sqrt(numel(state5spLT)) std(state5scLT)/sqrt(numel(state5scLT));
std(state6npLT)/sqrt(numel(state6npLT)) std(state6ncLT)/sqrt(numel(state6ncLT));  std(state6spLT)/sqrt(numel(state6spLT)) std(state6scLT)/sqrt(numel(state6scLT));
std(state7npLT)/sqrt(numel(state7npLT)) std(state7ncLT)/sqrt(numel(state7ncLT)); std(state7spLT)/sqrt(numel(state7spLT)) std(state7scLT)/sqrt(numel(state7scLT)); 
std(state8npLT)/sqrt(numel(state8npLT)) std(state8ncLT)/sqrt(numel(state8ncLT)); std(state8spLT)/sqrt(numel(state8spLT)) std(state8scLT)/sqrt(numel(state8scLT)); 
std(state9npLT)/sqrt(numel(state9npLT)) std(state9ncLT)/sqrt(numel(state9ncLT)); std(state9spLT)/sqrt(numel(state9spLT)) std(state9scLT)/sqrt(numel(state9scLT)); std(state10npLT)/sqrt(numel(state10npLT)) std(state10ncLT)/sqrt(numel(state10ncLT));
std(state10spLT)/sqrt(numel(state10spLT)) std(state10scLT)/sqrt(numel(state10scLT))];
    h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
%ylabel('Lifetime', 'Fontsize',10) % y-axis label
    set(h(1),'FaceColor',[1 0.8 0]);
    set(h(2),'FaceColor',[0.6 0 1]);
    set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})
    %I added the labels neutral and sad myself after making the figure
    legend({'rrMDD','Control'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
    title('Sad mood duration')
    box off
    
    %to make seperate plots (did not use this for figure 5 
    
    subplot(2,10,[1:5])
    colormap(winter)
    ydata = [mean(state1np) mean(state1cp); mean(state2np) mean(state2cp); mean(state3np) mean(state3cp); mean(state4np) mean(state4cp); mean(state5np) mean(state5cp); mean(state6np) mean(state6cp); mean(state7np) mean(state7cp); mean(state8np) mean(state8cp); mean(state9np) mean(state9cp); mean(state10np) mean(state10cp);];
    error=[std(state1np)/sqrt(numel(state1np)) std(state1cp)/sqrt(numel(state1cp)); std(state2np)/sqrt(numel(state2np)) std(state2cp)/sqrt(numel(state2cp)); std(state3np)/sqrt(numel(state3np)) std(state3cp)/sqrt(numel(state3cp)); std(state4np)/sqrt(numel(state4np)) std(state4cp)/sqrt(numel(state4cp)); std(state5np)/sqrt(numel(state5np)) std(state5cp)/sqrt(numel(state5cp)); std(state6np)/sqrt(numel(state6np)) std(state6cp)/sqrt(numel(state6cp)); std(state7np)/sqrt(numel(state7np)) std(state7cp)/sqrt(numel(state7cp)); std(state8np)/sqrt(numel(state8np)) std(state8cp)/sqrt(numel(state8cp)); std(state9np)/sqrt(numel(state9np)) std(state9cp)/sqrt(numel(state9cp)); std(state10np)/sqrt(numel(state10np)) std(state10cp)/sqrt(numel(state10cp))];
    h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
    set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})
    legend({'rrMDD','Control'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
   % ylabel('Probability', 'Fontsize',10) % y-axis label
    set(h(1),'FaceColor',[0 0.8 0]);
    set(h(2),'FaceColor',[0 0.4 1]);
    title('Neutral mood probability')
    text(3.8,0.09,'**','FontSize',28)
    box off
    

    %for LT PvsC, in neutral mood state 4 is sign different after multiple
    %comparisons. before MC, no other states are different.
    subplot(2,10,[6:10])
    colormap(winter)
    ydata = [mean(state1npLT) mean(state1ncLT); mean(state2npLT) mean(state2ncLT); mean(state3npLT) mean(state3ncLT); mean(state4npLT) mean(state4ncLT); mean(state5npLT) mean(state5ncLT); mean(state6npLT) mean(state6ncLT); mean(state7npLT) mean(state7ncLT); mean(state8npLT) mean(state8ncLT); mean(state9npLT) mean(state9ncLT); mean(state10npLT) mean(state10ncLT);];
    error=[std(state1npLT)/sqrt(numel(state1npLT)) std(state1ncLT)/sqrt(numel(state1ncLT)); std(state2npLT)/sqrt(numel(state2npLT)) std(state2ncLT)/sqrt(numel(state2ncLT)); std(state3npLT)/sqrt(numel(state3npLT)) std(state3ncLT)/sqrt(numel(state3ncLT)); std(state4npLT)/sqrt(numel(state4npLT)) std(state4ncLT)/sqrt(numel(state4ncLT)); std(state5npLT)/sqrt(numel(state5npLT)) std(state5ncLT)/sqrt(numel(state5ncLT)); std(state6npLT)/sqrt(numel(state6npLT)) std(state6ncLT)/sqrt(numel(state6ncLT)); std(state7npLT)/sqrt(numel(state7npLT)) std(state7ncLT)/sqrt(numel(state7ncLT)); std(state8npLT)/sqrt(numel(state8npLT)) std(state8ncLT)/sqrt(numel(state8ncLT)); std(state9npLT)/sqrt(numel(state9npLT)) std(state9ncLT)/sqrt(numel(state9ncLT)); std(state10npLT)/sqrt(numel(state10npLT)) std(state10ncLT)/sqrt(numel(state10ncLT))];
    h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
    set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})
    %legend({'rrMDD','Control'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
    %ylabel('Lifetime', 'Fontsize',10) % y-axis label
     set(h(1),'FaceColor',[0 0.8 0]);
     set(h(2),'FaceColor',[0 0.4 1]);
     text(3.8,2.9,'**','FontSize',28)
%      text(6.8,2.5,'+','FontSize',22)
%      text(10,2.5,'+','FontSize',22)
     box off
     
 
    
    % same for sad mood
    %for PvsC in sad mood: only 4 is significant before MC
    %for sad vs neutral mood within rrMDD: state 3 and 7 differ after MC.
    %before MC: 1, 4 and 8 and 10 also differ
    %for sad vs neutral mood within controls: no differences before and after MC.  
    
    
    subplot(2,10,[11:15])
    colormap(winter)
    ydata = [mean(state1nps) mean(state1cps); mean(state2nps) mean(state2cps); mean(state3nps) mean(state3cps); mean(state4nps) mean(state4cps); mean(state5nps) mean(state5cps); mean(state6nps) mean(state6cps); mean(state7nps) mean(state7cps); mean(state8nps) mean(state8cps); mean(state9nps) mean(state9cps); mean(state10nps) mean(state10cps);];
    error=[std(state1nps)/sqrt(numel(state1nps)) std(state1cps)/sqrt(numel(state1cps)); std(state2nps)/sqrt(numel(state2nps)) std(state2cps)/sqrt(numel(state2cps)); std(state3nps)/sqrt(numel(state3nps)) std(state3cps)/sqrt(numel(state3cps)); std(state4nps)/sqrt(numel(state4nps)) std(state4cps)/sqrt(numel(state4cps)); std(state5nps)/sqrt(numel(state5nps)) std(state5cps)/sqrt(numel(state5cps)); std(state6nps)/sqrt(numel(state6nps)) std(state6cps)/sqrt(numel(state6cps)); std(state7nps)/sqrt(numel(state7nps)) std(state7cps)/sqrt(numel(state7cps)); std(state8nps)/sqrt(numel(state8nps)) std(state8cps)/sqrt(numel(state8cps)); std(state9nps)/sqrt(numel(state9nps)) std(state9cps)/sqrt(numel(state9cps)); std(state10nps)/sqrt(numel(state10nps)) std(state10cps)/sqrt(numel(state10cps))];
    h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
    set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})
    legend({'rrMDD','Control'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
   % ylabel('Probability', 'Fontsize',10) % y-axis label
    title('Sad mood probability')
    set(h(1),'FaceColor',[1 0.8 0]);
    set(h(2),'FaceColor',[0.6 0 1]);
    text(4,0.09,'*','FontSize',28)
    text(2.8,0.14,'++','FontSize',20)
    text(6.8,0.085,'++','FontSize',20)
    text(3.8,0.14,'+','FontSize',20)
    text(7.8,0.085,'+','FontSize',20)
    box off
    
    %for LT P vs C in sad mood no difference after MC. Before MC: 10 and 7
    %differ
    %for LT neutral vs sad, for rrMDD state 4 is different after MC.
    %before: 1,3,4,7,8,10. No differences for controls before and after MC. 
    subplot(2,10,[16:20])
    colormap(winter)
    ydata = [mean(state1spLT) mean(state1scLT); mean(state2spLT) mean(state2scLT); mean(state3spLT) mean(state3scLT); mean(state4spLT) mean(state4scLT); mean(state5spLT) mean(state5scLT); mean(state6spLT) mean(state6scLT); mean(state7spLT) mean(state7scLT); mean(state8spLT) mean(state8scLT); mean(state9spLT) mean(state9scLT); mean(state10spLT) mean(state10scLT);];
    error=[std(state1spLT)/sqrt(numel(state1spLT)) std(state1scLT)/sqrt(numel(state1scLT)); std(state2spLT)/sqrt(numel(state2spLT)) std(state2scLT)/sqrt(numel(state2scLT)); std(state3spLT)/sqrt(numel(state3spLT)) std(state3scLT)/sqrt(numel(state3scLT)); std(state4spLT)/sqrt(numel(state4spLT)) std(state4scLT)/sqrt(numel(state4scLT)); std(state5spLT)/sqrt(numel(state5spLT)) std(state5scLT)/sqrt(numel(state5scLT)); std(state6spLT)/sqrt(numel(state6spLT)) std(state6scLT)/sqrt(numel(state6scLT)); std(state7spLT)/sqrt(numel(state7spLT)) std(state7scLT)/sqrt(numel(state7scLT)); std(state8spLT)/sqrt(numel(state8spLT)) std(state8scLT)/sqrt(numel(state8scLT)); std(state9spLT)/sqrt(numel(state9spLT)) std(state9scLT)/sqrt(numel(state9scLT)); std(state10spLT)/sqrt(numel(state10spLT)) std(state10scLT)/sqrt(numel(state10scLT))];
    h=barwitherr(error, ydata,1);%downloaded function called barwitherr, because I otherwise could not plot errors in correct location
    set(gca,'XTickLabel',{'1', '2', '3', '4','5', '6', '7','8','9', '10'})
    legend({'rrMDD','Control'}, 'Box', 'off', 'FontSize', 10, 'location','northeast')
    %ylabel('Lifetime', 'Fontsize',10) % y-axis label
    set(h(1),'FaceColor',[1 0.8 0]);
    set(h(2),'FaceColor',[0.6 0 1]);
    title('Sad mood duration')
     text(6.8,2.5,'*','FontSize',22)
     text(9.8,2.5,'*','FontSize',22)
     text(0.8,8.2,'+','FontSize',22)
     text(2.8,3.3,'+','FontSize',22)
     text(3.8,3.2,'+','FontSize',22)
     text(7,2.5,'+','FontSize',22)
     text(10,2.5,'+','FontSize',22)
    box off

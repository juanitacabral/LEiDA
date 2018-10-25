
load LEiDA_k_results Kmeans_results
n_Task=2;
n_Subjects=86;
Tmax=210;
TR=2;

%for 2 to 20 cluster, and for both tasks and all subjects, divide the time
%courses by group and task. 

%For every subject and task calculate probability of being in a given pattern. 
P=zeros(19,n_Subjects,2,20);

%For every subject and task calculate the mean lifetime of each pattern. 
LT=zeros(19,n_Subjects,2,20);

%For every subject and task calculate the switching probability. 
SwitchFreq=zeros(19,n_Subjects,2);

for k=2:20
        for task=1:n_Task
            for s=1:n_Subjects
                   
                if s<52 && task==1                      
                  T=Tmax*(s-1)*2+1:Tmax*s+Tmax*(s-1);
                elseif s<52 && task==2    
                   T=Tmax*s+Tmax*(s-1)+1:Tmax*s+Tmax*(s-1)+Tmax;        
                elseif s>51 && task==1       
                   T=Tmax*(s-1)*2+1:Tmax*s+Tmax*(s-1);
                elseif s>51 && task==2        
                   T=Tmax*s+Tmax*(s-1)+1:Tmax*s+Tmax*(s-1)+Tmax;       
                end
                          
                Ctime=Kmeans_results{k}.IDX(T);
                SwitchFreq(k-1,s,task)=mean(diff(Ctime)~=0)/TR;
                
                 for c=1:k
                % Probability
                P(k-1,s,task,c)=mean(Ctime==c);
             
                % Lifetimes
                Ctime_bin=Ctime==c;
                
                % Detect switches to (a) and from (b) state c
                a=find(diff(Ctime_bin)==1);
                b=find(diff(Ctime_bin)==-1);
                
                % We discard the cases where state sarts or ends ON
                if length(b)>length(a)
                    b(1)=[];
                elseif length(a)>length(b)
                    a(end)=[];
                elseif  ~isempty(a) && ~isempty(b) && a(1)>b(1)
                    b(1)=[];
                    a(end)=[];
                end
                
                if ~isempty(a) && ~isempty(b)            
                    C_Durations=b-a;
                else
                    C_Durations=0;     
                end              
                LT(k-1,s,task,c)=mean(C_Durations);
                LT(isnan(LT))=0;
                 end
            end
        end
end


%calculate probability of state occurence, differences between patients and
%controls and mood states.

pvalue_NvsS=zeros(2,19,20);
pvalue_PvsC=zeros(2,19,20);

LT_pval_NvsS=zeros(2,19,20);
LT_pval_PvsC=zeros(2,19,20);

Switch_pval_NvsS=zeros(2,19);
Switch_pval_PvsC=zeros(2,19);

for k=2:20
    disp(['Now running for ' num2str(k) ' clusters'])
    for c=1:k
        for t=1:2   % For each Mood, search differences between patients and controls
            a=P(k-1,1:51,t,c);  % Vector containing Prob of c in patients during mood t
            b=P(k-1,52:end,t,c); % Vector containing Prob of c in controls during mood t
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');          
            pvalue_PvsC(t,k-1,c)=min(stats.pvals);


            % SAME FOR LIFETIMES
            a=LT(k-1,1:51,t,c);  % Vector containing Prob of c in patients during mood t
            b=LT(k-1,52:end,t,c); % Vector containing Prob of c in controls during mood t
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');          
            LT_pval_PvsC(t,k-1,c)=min(stats.pvals);
            
            
            % SAME FOR Switchfreq
            a=SwitchFreq(k-1,1:51,t);  % Vector containing Switchfreq of c in Neutral mood in group s
            b=SwitchFreq(k-1,52:end,t);  % Vector containing Switchfreq of c in Sad mood in group s
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');          
            Switch_pval_PvsC(t,k-1)=min(stats.pvals);
                               
        end
        for group=1:2 % For each group, search differences between neutal or sad
            if group==1
                s=1:51;
            else
                s=52:86;
            end
            a=P(k-1,s,1,c);  % Vector containing Prob of c in Neutral mood in group s
            b=P(k-1,s,2,c);  % Vector containing Prob of c in Sad mood in group s
            stats=permutation_htest_np_paired([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');%paired sample T-test apadted from two sample T-test.          
            pvalue_NvsS(group,k-1,c)=min(stats.pvals);
            
            % SAME FOR LIFETIMES
            a=LT(k-1,s,1,c);  % Vector containing LT of c in Neutral mood in group s
            b=LT(k-1,s,2,c);  % Vector containing LT of c in Sad mood in group s
            stats=permutation_htest_np_paired([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');          
            LT_pval_NvsS(group,k-1,c)=min(stats.pvals);
    
            
            % SAME FOR Switchfreq
            a=SwitchFreq(k-1,s,1);  % Vector containing Switchfreq in Neutral mood in group s
            b=SwitchFreq(k-1,s,2);  % Vector containing Switchfreq in Sad mood in group s
            stats=permutation_htest_np_paired([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');          
            Switch_pval_NvsS(group,k-1)=min(stats.pvals);
end

end
end


% First correct p-values by number of K
% So we don't need to correct anymore by k in the following codes
for k=2:20
    pvalue_PvsC(:,k-1,:)=pvalue_PvsC(:,k-1,:)*k;
    pvalue_NvsS(:,k-1,:)=pvalue_NvsS(:,k-1,:)*k; 
    pvalue_NvsS2(:,k-1,:)=pvalue_NvsS2(:,k-1,:)*k;
    

    %%% Same for Lifetimes
    LT_pval_PvsC(:,k-1,:)=LT_pval_PvsC(:,k-1,:)*k;
    LT_pval_NvsS(:,k-1,:)=LT_pval_NvsS(:,k-1,:)*k;
    LT_pval_NvsS2(:,k-1,:)=LT_pval_NvsS2(:,k-1,:)*k;
    
    %%% We dont do this for Switching frequency
   
    
end

%% Count number of significant p_values:

% Probabilities:

P_PvsC=numel(find(pvalue_PvsC>0 & pvalue_PvsC<0.05));
disp(['Patients and Controls differ in prob of ' num2str(P_PvsC) ' FC patterns']) 
P_NvsS=numel(find(pvalue_NvsS>0 & pvalue_NvsS<0.05));
disp(['Neutral and Sad mood differ in prob of ' num2str(P_NvsS) ' FC patterns']) 

% Lifetimes

LT_PvsC=numel(find(LT_pval_PvsC>0 & LT_pval_PvsC<0.05));
disp(['Patients and Controls differ in lifetime of ' num2str(LT_PvsC) ' FC patterns']) 
LT_NvsS=numel(find(LT_pval_NvsS2>0 & LT_pval_NvsS2<0.05));
disp(['Neutral vs Sad differ in lifetime of ' num2str(LT_NvsS) ' FC patterns']) 

%switching frequency
Switch_PvsC=numel(find(Switch_pval_PvsC>0 & Switch_pval_PvsC<0.05));
disp(['Patients and Controls differ in frequency of ' num2str(Switch_PvsC) 'clustering solutions']) 
Switch_NvsS=numel(find(Switch_pval_NvsS>0 & Switch_pval_NvsS<0.05));
disp(['Neutral vs Sad differ in lifetime of ' num2str(Switch_NvsS) ' clustering solution'])

save ClusterStats pvalue_NvsS pvalue_PvsC P LT_pval_NvsS LT_pval_PvsC LT
% Note that the stats saved are already corrected by the number k

load AAL_labels.mat label90
Order=[1:2:89 90:-2:2];

% Patients vs Controls, plot most significant networks
%these are the right values that we calculated with 10.000 permutations
%after revision for HMB:

% load Between_results_sept_2018 pvalue_PvsC LT_pval_PvsC P LT
% load Within_results_okt_2018 LT_pval_NvsS pvalue_NvsS 


%For probability
for t=1:2
    if t==1
        task='Neutral';
    else
        task='Sad';
    end
    for k=2:20
        for c=1:k
            if pvalue_PvsC(t,k-1,c)<(0.05) && pvalue_PvsC(t,k-1,c)>0
                V=Kmeans_results{k}.C(c,:);
                figure('Name',['Patients vs Controls under ' task ' mood K=' num2str(k) ' c=' num2str(c) ])
                subplot(2,3,1)
                plot_nodes_in_cortex(V);
                subplot(2,3,[2 5])
                barh(V(Order),'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
                set(gca,'YTick',1:90,'Fontsize',6)
                set(gca,'YTickLabel',label90(Order,:))
                ylim([0 91])
                xlim([-.15 .15])
                grid on
                subplot(2,3,3)
                pat=P(k-1,1:51,t,c);
                con=P(k-1,52:end,t,c);
                bar([mean(pat) mean(con)])
                hold on
                errorbar([mean(pat) mean(con)],[std(pat)/sqrt(numel(pat)) std(con)/sqrt(numel(con))],'LineStyle','none')
                set(gca,'XTickLabel',{'Patients', 'Controls'})
                ylabel('Probability')
                box off
                subplot(2,3,4)
                FC=V'*V;
                imagesc(FC(Order,Order))
                axis square
            end
        end
    end
end

% Neutral vs Sad 
%probability plot most significant networks
for s=1:2
    if s==1
        group='Patients';
        sample=1:51;
    else
        group='Controls';
        sample=52:86;
    end
    for k=2:20

        for c=1:k
            if pvalue_NvsS(s,k-1,c)<(0.05) && pvalue_NvsS(s,k-1,c)>0
                V=Kmeans_results{k}.C(c,:);
                figure('Name',['Neutral vs Sad for ' group ' K=' num2str(k) ' c=' num2str(c) ])
                subplot(2,3,1)
                plot_nodes_in_cortex(V);
                subplot(2,3,[2 5])
                barh(V(Order),'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
                set(gca,'YTick',1:90,'Fontsize',6)
                set(gca,'YTickLabel',label90(Order,:))
                ylim([0 91])
                xlim([-.15 .15])
                grid on
                subplot(2,3,3)
                neutral=P(k-1,sample,1,c);
                sad=P(k-1,sample,2,c);
                bar([mean(neutral) mean(sad)])
                hold on
                errorbar([mean(neutral) mean(sad)],[std(neutral)/sqrt(numel(neutral)) std(sad)/sqrt(numel(sad))],'LineStyle','none')
                set(gca,'XTickLabel',{'Neutral', 'Sad'})
                ylabel('Probability')
                box off
                subplot(2,3,4)
                FC=V'*V;
                imagesc(FC(Order,Order))
                axis square
            end
        end
    end
end


%pt vs controls
%LT plot most significant networks
for t=1:2
    if t==1
        task='Neutral';
    else
        task='Sad';
    end
    for k=2:20
        for c=1:k
            if LT_pval_PvsC(t,k-1,c)<(0.05) && LT_pval_PvsC(t,k-1,c)>0
                V=Kmeans_results{k}.C(c,:);
                figure('Name',['Patients vs Controls under ' task ' mood K=' num2str(k) ' c=' num2str(c) ])
                subplot(2,3,1)
                plot_nodes_in_cortex(V);
                subplot(2,3,[2 5])
                barh(V(Order),'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
                set(gca,'YTick',1:90,'Fontsize',6)
                set(gca,'YTickLabel',label90(Order,:))
                ylim([0 91])
                xlim([-.15 .15])
                grid on
                subplot(2,3,3)
                pat=LT(k-1,1:51,t,c);
                con=LT(k-1,52:end,t,c);
                bar([nanmean(pat) nanmean(con)])
                hold on
                errorbar([mean(pat) mean(con)],[std(pat)/sqrt(numel(pat)) std(con)/sqrt(numel(con))],'LineStyle','none')
                set(gca,'XTickLabel',{'Patients', 'Controls'})
                ylabel('Lifetime')
                box off
                subplot(2,3,4)
                FC=V'*V;
                imagesc(FC(Order,Order))
                axis square
            end
        end
    end
end

% Neutral vs Sad
%LT plot most significant networks

for s=1:2
    if s==1
        group='Patients';
        sample=1:51;
    else
        group='Controls';
        sample=52:86;
    end
    for k=2:20

        for c=1:k
            if LT_pval_NvsS(s,k-1,c)<(0.05) && LT_pval_NvsS(s,k-1,c)>0
                V=Kmeans_results{k}.C(c,:);
                figure('Name',['Neutral vs Sad lifetime for ' group ' K=' num2str(k) ' c=' num2str(c) ])
                subplot(2,3,1)
                plot_nodes_in_cortex(V);
                subplot(2,3,[2 5])
                barh(V(Order),'FaceColor',[0.4 0.4 0.4],'EdgeColor','none')
                set(gca,'YTick',1:90,'Fontsize',6)
                set(gca,'YTickLabel',label90(Order,:))
                ylim([0 91])
                xlim([-.15 .15])
                grid on
                subplot(2,3,3)
                neutral=LT(k-1,sample,1,c);
                sad=LT(k-1,sample,2,c);
                bar([mean(neutral) mean(sad)])
                hold on
                errorbar([mean(neutral) mean(sad)],[std(neutral)/sqrt(numel(neutral)) std(sad)/sqrt(numel(sad))],'LineStyle','none')
                set(gca,'XTickLabel',{'Neutral', 'Sad'})
                ylabel('Lifetime')
                box off
                subplot(2,3,4)
                FC=V'*V;
                imagesc(FC(Order,Order))
                axis square
            end
        end
    end
end

%Switchfreq patients vs controls
for t=1:2
    if t==1
        task='Neutral';
    else
        task='Sad';
    end
    for k=2:20
            if Switch_pval_PvsC(t,k-1)<(0.05) && Switch_pval_PvsC(t,k-1)>0
                figure('Name',['Patients vs Controls under ' task ' mood K=' num2str(k) ' c=' num2str(c) ])
                subplot(1,1,1)
                pat=SwitchFreq(k-1,1:51,t);
                con=SwitchFreq(k-1,52:end,t);
                bar([nanmean(pat) nanmean(con)])
                hold on
                errorbar([mean(pat) mean(con)],[std(pat)/sqrt(numel(pat)) std(con)/sqrt(numel(con))],'LineStyle','none')
                set(gca,'XTickLabel',{'Patients', 'Controls'})
                ylabel('Switching frequency')
                box off
            end
        end
end 

%Switchfreq neutral vs sad

for s=1:2
    if s==1
        group='Patients';
        sample=1:51;
    else
        group='Controls';
        sample=52:86;
    end
    for k=2:20
            if Switch_pval_NvsS(s,k-1)<(0.05) && Switch_pval_NvsS(s,k-1)>0
                figure('Name',['Patients vs Controls under ' task ' mood K=' num2str(k) ' c=' num2str(c) ])
                subplot(1,1,1)
                Neut=SwitchFreq(k-1,1:51,t);
                con=SwitchFreq(k-1,52:end,t);
                bar([nanmean(pat) nanmean(con)])
                hold on
                errorbar([mean(pat) mean(con)],[std(pat)/sqrt(numel(pat)) std(con)/sqrt(numel(con))],'LineStyle','none')
                set(gca,'XTickLabel',{'Patients', 'Controls'})
                ylabel('Switching frequency')
                box off
            end
        end
end


%make figure for switching frequencies
for k=2:20
      meanFreqpatneut(k)=mean(SwitchFreq(k-1,1:52,1));
      meanFreqcontneut(k)=mean(SwitchFreq(k-1,52:end,1)); 
      meanFreqpatsad(k)=mean(SwitchFreq(k-1,1:52,2));
      meanFreqcontsad(k)=mean(SwitchFreq(k-1,52:end,2));
end

figure
subplot(2,2,[1 2])
bar([meanFreqpatneut' meanFreqcontneut'])
title('Switching frequency neutral mood')
legend('remitted-MDD','controls')
ylabel('Frequency (Hz)')
xlabel('K')

subplot(2,2,[3 4])
bar([meanFreqpatsad' meanFreqcontsad'])
title('Switching frequency sad mood')
legend('remitted-MDD','controls')
ylabel('Frequency (Hz)')
xlabel('K')
 

% FIGURE SHOWING ALL CLUSTERING SOLUTIONS INTO K FC PATTERNS
% FC patterns are sorted according to their probability of occurence

% Add symbol for the FC states that are significantly different 
% (significance is p<0.05/K

% * Patients vs Controls (NEUTRAL)
% + Patients vs Controls (SAD)
% o Neutral vs Sad (PATIENTS)
% # Neutral vs Sad (CONTROLS)

figure
colormap(jet);
for k=2:20
    [~, ind_sort]=sort(Kmeans_results{k}.SUMD,'descend');   
    % I added this correction to sort each state according to its
    % probability of occurrence
    for c=1:k    
      V=Kmeans_results{k}.C(ind_sort(c),:);
      u=sub2ind([20 19],c,k-1);
      subplot(19,20,u)
            FC=V'*V;
            imagesc(FC(Order,Order))
            %axis square
            set(gca,'YTick',[])
            set(gca,'XTick',[])
            Add_Symbol=[];
            if pvalue_PvsC(1,k-1,ind_sort(c))<(0.05) && pvalue_PvsC(1,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '*p'];  % ONLY IN NEUTRAL MOOD
            end
            if pvalue_PvsC(2,k-1,ind_sort(c))<(0.05) && pvalue_PvsC(2,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '+p']; % ONLY IN SAD MOOD
            end
            if pvalue_NvsS(1,k-1,ind_sort(c))<(0.05) && pvalue_NvsS(1,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '#p'];  % ONLY IN PATIENTS
            end
            if pvalue_NvsS(2,k-1,ind_sort(c))<(0.05) && pvalue_NvsS(2,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '!p']; % ONLY IN CONTROLS
            end
             if LT_pval_PvsC(1,k-1,ind_sort(c))<(0.05) && pvalue_PvsC(1,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '*L'];  % ONLY IN NEUTRAL MOOD
            end
            if LT_pval_PvsC(2,k-1,ind_sort(c))<(0.05) && pvalue_PvsC(2,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '+L']; % ONLY IN SAD MOOD
            end
            if LT_pval_NvsS(1,k-1,ind_sort(c))<(0.05) && LT_pval_NvsS(1,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '#L'];  % ONLY IN PATIENTS
            end
            if LT_pval_NvsS(2,k-1,ind_sort(c))<(0.05) && LT_pval_NvsS(2,k-1,ind_sort(c))>0
            Add_Symbol=[Add_Symbol '!L']; % ONLY IN CONTROLS
            end
            ylabel(Add_Symbol)
%             xlabel('Number of Clusters (K)') % x-axis label 

    end
end


% Select the Eigenvectors with significant differences between controls and patients
% In both mood conditions

Sig_Neutral=squeeze(pvalue_PvsC(1,:,:));
Sig_Sad=squeeze(pvalue_PvsC(2,:,:));

Sign_neutral_LT=squeeze(LT_pval_PvsC(1,:,:));
Sign_sad_LT=squeeze(LT_pval_PvsC(2,:,:));

Sig_NvsS_P= squeeze(pvalue_NvsS(1,:,:));
Sig_NvsS_LT= squeeze(LT_pval_NvsS(1,:,:));
% Find indices of the networks that are significant in both conditions
[k_sig, c_sig]=find(Sig_Neutral<0.05 & Sig_Neutral>0 & Sig_Sad<0.05 & Sig_Sad>0);
[neutral_sig, cneutral_sig]=find(Sig_Neutral<0.05 & Sig_Neutral>0);

figure('Name','Patients vs Controls')
colormap(jet)
for n=1:length(k_sig)
    subplot(2,length(k_sig),n)
    V(n,:)=Kmeans_results{k_sig(n)+1}.C(c_sig(n),:);
    FC=V(n,:)'*V(n,:);
    imagesc(FC(Order,Order))
    title(['One of ' num2str(k_sig(n)+1) ' clusters']) 
    %xlabel(['p_N=' num2str(Sig_Neutral(x(n),y(n))) ' p_S=' num2str(Sig_Sad(x(n),y(n))) ])
    subplot(2,length(k_sig),n+length(k_sig))
    plot_nodes_in_cortex(V(n,:));    
    
end
% We save a 90x1 vector containing the pattern that significantly
% differs between patients and controls in both neutral conditions
save Vector_PvsC V

for n=1:length(neutral_sig)
    subplot(2,length(neutral_sig),n)
    V(n,:)=Kmeans_results{neutral_sig(n)+1}.C(cneutral_sig(n),:);
    FC=V(n,:)'*V(n,:);
    imagesc(FC(Order,Order))
    %title(['One of ' num2str(k_sig(n)+1) ' clusters']) 
    %xlabel(['p_N=' num2str(Sig_Neutral(x(n),y(n))) ' p_S=' num2str(Sig_Sad(x(n),y(n))) ])
    subplot(2,length(neutral_sig),n+length(neutral_sig))
    plot_nodes_in_cortex(V(n,:));    
end


VC= corrcoef(V'); % All correlate above 0.84. Same network. 

%to calculate correlateions between other networks:

a=Kmeans_results{4}.C(1,:);
b=Kmeans_results{5}.C(1,:);
c=Kmeans_results{6}.C(6,:);
d=Kmeans_results{14}.C(6,:);

A=[a' b' c' d']

VD=corrcoef(A)

 
% to look up SEM and mean:
% 
% std(LT(9,1:51,1,8))/sqrt(numel(LT(9,1:51,1,8)))
% std(LT(9,52:end,1,8))/sqrt(numel(LT(9,52:end,1,8)))
% std(P(9,1:51,1,8))/sqrt(numel(P(9,1:51,1,8)))
% std(P(9,52:end,1,8))/sqrt(numel(P(9,52:end,1,8)))
% 
% I can switch colors in plot_nodes_in_cortex
% 
% std(SwitchMatrix(1:51,1,5,4)/sqrt(numel(SwitchMatrix(1:51,1,5,4))
% p_NvsS_pat= squeeze(pvalue_(1,:,:));

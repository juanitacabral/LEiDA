%function FCSwitching

%%%%%%%%%%%%%%%%%%%%%%
%
%  Calculates the Switching matrix for each subject and each condition
%  Plots the mean switching matrix for each condition
%  Computes statistics and add asterics to the switches that are
%  significantly different between conditions.
%
%%%%%%%%%%%%%%%%%%%%%%%%



load LEiDA_k_results Kmeans_results 

n_Task=2;
n_Subjects=86;
Tmax=210;
TR = 2;

k=10;

[~, ind_sort]=sort(Kmeans_results{k}.SUMD,'descend');

SwitchMatrix=zeros(n_Subjects,2,k,k);
 
Ctime=zeros(1,Tmax);

% Select the time interval per subject and task

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
        
        Ctime(1)=find(ind_sort==Kmeans_results{k}.IDX(T(1))); %sort according to probability of occurence       
              
        % for every subject we now have the 210 timepoints(TR) for this subject (in total we have 210*86*2 timepoints)
        %ctime is which cluster is active for every time point (210)
        
        for t=2:length(T) %for T=2:210
            Ctime(t)=find(ind_sort==Kmeans_results{k}.IDX(T(t)));
            if Ctime(t)~= Ctime(t-1) %if cluster at that timepoint is not the same as cluster at next timepoint
                SwitchMatrix(s,task,Ctime(t-1),Ctime(t))=SwitchMatrix(s,task,Ctime(t-1),Ctime(t))+1; %Switchmatrix of subject during task, display cluster of follow timepoint, cluster of this timpoint, switchmatrix gets one extra value
            end
            
        end       
    end
end

%meanswitching for the rrMDD and controls for sad and neutral 
MeanSwitchpatneutral=squeeze(mean(squeeze(SwitchMatrix(1:51,1,:,:))));
MeanSwitchconneutral=squeeze(mean(squeeze(SwitchMatrix(52:end,1,:,:))));

MeanSwitchpatsad=squeeze(mean(squeeze(SwitchMatrix(1:51,2,:,:))));
MeanSwitchconsad=squeeze(mean(squeeze(SwitchMatrix(52:end,2,:,:))));

%for figure 6
% here I also want to average the average of sad and neutral for the groups
% seperately

%for supplementary figure
% here I also want to average the average of sad and neutral for the groups
% seperately
% First I concatenate the switch matrices in neutral and sad mood, and only
% average after

MeanSwitchTotalPat=squeeze(mean(cat(1,squeeze(SwitchMatrix(1:51,1,:,:)),squeeze(SwitchMatrix(1:51,2,:,:)))));
MeanSwitchTotalCon=squeeze(mean(cat(1,squeeze(SwitchMatrix(52:end,1,:,:)),squeeze(SwitchMatrix(52:end,2,:,:)))));

%total neutral and sad averaged over both group
MeanSwitchTotalneutral=squeeze(mean(squeeze(SwitchMatrix(:,1,:,:))));
MeanSwitchTotalsad=squeeze(mean(squeeze(SwitchMatrix(:,2,:,:))));

% To consider the probability of switching instead of number of switches

for c=1:k
        MeanSwitchpatneutral(c,:) = MeanSwitchpatneutral(c,:)/sum(MeanSwitchpatneutral(c,:));  
        MeanSwitchpatsad(c,:)     = MeanSwitchpatsad(c,:)/sum(MeanSwitchpatsad(c,:));   
        MeanSwitchconneutral(c,:) = MeanSwitchconneutral(c,:)/sum(MeanSwitchconneutral(c,:));  
        MeanSwitchconsad(c,:)     = MeanSwitchconsad(c,:)/sum(MeanSwitchconsad(c,:)); 
        MeanSwitchTotalneutral(c,:) = MeanSwitchTotalneutral(c,:)/sum(MeanSwitchTotalneutral(c,:));
        MeanSwitchTotalsad(c,:) = MeanSwitchTotalsad(c,:)/sum(MeanSwitchTotalsad(c,:));
        MeanSwitchTotalPat(c,:)= MeanSwitchTotalPat(c,:)/sum(MeanSwitchTotalPat(c,:));
        MeanSwitchTotalCon(c,:)= MeanSwitchTotalCon(c,:)/sum(MeanSwitchTotalCon(c,:));
end

% PLOT SWITCHING MATRICES

% Find the maximum value of all matrices to scale all colorbars
max_value=max(cat(1,MeanSwitchpatneutral(:),MeanSwitchpatsad(:),MeanSwitchconneutral(:),...
    MeanSwitchpatneutral(:)));




%figure S7
figure
colormap(jet)
subplot(2,2,1)
imAlpha=ones(size(MeanSwitchpatneutral));
imAlpha(MeanSwitchpatneutral==0)=0;
imagesc(MeanSwitchpatneutral,'AlphaData',imAlpha, [0 max_value]);
title('rrMDD Neutral','FontSize',14, 'FontWeight','Bold')  
ylabel('From FC pattern')
xlabel('To FC pattern')
axis square
colorbar


subplot(2,2,2)
imAlpha=ones(size(MeanSwitchconneutral));
imAlpha(MeanSwitchconneutral==0)=0;
imagesc(MeanSwitchconneutral,'AlphaData',imAlpha, [0 max_value]);
title('Control neutral','FontSize',14, 'FontWeight','Bold')   
axis square
ylabel('From FC pattern')
xlabel('To FC pattern')
colorbar 

subplot(2,2,3)

imAlpha=ones(size(MeanSwitchpatsad));
imAlpha(MeanSwitchpatsad==0)=0;
imagesc(MeanSwitchpatsad,'AlphaData',imAlpha, [0 max_value]);
title('Remitted-MDD Sad','FontSize',14, 'FontWeight','Bold')  
axis square
ylabel('From FC pattern')
xlabel('To FC pattern')
colorbar

subplot(2,2,4)
imAlpha=ones(size(MeanSwitchconsad));
imAlpha(MeanSwitchconsad==0)=0;
imagesc(MeanSwitchconsad,'AlphaData',imAlpha, [0 max_value]);
title('Control Sad','FontSize',14, 'FontWeight','Bold')  
axis square
ylabel('From FC pattern')
xlabel('To FC pattern')
colorbar



for c_out=1:k
    for c_in=1:k

        %a=squeeze(mean(squeeze(SwitchMatrix(:,1,:,:))));
        %b=squeeze(mean(squeeze(SwitchMatrix(:,2,:,:))));
        
        % Patients vs Controls, task 1 (neutral)
        a=SwitchMatrix(1:51,1,c_out,c_in)';   % Patients neutral
        b=SwitchMatrix(52:end,1,c_out,c_in)'; % Controls neutral
        c=SwitchMatrix(1:51,2,c_out,c_in)';   % Patients Sad
        d=SwitchMatrix(52:end,2,c_out,c_in)'; % controls Sad
        
        stats_ab=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');
        Switch_pval_ab(c_out,c_in)=min(stats_ab.pvals);
  
        stats_ac=permutation_htest_np_paired([a,c],[ones(1,numel(a)) 2*ones(1,numel(c))],10000,0.05,'ttest');
        Switch_pval_ac(c_out,c_in)=min(stats_ac.pvals);    
  
            
        stats_bd=permutation_htest_np_paired([b,d],[ones(1,numel(b)) 2*ones(1,numel(d))],10000,0.05,'ttest');
        Switch_pval_bd(c_out,c_in)=min(stats_bd.pvals);
        
        stats_cd=permutation_htest2_np([c,d],[ones(1,numel(c)) 2*ones(1,numel(d))],10000,0.05,'ttest');
        Switch_pval_cd(c_out,c_in)=min(stats_cd.pvals);
            
            if min(stats_ab.pvals) < 0.05 
                subplot(2,2,1)
                text(c_in-0.4,c_out,'*','FontSize',14, 'FontWeight','Bold','Color', 'white')
                w=text(c_in-0.2,c_out,num2str(mean(MeanSwitchpatneutral(c_out,c_in)),'% 10.2f'))
                set(w, 'FontSize', 11, 'Color', 'white')
                subplot(2,2,2)
                text(c_in-0.4,c_out,'*','FontSize',14, 'FontWeight','Bold', 'Color', 'white')
                y=text(c_in-0.2,c_out,num2str(mean(MeanSwitchconneutral(c_out,c_in)),'% 10.2f'))
                set(y, 'FontSize', 11, 'Color', 'white')
            end
           
            if min(stats_ac.pvals) < 0.05 
                subplot(2,2,1)
                text(c_in-.4,c_out+0.1,'+','FontSize',8, 'FontWeight','Bold', 'Color', 'r')
            end
            if min(stats_ac.pvals) < 0.05 && min(stats_ab.pvals) > 0.05 
                subplot(2,2,1)
                z=text(c_in-0.2,c_out,num2str(mean(MeanSwitchpatneutral(c_out,c_in)),'% 10.2f'))
                set(z, 'FontSize', 11, 'Color', 'white')
            end
            if min(stats_ac.pvals) < 0.05 && min(stats_ac.pvals) > 0 
                subplot(2,2,3)
                text(c_in-.4,c_out,'+','FontSize',10, 'FontWeight','Bold', 'Color', 'r')
                e=text(c_in,c_out,num2str(mean(MeanSwitchpatsad(c_out,c_in)),'% 10.2f'))
                set(e, 'FontSize', 11, 'Color', 'white')
            end
             if min(stats_bd.pvals) < 0.05 && min(stats_ac.pvals) > 0
                subplot(2,2,2)
                text(c_in-.4,c_out,'+','FontSize',10, 'FontWeight','Bold', 'Color', 'white')
                v=text(c_in-0.2,c_out,num2str(mean(MeanSwitchconneutral(c_out, c_in)),'% 10.2f'))
                set(v, 'FontSize', 11, 'Color', 'white')
                subplot(2,2,4)
                text(c_in-.4,c_out,'+','FontSize',10, 'FontWeight','Bold', 'Color', 'white')
                f=text(c_in-0.2,c_out,num2str(mean(MeanSwitchconsad(c_out, c_in)),'% 10.2f'))
                set(f, 'FontSize', 11, 'Color', 'white')
             end    
             if min(stats_cd.pvals) < 0.05
                subplot(2,2,3)
                text(c_in-0.2,c_out,'*','FontSize',10, 'FontWeight','Bold','Color', 'r')
             end
             if min(stats_cd.pvals) < 0.05 %&& min(stats_ac.pvals) > 0.05
                subplot(2,2,3)
                g=text(c_in,c_out,num2str(mean(MeanSwitchpatsad(c_out,c_in)),'% 10.2f'))
                set(g, 'FontSize', 11, 'Color', 'white')
             end
                if min(stats_cd.pvals) < 0.05 
                subplot(2,2,4)
                text(c_in-.4,c_out+0.3,'*','FontSize',8, 'FontWeight','Bold', 'Color', 'r')
                end
                if min(stats_cd.pvals) < 0.05  %&& min(stats_bd.pvals) >0.05
                h=text(c_in-0.2,c_out,num2str(mean(MeanSwitchconsad(c_out, c_in)),'% 10.2f'))
                set(h, 'FontSize', 11, 'Color', 'white')
             end  
    
    end
end

%to correct for the number of states
Switchcorrected_ab=Switch_pval_ab*10
Switchcorrected_ac=Switch_pval_ac*10
Switchcorrected_bd=Switch_pval_bd*10
Switchcorrected_cd=Switch_pval_cd*10

%figure 6 for total
figure (2)
colormap(jet)
subplot(2,2,1)
imAlpha=ones(size(MeanSwitchTotalneutral));
imAlpha(MeanSwitchTotalneutral==0)=0;
imagesc(MeanSwitchTotalneutral,'AlphaData',imAlpha, [0 0.6]);
title('Whole group Neutral','FontSize',14, 'FontWeight','Bold')  
ylabel('From FC pattern')
xlabel('To FC pattern')
axis square
colorbar


subplot(2,2,2)
imAlpha=ones(size(MeanSwitchTotalsad));
imAlpha(MeanSwitchTotalsad==0)=0;
imagesc(MeanSwitchTotalsad,'AlphaData',imAlpha, [0 0.6]);
title('Whole group sad','FontSize',14, 'FontWeight','Bold')   
axis square
ylabel('From FC pattern')
xlabel('To FC pattern')
colorbar 

%figure S8 for total pat and controls
figure (3)
colormap(jet)
subplot(2,2,1)
imAlpha=ones(size(MeanSwitchTotalPat));
imAlpha(MeanSwitchTotalPat==0)=0;
imagesc(MeanSwitchTotalPat,'AlphaData',imAlpha, [0 0.6]);
title('rrMDD mean mood','FontSize',14, 'FontWeight','Bold')  
ylabel('From FC pattern')
xlabel('To FC pattern')
axis square
colorbar


subplot(2,2,2)
imAlpha=ones(size(MeanSwitchTotalCon));
imAlpha(MeanSwitchTotalCon==0)=0;
imagesc(MeanSwitchTotalCon,'AlphaData',imAlpha, [0 0.6]);
title('Control mean mood','FontSize',14, 'FontWeight','Bold')   
axis square
ylabel('From FC pattern')
xlabel('To FC pattern')
colorbar 



%to calculate SE
%std(SwitchMatrix(52:end,1,5,4)/sqrt(numel(SwitchMatrix(52:end,1,5,4))))
%0.1767/(sum(SwitchMatrix(52:end,1,5,4)))



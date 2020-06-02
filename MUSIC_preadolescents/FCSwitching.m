function FCSwitching

%%%%%%%%%%%%%%%%%%%%%%
%
%  Calculates the Switching matrix for each subject and each condition
%  Plots the mean switching matrix for each condition
%  Computes statistics and add asterics to the switches that are
%  significantly different between conditions.
%
%%%%%%%%%%%%%%%%%%%%%%%%

load LEiDA_results_v6.mat Kmeans_results Time_all rangeK
load T_paradigm.mat T_paradigm

k=7;
n_Subjects=17;

IDX=Kmeans_results{rangeK==k}.IDX;
clear Kmeans_results

% We save one kxk switching matrix for each subject and each task
SwitchMatrix=zeros(n_Subjects,2,k,k);

% SILENCE SWITCHING MATRICES
for task=1:2
    
    disp([ 'Task ' num2str(task)])
    
    for s=1:n_Subjects
        
        disp([ 'Subject ' num2str(s)])
        
        Ctime=IDX(Time_all==s);
        
        if task ==1  % Detence Silence Epochs
            T_task=find(T_paradigm==0);
        elseif task ==2 % Detence Music Epochs
            T_task=find(T_paradigm>0);
        end   
        
        % Select the time intervals representing SILENCE or MUSIC
        
        for tr=2:length(T_task)
            t=T_task(tr);
            if (T_task(tr)-T_task(tr-1))==1
                if Ctime(t)~= Ctime(t-1)
                    disp(['Detected a switch from state ' num2str(Ctime(t-1)) ' to ' num2str(Ctime(t)) ' at time ' num2str(t)])
                    SwitchMatrix(s,task,Ctime(t-1),Ctime(t))=SwitchMatrix(s,task,Ctime(t-1),Ctime(t))+1;
                    %SwitchMatrix(s,task,Ctime(t),Ctime(t-1))=SwitchMatrix(s,task,Ctime(t),Ctime(t-1))+1;
                end
            end
        end
        for c=1:k
            if sum(SwitchMatrix(s,task,c,:))>0
                SwitchMatrix(s,task,c,:)=SwitchMatrix(s,task,c,:)/sum(SwitchMatrix(s,task,c,:));
            else
                SwitchMatrix(s,task,c,:)=0;
            end
        end
    end
end

MeanSwitchSilence=squeeze(mean(squeeze(SwitchMatrix(:,1,:,:))));
MeanSwitchMusic=squeeze(mean(squeeze(SwitchMatrix(:,2,:,:))));

% PLOT SWITCHING MATRICES

WhiteDiag=not(diag(ones(1,k)));


figure
colormap(jet)
subplot(1,3,1)
imagesc(MeanSwitchSilence,'AlphaData',WhiteDiag,[0 max(cat(1,MeanSwitchSilence(:),MeanSwitchMusic(:)))])
title({'SILENCE','Switching Probability'})
ylabel('From PC state')
xlabel('To PC state')
colorbar
axis square 

subplot(1,3,2)
imagesc(MeanSwitchMusic,'AlphaData',WhiteDiag,[0 max(cat(1,MeanSwitchSilence(:),MeanSwitchMusic(:)))])
title({'MUSIC','Switching Probability'})   
ylabel('From PC state')
xlabel('To PC state')
colorbar 
axis square 

subplot(1,3,3)
imagesc(MeanSwitchMusic-MeanSwitchSilence,'AlphaData',WhiteDiag)
title('Music - Silence')   
ylabel('From PC state')
xlabel('To PC state')
colorbar 

for c_out=1:k
    for c_in=1:k
        if not(c_out==c_in)
            a=SwitchMatrix(:,1,c_out,c_in)'; % Silence
            b=SwitchMatrix(:,2,c_out,c_in)'; % Music
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],10000,0.05,'ttest');
            SwitchSignif(c_out,c_in)=min(stats.pvals);
            % Add the mean value in the matrix
            % Round the values at 2 decimals:
            value=round(MeanSwitchSilence(c_out,c_in)*100)/100;
            subplot(1,3,1)
            hold on
            text(c_in-.3,c_out+.1,num2str(value))
            subplot(1,3,2)
            hold on
            value=round(MeanSwitchMusic(c_out,c_in)*100)/100;
            text(c_in-.3,c_out+.1,num2str(value))
            subplot(1,3,3)
            hold on
            text(c_in-.3,c_out+.1,num2str(SwitchSignif(c_out,c_in)))
            % When the switching is significantly different between groups,
            if SwitchSignif(c_out,c_in)<0.05
                subplot(1,3,1)
                text(c_in,c_out-.1,'**','Fontsize',14)
                subplot(1,3,2)
                text(c_in,c_out-.1,'**','Fontsize',14)
                subplot(1,3,3)
                text(c_in,c_out-.1,'**','Fontsize',14)
            elseif SwitchSignif(c_out,c_in)<0.1
                subplot(1,3,1)
                text(c_in,c_out-.1,'*','Fontsize',14)
                subplot(1,3,2)
                text(c_in,c_out-.1,'*','Fontsize',14)
                subplot(1,3,3)
                text(c_in,c_out-.1,'*','Fontsize',14)
            end
        end
    end
end

save SwitchingData SwitchSignif MeanSwitchMusic MeanSwitchSilence SwitchMatrix
  



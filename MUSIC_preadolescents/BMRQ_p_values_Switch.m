%function BMRQ_pvalues_Switch

load LEiDA_results_v6.mat P
load BMRQ_scores BMRQforanalysisOct
load SwitchingData

rangeK=3:15;
k=7;


for sc=1:size(BMRQforanalysisOct,2)
    
    scores=BMRQforanalysisOct(:,sc);
    for c_out=1:k
        for c_in=1:k
            
            if not(c_out==c_in)
                
                P_switch_Silence=SwitchMatrix(:,1,c_out,c_in);
                P_switch_Music=SwitchMatrix(:,2,c_out,c_in);
                
                P_switch_MusicRest=P_switch_Music-P_switch_Silence;
                
                [cc, p]=corrcoef(scores,P_switch_MusicRest);
                pval_Switch_BMRQ(sc,c_out,c_in)=p(2);
            end
            
        end
        
    end
end
  
figure
colormap(jet)
for sc=1:6
    subplot(2,3,sc)
    
    imagesc(MeanSwitchMusic-MeanSwitchSilence)
    title('Music - Silence')
    ylabel('From PC state')
    xlabel('To PC state')
    colorbar
    for c_out=1:k
        for c_in=1:k
            if not(c_out==c_in)   
                if pval_Switch_BMRQ(sc,c_out,c_in)<0.05
                    text(c_in,c_out-.1,'*','Fontsize',14)
                    text(c_in-.3,c_out+.1,num2str(pval_Switch_BMRQ(sc,c_out,c_in)))
                
                end
            end
        end
    end
    title(['BMRQ score' num2str(sc)])

end

load LEiDA_results_v6.mat Kmeans_results

% COMPARE BMRQ 6 with k=14 c=9
figure
sc=6;
c_out=3;
c_in=2;

subplot(1,3,1)
scores=BMRQforanalysisOct(:,sc);
P_switch_Silence=SwitchMatrix(:,1,c_out,c_in);
P_switch_Music=SwitchMatrix(:,2,c_out,c_in);
P_switch_MusicRest=P_switch_Music-P_switch_Silence;
plot(P_switch_MusicRest,scores,'*')
xlabel(['Diff Switch Prob from k=7 PC ' num2str(c_out) ' to ' num2str(c_in)])
ylabel(['BMRQ score ' num2str(sc)])
[cc, p]=corrcoef(scores,P_switch_MusicRest);
title(['cc= ' num2str(cc(2)) ' p=' num2str(p(2))])

Vc=Kmeans_results{rangeK==k}.C;

subplot(1,3,2)
plot_nodes_in_cortex(Vc(c_out,:))

subplot(1,3,3)
plot_nodes_in_cortex(Vc(c_in,:))
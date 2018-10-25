%for figure 1 to find when there are switches for the individual
%participants. 
load LEiDA_k_results Kmeans_results
n_Task=2;
n_Subjects=86;
Tmax=210;

%for 2 to 20 cluster, and for both tasks and all subjects, divide the time
%courses by group and task. 

k=10;
task=1;
% figure
%             for s=52:n_Subjects
%                    
%                 if s<52 && task==1                      
%                   T=Tmax*(s-1)*2+1:Tmax*s+Tmax*(s-1);
%                 elseif s<52 && task==2    
%                    T=Tmax*s+Tmax*(s-1)+1:Tmax*s+Tmax*(s-1)+Tmax;        
%                 elseif s>51 && task==1       
%                    T=Tmax*(s-1)*2+1:Tmax*s+Tmax*(s-1);
%                 elseif s>51 && task==2        
%                    T=Tmax*s+Tmax*(s-1)+1:Tmax*s+Tmax*(s-1)+Tmax;       
%                 end
%                         
%                 Ctime=Kmeans_results{k}.IDX(T);
%                 plot(Ctime)
%                 title(['subject' num2str(s)]) 
%                 pause
%                 
%             end

            
s=53;
time_interval=2:210;
            
load Total_DELTA_External_SC_FC_included.mat
BOLD = Total_DELTA_External_SC_FC_included(:, [8]);
signal=(BOLD{s,:});
clear Total_DELTA_External_SC_FC_included BOLD

[N, Tmax]=size(signal);


% FILTER SETTINGS
TR=2;
f_low = .02;                    % lowpass frequency of filter
f_high = 0.08;                    %highpass


for seed=1:N
    ts=demean(detrend(signal(seed,:)));
    signal_act(seed,:) = bandpass(ts, [f_low f_high], 1/TR, 0);
    Phase_BOLD(seed,:) = angle(hilbert(signal_act(seed,:)));
end

signal_act=signal_act(:,time_interval);
load LEiDA_data FCD_iFC

FCD_iFCneut=FCD_iFC(:, [1]);

a=time_interval(1):time_interval(end);

figure
colormap(jet)

% Plot the BOLD signals in the top plot
subplot(4,10,[1:10])
plot(a*TR,Phase_BOLD(:,a)','*-')
xlim([(time_interval(1))*TR (time_interval(end))*TR])
xlabel('Time (seconds)')
ylabel('BOLD signal')
set(gca,'xaxisLocation','top')

u=1;
     Order=[1:2:90 90:-2:1];

for t=time_interval(60:2:end)  % In this way it will jump every .. points
    
    %Calculate the Instantaneous FC (BOLD Phase Coherence)
    iFC=zeros(N);
    
    for n=1:N
        for p=1:N
            iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));  % Actually the adif function is useless
        end
    end

    % Get the leading eigenvector
    [V1,~]=eigs(iFC,1);
    Leading_Eig(t,:)=V1;
    iFC_values(t,:)=iFC(triu(ones(N),1)>0);
    
    subplot(4,10,u+10)    % To start in the 3rd line
    imagesc(iFC(Order,Order),[-1 1])
    title(['t=' num2str(t*TR) 's'] )
    
    V1_mat=V1*V1';

    subplot(4,10,u+20)     % To start in the second line
    imagesc(V1_mat(Order,Order))
    title(['t=' num2str(t*TR) 's'])
    
   subplot(4,10,u+30)          % to start in the first line
    plot_nodes_in_cortex(V1);
    title(['t=' num2str(t*TR) 's'])

    
    u=u+1;
   end



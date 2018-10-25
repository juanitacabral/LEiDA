%for figure 1D

load Total_DELTA_External_SC_FC_included.mat
tc_aal = Total_DELTA_External_SC_FC_included(:, [8 10]);
clear Total_DELTA_External_SC_FC_included

load LEiDA_k_results Kmeans_results

[N_areas, Tmax]=size(tc_aal{1,1});

% Select here subject, task
s=20;
task=1;

% FILTER SETTINGS
TR=2;                 
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                    %highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn); 
      
signal = tc_aal{s,task};
        
% Get the BOLD signal using the Hilbert transform
for seed=1:N_areas
    ts=demean(detrend(signal(seed,:)));
    signal_filt(seed,:) =filtfilt(bfilt,afilt,ts);
end

%select the clustering solution
k=4;
[~, ind_sort]=sort(Kmeans_results{k}.SUMD,'descend');


%get Ctime

for cluster=k
        for task=1
            for subject=s
                   
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
            end
        end
end


% Figure BOLD with shaded State Activations

figure

subplot(2,1,1)
hold on
ymax=150;

% Blue Blocks
NetBlue=Ctime==ind_sort(4); %1 or 0
NetON=find(diff(NetBlue)==1); %when it goes ON, value is 1
NetOF=find(diff(NetBlue)==-1); %when it goes OF, value is 0 thus -1

if length(NetON)< length(NetOF) %it was giving an error when lenght of NET ON was greater then NET OFF
    chooselength= length(NetON);
else
    chooselength= length(NetOF);
end

for t=1:chooselength  
    ON=NetON(t); 
    if NetOF(t)>ON  
        OFF=NetOF(t); 
    else
        OFF=NetOF(t+1); 
    end
    x = [ON OFF OFF ON];
    y = [-ymax -ymax ymax ymax];
    p=patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',[0 0 1],'FaceAlpha',0.2);
end

% Red Blocks
NetRed=Ctime==ind_sort(2);
NetON=find(diff(NetRed)==1);
NetOF=find(diff(NetRed)==-1);

if length(NetON)< length(NetOF)
    chooselength= length(NetON);
else
    chooselength= length(NetOF);
end

for t=1:chooselength
    ON=NetON(t);
    if NetOF(t)>ON
        OFF=NetOF(t);
    else
        OFF=NetOF(t+1);
    end
    x = [ON OFF OFF ON];
    y = [-ymax -ymax ymax ymax];
    p=patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',[1 0 0],'FaceAlpha',0.2); 
end


% Green Blocks
NetGreen=Ctime==ind_sort(3);
NetON=find(diff(NetGreen)==1);
NetOF=find(diff(NetGreen)==-1);

if length(NetON)< length(NetOF)
   chooselength= length(NetON);
else
   chooselength= length(NetOF);
end

for t=1:chooselength
    ON=NetON(t);
    if NetOF(t)>ON
        OFF=NetOF(t);
    else
        OFF=NetOF(t+1);
    end
    x = [ON OFF OFF ON];
    y = [-ymax -ymax ymax ymax];
    p=patch(x,y,'r');
    set(p,'LineStyle','none','FaceColor',[0 1 0],'FaceAlpha',0.2);
end

% Plot the BOLD signals in the end
plot(signal_filt')
hold on
%stairs(signal_filt','k')
ylabel('BOLD')
%xlabel('Time (TR)')
ylim([-ymax ymax])
xlim([0 210])
axis off
box off



subplot(2,1,2)
stairs(1+0.8*(Ctime==ind_sort(4)),'b','LineWidth',2)
hold on
stairs(2+0.8*(Ctime==ind_sort(2)),'r','LineWidth',2)
stairs(3+0.8*(Ctime==ind_sort(3)),'g','LineWidth',2)
set(gca,'YTick',1:3)
set(gca,'YTickLabel',['State 1';'State 2';'State 3'])
ylim([0.5 4])
xlim([0 210])
axis off
box off


% %get matrix of state ..
% 
% Order=[1:2:89 90:-2:2];
% 
% load AllVs Vs
% 
% V=Vs(ind_sort(8),:);
% 
% subplot(1,1,1)
% colormap(jet)
%                 FC=V'*V;
%                 imagesc(FC(Order,Order))
%                 axis square
%                 box off
%                 axis off

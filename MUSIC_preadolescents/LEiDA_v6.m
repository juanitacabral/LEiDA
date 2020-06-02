function LEiDA_v6(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
%
% This function processes, clusters and analyses BOLD data using LEiDA.
%
% Here 17 Children aged 10-11 years listening to music inside the scanner
%
% NOTE: Step 4 can be run independently once data is saved by calling
%       LEiDA('LEiDA_results.mat')
%
% 1 - Read the BOLD data from the folders and computes the BOLD phases
%   - Calculate the instantaneous BOLD synchronization matrix
%   - Compute the Leading Eigenvector at each frame from all fMRI scans
% 2 - Cluster the Leading Eigenvectors
% 3 - Compute the probability and lifetimes each cluster in each session

%   - Calculate signigifance between tasks
%   - Saves the Eigenvectors, Clusters and statistics into LEiDA_results.mat
%
% 4 - Plots FC states and errorbars for each clustering solution
%   - Adds an asterisk when results are significantly different between
%   tasks
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Joana Cabral Feb 2018
% joana.cabral@psych.ox.ac.uk
%
% First use in
% Cabral, et al. 2017 Scientific reports 7, no. 1 (2017): 5135.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(varargin) % If no input is given, compute LEiDA
    
    %% 1 - Compute the Leading Eigenvectors from the BOLD datasets
    
    disp('Processing the eigenvectors from BOLD data')
    % Load here the BOLD data (which may be in different formats)
    % Here the BOLD time courses in AAL parcellation are organized as cells,
    % where tc_aal{1,1} corresponds to the BOLD data from subject 1 in
    % condition 1 and contains a matrix with lines=N_areas and columns=Tmax.
    
    Directory='/Users/au548116/Desktop/MusicTraining/Aarhus/ANALYSIS/MusicKids/MusicKids/ANALYSIS_240518/';
    names=dir([Directory '*.txt']);
    n_Subjects=size(names,1);
    n_Task=3; % REST, LEARN, Non-LEARN 
    BOLD=load([Directory '/' names(2).name]);
    [N_areas, Tmax]=size(BOLD);
    
    % Parameters of the data
    TR=1.3;  % Repetition Time (seconds)
    
    % Preallocate variables to save FC patterns and associated information
    Leading_Eig=zeros((Tmax-2)*n_Subjects,N_areas); % All leading eigenvectors
    Time_all=zeros(1, n_Subjects*(Tmax-2)); % vector with subject nr and task at each t
    t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)
    
    % Bandpass filter settings
    fnq=1/(2*TR);                 % Nyquist frequency
    flp = .02;                    % lowpass frequency of filter (Hz)
    fhi = 0.1;                    % highest BOLD frequency conssidered (Hz)
    % In rest studies this is always max 0.1Hz
    Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k=2;                          % 2nd order butterworth filter
    [bfilt,afilt]=butter(k,Wn);   % construct the filter
    clear fnq flp fhi Wn k
    
    for s=1:length(names)
        
        % Get the BOLD signals from subject s
        BOLD = load([Directory '/' names(s).name]);
        % [Tmax]=size(BOLD,2); Get Tmax here, if it changes between scans
        Phase_BOLD=zeros(N_areas,Tmax);
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
            BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
            signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
        
        for t=2:Tmax -1
            
            %Calculate the Instantaneous FC (BOLD Phase Synchrony)
            iFC=zeros(N_areas);
            for n=1:N_areas
                for p=1:N_areas
                    iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
            
            % Get the leading eigenvector
            [V1,~]=eigs(iFC,1);
            % Make sure the largest component is negative
%             if mean(V1)>0
%                 V1=-V1;
%             end           
%                         
            if mean(V1>0)>.5
                V1=-V1;
            elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
                V1=-V1;
            end
            
            % Save V1 from all frames in all fMRI sessions in Leading eig
            t_all=t_all+1; % Update time
            Leading_Eig(t_all,:)=V1;
            Time_all(:,t_all)=s; % Information that at t_all, V1 corresponds to subject s in a given task
        end
    end
    clear BOLD tc_aal signal_filt iFC V1 Phase_BOLD
    
    %% 2 - Cluster the Leading Eigenvectors
    
    disp('Clustering the eigenvectors into')
    % Leading_Eig is a matrix containing all the eigenvectors:
    % Collumns: N_areas are brain areas (variables)
    % Rows: Tmax*n_Subjects are all time points (independent observations)
    
    % Set maximum/minimum number of clusters
    % There is no fixed number of states the brain can display
    % Extend the range depending on the hypothesis of each work
    maxk=15;
    mink=3;
    rangeK=mink:maxk;
    
    % Set the parameters for Kmeans clustering
    Kmeans_results=cell(size(rangeK));
    
    for k=1:length(rangeK)
        disp(['- ' num2str(rangeK(k)) ' clusters'])
        [IDX, C, SUMD, D]=kmeans(Leading_Eig,rangeK(k),'Replicates',500,'MaxIter',5000,'Display','final','Options',statset('UseParallel',0));
        [~, ind_sort]=sort(hist(IDX,1:rangeK(k)),'descend');
        [~,idx_sort]=sort(ind_sort,'ascend');
        Kmeans_results{k}.IDX=idx_sort(IDX);   % Cluster time course - numeric collumn vectors
        Kmeans_results{k}.C=C(ind_sort,:);       % Cluster centroids (FC patterns)
        Kmeans_results{k}.SUMD=SUMD(ind_sort); % Within-cluster sums of point-to-centroid distances
        Kmeans_results{k}.D=D(ind_sort);       % Distance from each point to every centroid   end
    end
    
    save LEiDA_results_v6.mat Leading_Eig Time_all Kmeans_results Tmax n_Subjects n_Task mink maxk 
    
    %% 3 - Analyse the Clustering results

    load VecLvsNL vectorlearnednotlearned
    
    TR_start=1;
    TR_end=3;
    
    T_paradigm=zeros(1,(Tmax-2));
    T_paradigm(64+TR_start:240+TR_end)=1; % First Music Piece
    T_paradigm(264+TR_start:440+TR_end)=2; % Second Music Piece   
    T_paradigm(18+TR_start:40+TR_end)=-1; % Set as -1 all the volumes you want to discard
    
    % For every fMRI scan calculate probability and lifetimes of each pattern c.
    P=zeros(n_Task,n_Subjects,maxk-mink+1,maxk);
    %LT=zeros(n_Task,n_Subjects,maxk-mink+1,maxk);
    
    for k=1:length(rangeK)
        for s=1:n_Subjects
            
            T=(Time_all==s);
            Ctime=Kmeans_results{k}.IDX(T); % 
            
            for task=1:n_Task   % 0, SILENCE, 1 - LEARNED, 2 NON-LEARNED
                % Select the time points representing this subject and task  
                if task==1              
                    Ctime_task=Ctime(T_paradigm==0);
                elseif task>1
                    Piece=find(vectorlearnednotlearned(:,s)==(task-1));
                    Ctime_task=Ctime(T_paradigm==Piece);
                end
                
                for c=1:rangeK(k)
                    % Probability
                    P(task,s,k,c)=mean(Ctime_task==c);
                end
            end
        end
    end
    
    P_pval_RESTvsMUSIC=zeros(maxk-mink+1,maxk);
    P_pval_LEARNvsNONL=zeros(maxk-mink+1,maxk);
    
    disp('Test significance between Silence and Music Learn and non-learn')   
    for k=1:length(rangeK)

        disp(['Now running for ' num2str(rangeK(k)) ' clusters'])
        for c=1:rangeK(k)
            % Compare Probabilities
            a=squeeze(P(1,:,k,c));  % Vector containing Prob of c in Baseline
            b=squeeze(P(2,:,k,c));  % Vector containing Prob of c in Learned
            d=squeeze(P(3,:,k,c));  % Vector containing Prob of c in Non-learned 
            e=[b d];
            stats=permutation_htest2_np([a,e],[ones(1,numel(a)) 2*ones(1,numel(e))],10000,0.05,'ttest');
            P_pval_RESTvsMUSIC(k,c)=min(stats.pvals);     
            stats=permutation_htest2_np([b,d],[ones(1,numel(b)) 2*ones(1,numel(d))],10000,0.05,'ttest');
            P_pval_LEARNvsNONL(k,c)=min(stats.pvals);   
        end
    end 
    disp('%%%%% LEiDA SUCCESSFULLY COMPLETED %%%%%%%')
    disp('Saving LEiDA results')
    save LEiDA_results_v6.mat Leading_Eig Time_all Kmeans_results P P_pval_LEARNvsNONL P_pval_RESTvsMUSIC rangeK
    
else
    load(varargin{1})
end

%% 4 - Plot FC patterns and stastistics between groups

load AAL_labels label90
disp(' ')
disp('%%% PLOTS %%%%')
disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])

K = input('Number of clusters: ');

Best_Clusters=Kmeans_results{rangeK==K};
k=find(rangeK==K);

% Get the K patterns
V=Best_Clusters.C;
[~, N]=size(Best_Clusters.C);
Order=[1:2:N N:-2:2];

cmap=[0 0 1;  1 0 0 ; .7 .7 .7 ; 1 0.5 0; 0 1 1 ; 1 0 1; 1 1 0; 0 1 0];
    
P_pval=P_pval_RESTvsMUSIC;

% % % load BMRQ_scores.mat Scores
% % % 
figure
colormap(jet) 
% % Pannel A - Plot the FC patterns over the cortex 
% % Pannel B - Plot the FC patterns in matrix format
% % Pannel C - Plot the probability of each state in each condition
% % Pannel D - Plot the lifetimes of each state in each condition
%    
for c=1:K
    subplot(3,K,c)
    % This needs function plot_nodes_in_cortex.m and aal_cog.m
    plot_nodes_in_cortex(V(c,:))
    title({['State #' num2str(c)]})
    
    subplot(3,K,c+K)
    FC_V=V(c,:)'*V(c,:);  
    li=max(abs(FC_V(:)));
    imagesc(FC_V(Order,Order),[-li li])   
    axis square
    title('FC pattern') 
    ylabel('Brain area #')
    xlabel('Brain area #')   
    
    subplot(3,K,c+K*2)  
            Rest=squeeze(P(1,:,k,c));
            Music=[squeeze(P(2,:,k,c)) squeeze(P(3,:,k,c))];
            bar([mean(Rest) mean(Music)],'EdgeColor','none','FaceColor',[.5 .5 .5])
            hold on
            % Error bar containing the standard error of the mean
            errorbar([mean(Rest) mean(Music)],[std(Rest)/sqrt(numel(Rest)) std(Music)/sqrt(numel(Music))],'LineStyle','none','Color','k')
            set(gca,'XTick',[1 2])            
            set(gca,'XTickLabel',{'Silence','Music'})
            set(gca,'XTickLabelRotation',45)
            if P_pval(k,c)<0.05/K
                plot(1.5,max([mean(Rest) mean(Music)])+.05,'*k')
            end      
            if c==1
                ylabel('Probability')
                ylim([0 0.45])
            else
                ylim([0 0.25])
            end
            box off
            xlim([0 3])

end


N_areas=90;
Order=[1:2:N_areas N_areas:-2:2];
label90=label90(Order,:);

figure
for c=1:K
    subplot(1,K,c)
    Vc=V(c,Order);
    hold on
    barh(Vc.*(Vc<0),'FaceColor',[0.2  .2  1],'EdgeColor','none','Barwidth',.5)
    barh(Vc.*(Vc>=0),'FaceColor',[1 .2 .2],'EdgeColor','none','Barwidth',.5)
    ylim([0 91])
    xlim([-.15 .15])
    set(gca,'YTick',1:N_areas,'Fontsize',8)
    if c==1
        set(gca,'YTickLabel',label90(end:-1:1,:),'Fontsize',6)
    else
        set(gca,'YTickLabel',[])
    end
    ylim([0 91])
    xlim([-.15 .15])
    set(gca,'Ydir','reverse')
    title(['State ' num2str(c)])
    grid on
end


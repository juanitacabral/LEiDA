function LEiDA(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)
%
% This function processes, clusters and analyses BOLD data using LEiDA.
% Here the example_BOLD is a dataset containing rest and task conditions
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

% Joana Cabral Oct 2017
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
    
    load example_BOLD.mat tc_aal
    [n_Subjects, n_Task]=size(tc_aal);
    [N_areas, Tmax]=size(tc_aal{1,1});
    
    % Parameters of the data
    TR=2;  % Repetition Time (seconds)
    
    % Preallocate variables to save FC patterns and associated information
    Leading_Eig=zeros(Tmax*n_Subjects,N_areas); % All leading eigenvectors
    Time_all=zeros(2, n_Subjects*Tmax); % vector with subject nr and task at each t
    t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)
    
    % Bandpass filter settings
    fnq=1/(2*TR);                 % Nyquist frequency
    flp = .02;                    % lowpass frequency of filter (Hz)
    fhi = 0.1;                    % highpass
    Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
    k=2;                          % 2nd order butterworth filter
    [bfilt,afilt]=butter(k,Wn);   % construct the filter
    clear fnq flp fhi Wn k
    
    for s=1:n_Subjects
        for task=1:n_Task
            
            % Get the BOLD signals from this subject in this task
            BOLD = tc_aal{s,task};
            % [Tmax]=size(BOLD,2); Get Tmax here, if it changes between scans
            Phase_BOLD=zeros(N_areas,Tmax);
            
            % Get the BOLD phase using the Hilbert transform
            for seed=1:N_areas
                BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
                signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
                Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
            end
            
            for t=1:Tmax
                
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
                % This step is important because the same eigenvector can 
                % be returned either as V or its symmetric -V and we need
                % to make sure it is always the same (so we choose always
                % the most negative one)
                if mean(V1>0)>.5
                    V1=-V1;
                elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
                    V1=-V1;
                end
                
                % Save V1 from all frames in all fMRI sessions in Leading eig
                t_all=t_all+1; % Update time
                Leading_Eig(t_all,:)=V1;
                Time_all(:,t_all)=[s task]; % Information that at t_all, V1 corresponds to subject s in a given task
            end
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
    maxk=12;
    mink=3;
    rangeK=mink:maxk;
    
    % Set the parameters for Kmeans clustering
    Kmeans_results=cell(size(rangeK));
    
    for k=1:length(rangeK)      
        disp(['- ' num2str(rangeK(k)) ' clusters'])
        [IDX, C, SUMD, D]=kmeans(Leading_Eig,rangeK(k),'Replicates',10,'MaxIter',1000,'Display','off','Options',statset('UseParallel',1));
        Kmeans_results{k}.IDX=IDX;   % Cluster time course - numeric collumn vectos
        Kmeans_results{k}.C=C;       % Cluster centroids (FC patterns)
        Kmeans_results{k}.SUMD=SUMD; % Within-cluster sums of point-to-centroid distances
        Kmeans_results{k}.D=D;       % Distance from each point to every centroid
    end
    
    save LEiDA_results.mat Leading_Eig Time_all Kmeans_results
    
    %% 3 - Analyse the Clustering results
    
    % For every fMRI scan calculate probability and lifetimes of each pattern c.
    P=zeros(n_Task,n_Subjects,maxk-mink+1,maxk);
    LT=zeros(n_Task,n_Subjects,maxk-mink+1,maxk);
    
    for k=1:length(rangeK)
        for task=1:n_Task   % 1, Baselineline, 2, baby face task
            for s=1:n_Subjects
                
                % Select the time points representing this subject and task               
                T=((Time_all(1,:)==s)+(Time_all(2,:)==task))>1;
                Ctime=Kmeans_results{k}.IDX(T);
                
                for c=1:rangeK(k)
                    % Probability
                    P(task,s,k,c)=mean(Ctime==c);
                    
                    % Mean Lifetime
                    Ctime_bin=Ctime==c;
                    
                    % Detect switches in and out of this state
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
                    LT(task,s,k,c)=mean(C_Durations)*TR;
                end                
            end
        end
    end   
    
    P_pval=zeros(maxk-mink+1,maxk);
    LT_pval=zeros(maxk-mink+1,maxk);
    
    disp('Test significance between Rest and Task')   
    for k=1:length(rangeK)

        disp(['Now running for ' num2str(k) ' clusters'])
        for c=1:rangeK(k)
            % Compare Probabilities
            a=squeeze(P(1,:,k,c));  % Vector containing Prob of c in Baseline
            b=squeeze(P(2,:,k,c));  % Vector containing Prob of c in Task
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
            P_pval(k,c)=min(stats.pvals);
            
            % Compare Lifetimes
            a=squeeze(LT(1,:,k,c));  % Vector containing Lifetimes of c in Baseline
            b=squeeze(LT(2,:,k,c));  % Vector containing Lifetimes of c in Task
            stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest');
            LT_pval(k,c)=min(stats.pvals);            
        end
    end 
    disp('%%%%% LEiDA SUCCESSFULLY COMPLETED %%%%%%%')
    disp('Saving LEiDA results')
    save LEiDA_results.mat Leading_Eig Time_all Kmeans_results P LT P_pval LT_pval rangeK
else
    load(varargin{1})
end

%% 4 - Plot FC patterns and stastistics between groups

disp(' ')
disp('%%% PLOTS %%%%')
disp(['Choose number of clusters between ' num2str(rangeK(1)) ' and ' num2str(rangeK(end)) ])
Pmin_pval=min(P_pval(P_pval>0));
LTmin_pval=min(LT_pval(LT_pval>0));
if Pmin_pval<LTmin_pval
   [k,~]=ind2sub([length(rangeK),max(rangeK)],find(P_pval==Pmin_pval));
else
   [k,~]=ind2sub([length(rangeK),max(rangeK)],find(LT_pval==LTmin_pval));
end
disp(['Note: The most significant difference is detected with K=' num2str(rangeK(k)) ' (p=' num2str(min(Pmin_pval,LTmin_pval)) ')'])  

% To correct for multiple comparisons, you can divide p by the number of
% clusters considered

K = input('Number of clusters: ');
Best_Clusters=Kmeans_results{rangeK==K};
k=find(rangeK==K);

% Clusters are sorted according to their probability of occurrence
ProbC=zeros(1,K);
for c=1:K
    ProbC(c)=mean(Best_Clusters.IDX==c);
end
[~, ind_sort]=sort(ProbC,'descend'); 

% Get the K patterns
V=Best_Clusters.C(ind_sort,:);
[~, N]=size(Best_Clusters.C);
Order=[1:2:N N:-2:2];

figure
colormap(jet) 
% Pannel A - Plot the FC patterns over the cortex 
% Pannel B - Plot the FC patterns in matrix format
% Pannel C - Plot the probability of each state in each condition
% Pannel D - Plot the lifetimes of each state in each condition
   
for c=1:K
    subplot(4,K,c)
    % This needs function plot_nodes_in_cortex.m and aal_cog.m
    plot_nodes_in_cortex(V(c,:))
    title({['State #' num2str(c)]})
    subplot(4,K,K+c)
    FC_V=V(c,:)'*V(c,:);  
    li=max(abs(FC_V(:)));
    imagesc(FC_V(Order,Order),[-li li])   
    axis square
    title('FC pattern') 
    ylabel('Brain area #')
    xlabel('Brain area #')   
    
    subplot(4,K,2*K+c)  
            Rest=squeeze(P(1,:,k,ind_sort(c)));
            Task=squeeze(P(2,:,k,ind_sort(c)));
            bar([mean(Rest) mean(Task)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            % Error bar containing the standard error of the mean
            errorbar([mean(Rest) mean(Task)],[std(Rest)/sqrt(numel(Rest)) std(Task)/sqrt(numel(Task))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'Rest', 'Task'})
            if P_pval(k,ind_sort(c))<0.05
                plot(1.5,max([mean(Rest) mean(Task)])+.01,'*k')
            end             
            if c==1
                ylabel('Probability')
            end
            box off
            
     subplot(4,K,3*K+c)  
            Rest=squeeze(LT(1,:,k,ind_sort(c)));
            Task=squeeze(LT(2,:,k,ind_sort(c)));
            bar([mean(Rest) mean(Task)],'EdgeColor','w','FaceColor',[.5 .5 .5])
            hold on
            errorbar([mean(Rest) mean(Task)],[std(Rest)/sqrt(numel(Rest)) std(Task)/sqrt(numel(Task))],'LineStyle','none','Color','k')
            set(gca,'XTickLabel',{'Rest', 'Task'})
            if LT_pval(k,ind_sort(c))<0.05
                plot(1.5,max([mean(Rest) mean(Task)])+.01,'*k')
            end             
            if c==1
                ylabel('Lifetime (seconds)')
            end
            box off           
end
function LEiDA_data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEADING EIGENVECTOR DYNAMICS ANALYSIS
%
% This function processes recurrent MDD data for LEiDA
%
% - Reads the BOLD data from the folders
% - Computes the BOLD phases
% - Calculates the Order Parameter, Mestastability and Synchrony
% - Calculates the instantaneous BOLD synchronization matrix
% - Calculates the instantaneous Leading Eigenvector
% - Calculates de FCD matrices
%
%	The data include 86 participants
%       -51 patients; 35 healthy controls
%       -Patients additionally include relapsed or no-relapse
%   Scans include 210 TRs each (TR = )
%
%   Order of variables in data structure:
%   -1 subjnr
%   -2 BaselineHDRS
%   -3 PT_HC **
%   -4 HC_NoRelapse_Relapse
%   -5 TimeToEvent
%   -6 Prev.Episodes before Bsl
%   -7 SC (LR-oriented)
%   -8 FC_TC_NoMood (LR-oriented) **
%   -9 FC-z_NoMood (LR-oriented)
%   -10 FC_TC_Mood (LR-oriented) **
%   -11 FC_z_Mood (LR-oriented)
%
% Saves the Leading_Eigenvectors and FCD matrices into LEiDA_data.mat
%
% Written by Joana Cabral, May 2016 (joana.cabral@psych.ox.ac.uk)
% Modified by Kristina Rapuano, October 2016
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
load Total_DELTA_External_SC_FC_included.mat
tc_aal = Total_DELTA_External_SC_FC_included(:, [8 10]);

[n_Subjects, n_Task]=size(tc_aal);
[N_areas, Tmax]=size(tc_aal{1,1});

% Create empty variable to save patterns 
Leading_Eig=cell(n_Subjects, n_Task);
FCD_eig=cell(n_Subjects, n_Task);
FCD_iFC=cell(n_Subjects, n_Task);
iFC_values=zeros(Tmax,(N_areas*(N_areas-1)/2));

% FILTER SETTINGS
TR=2;                 
fnq=1/(2*TR);                 % Nyquist frequency
flp = .02;                    % lowpass frequency of filter
fhi = 0.1;                    %highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
f_low = 0.04;                 % narrow bandpass filter
f_high = 0.08;

for s=1:n_Subjects
    
    for task=1:n_Task
        
        signal = tc_aal{s,task};
        
        % Get the BOLD phase using the Hilbert transform
        for seed=1:N_areas
            signal_act(seed,:) = bandpass(signal(seed,:), [f_low f_high], 1/TR, 0);
            ts=demean(detrend(signal(seed,:)));
            signal_filt =filtfilt(bfilt,afilt,ts);
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
        end
        
        %Calculate the Kuramoto Order Parameter
        OP=abs(sum(exp(1i*Phase_BOLD))/N_areas);
        Metasta(s,task)=std(OP);
        Synchro(s,task)=mean(OP);
        
        for t=1:Tmax
            
            %Calculate the Instantaneous FC (BOLD Phase Coherence)
            iFC=zeros(N_areas);
            
            for n=1:N_areas
                for p=1:N_areas
                    iFC(n,p)=cos(adif(Phase_BOLD(n,t),Phase_BOLD(p,t)));
                end
            end
            
            % Get the leading eigenvector
            [Leading_Eig{s,task}(t,:),~]=eigs(iFC,1);
            iFC_values(t,:)=iFC(triu(ones(N_areas),1)>0);
        end
        
        % Calculate the FCD
        for t=1:Tmax
            eig1=squeeze(Leading_Eig{s,task}(t,:));
            iFC1=iFC_values(t,:);
            for t2=1:Tmax
                eig2=squeeze(Leading_Eig{s,task}(t2,:));
                iFC2=iFC_values(t2,:);
                % Cosine similarity between vectors at t1 and t2
                FCD_eig{s,task}(t,t2)=dot(eig1,eig2)/norm(eig1)/norm(eig2);
                FCD_iFC{s,task}(t,t2)=dot(iFC1,iFC2)/norm(iFC1)/norm(iFC2);
            end
        end
    end
end

save('LEiDA_data','FCD_eig','FCD_iFC','Metasta','Synchro','Leading_Eig')






LEADING EIGENVECTOR DYNAMICS ANALYSIS (LEiDA)

The LEiDA consists of computing instantaneous BOLD phase coherence matrices 
and clustering the corresponding leading eigenvectors into a set of 
patterns that can be visualized on the cortical surface or in matrix format.
  
This repository includes the codes and data to replicate the analysis
reported in the manuscript:

Fine-grained analysis of functional connectivity dynamics links cognitive 
performance in healthy aging to spontaneous switching between brain states 
J Cabral, D Vidaurre, P Marques, R Magalhães, P Silva Moreira, JM Soares,
G Deco, N Sousa and ML Kringelbach

The first function is LEiDA_data.m that loads the BOLD data from 
Aging_data.mat, computes the  instantaneous phases and saves the leading 
eigenvector at each time point.

On a second step, LEiDA_Cluster.m applies a kmeans clustering algorithm to 
the leading eigenvectors and saves the optimal solution.

Finally, LEiDA_analysis.m loads the results and plots the optimal set of 
vectors on the cortical surface and compares the FC patterns with the 
average BOLD FC in Static_FC.mat 

LEiDA needs the SPM functions spm_vol() and spm_slice_vol(), so SPM needs 
be added to the Matlab pathway before running

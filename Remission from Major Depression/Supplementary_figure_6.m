% Patients vs Controls, plot most significant networks

%For probability and lifetime k=10, c=8

load LEiDA_k_results Kmeans_results
load Between_results_sept_2018 P pvalue_PvsC LT LT_pval_PvsC

load AAL_labels.mat label90
Order=[1:2:89 90:-2:2];


            %%suplemental figure 6 neutral vs sad mood rrMDD plot the matrices matrix
            
            for k=10
        for c=1:10
  
            
                
                figure('Name',['K=' num2str(k) ' c=' num2str(c) ])
                colormap(jet)
                V=Kmeans_results{k}.C(c,:);
             
                subplot(2,2,1)               
                FC=V'*V;
                imagesc(FC(Order,Order))
                box off
            
        end
            end
            
            %plot the networks from k=1 to 10 for every c at a time 
            for k=10
        for c=1
         
                plot_nodes_in_cortex_new(V);      
        end
            end
            
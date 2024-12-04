# -*- coding: utf-8 -*-
"""

This script generates a lens function and creates a Mapper graph for each specific parameter combination on the METABRIC discovery dataset. 
These graphs are searched for hotspots, and survival analysis is performed between the hotspot group and the neighbourhood. 
If the logrank p-value < 0.01, the search is considered successful and ends. 

For all hotspots in a single lens function that can be generated by multiple interval and overlap combinations, 
these hotspots are ranked by p-value and the hotspot with the lowest p-value is accepted as the final parameter selection.

The accepted graph parameters from the hotspot search are saved in 
'outputs/hotspot_search/discovery/"
 directory, but the succesful parameters identified for the paper are saved in 
 'outputs/hotspot_search/discovery_final_results/" to prevent overwriting. 

"""


import hot_mapper as hm
import numpy as np
import pandas as pd
import hdbscan
from lifelines.statistics import logrank_test




#### SET UP EXPERIMENT
project_directory = "..."
output_path = f"{project_directory}/output/hotspot_search/discovery"
dataset_name = "metabric"


#### READ IN FILES
X = pd.read_csv(f"{project_directory}/output/processed_data/{dataset_name}_dct.csv", index_col = 0)
survival_df = pd.read_csv(f"{project_directory}/output/processed_data/{dataset_name}_survival.csv", index_col = 0)


#### USE SURVIVAL FOR ATTRIBUTE FUNCTION
# We want to search for patients experiencing survival event before 10 years
survival_df.index = list(survival_df["patient_id"])
survival_df = survival_df.drop("patient_id", axis = 1)
rfs = np.array((survival_df['time']<=120) & (survival_df['event'] == 1)).astype(int)





#### SET UP PARAMETERS
#initialise the search class from the hotmapper module
search = hm.automated_parameter_search.Search(np.array(X))

#select the parameter options for the search 
parameters = {"predefined_lens" : None, # This parameter is only used for the validation set when we have found a lens
               "non_zero_lens_features" : int(X.shape[1]/2), #50\% of non-zero features in the lens function
              "interval_list" : range(10,32,2) , #no. of interval options
              "overlap_list" : np.linspace(0.1,0.45,8), #percentage of overlap options
              "clustering_algorithm" : hdbscan.HDBSCAN(), #keep clustering algorithm consistent 
              "attribute_function" : rfs, #patients who have a survival event before 10 years
              "epsilon" : 0.1, #minimum difference in attribute between hotspot and neighbourhood
              'min_samples' : 30, #minimum sample size for hotspot
              'extreme': "higher"} #hotspots with higher occurence of event before 10 years 


#### RUN SEARCH
#specify the number of times to run a search for a hotspot with signficant survival difference 
runs = 100 
signficance = False
count = 0

while count < runs: 
    if signficance == False:
        print("\nSearching lens space...")
        print(F"{count} / {runs}")
        
        #### HOTSPOT SEARCH
        #mapper graphs are built for each lens and tested for the presence of a hotspot
        search.build_graphs(parameters)
        
        #the output is any sucessful parameters creating a graph containing hotspots
        p_success = search.parameters

        #the attributes to build the lens function are in parameter_lens
        weights = search.parameter_lens["weights"]
        feature_list = search.parameter_lens["feature_list"]
        
        
        #### SURVIVAL ANALYSIS 
        survival_results = []
        hotspot_id = []
        columns = ["interval", "overlap", "nodes", "size", "logrank"]
        
        # multiple graphs are generated from the different successful parameter options 
        # the hotspots in each graph are tested individually for significant survival 
        for ps, collection in search.parameter_samples.items():
            for i, hotspot_samples in enumerate(collection):
                hotspot_results = []
                
                #seperate the er+ cohort into hotspot and global neighbourhood
                y = pd.DataFrame([0] * len(X.index), columns = ["hotspot"], index = X.index)
                X_H = X.iloc[hotspot_samples]
                y.loc[X_H.index,"hotspot"] = 1
                
                #apply log rank test to each hotspot division 
                hot = (y["hotspot"] == 1)
                T = survival_df["time"]
                E = survival_df["event"]
        
                lr = logrank_test(T[hot], T[~hot], event_observed_A=E[hot], event_observed_B=E[~hot], alpha=99)
                pvalue = lr.p_value
                
                
                # if any hotspots have lower p-value than 0.01 then results are saved 
                if pvalue < 0.01:
                    
                    # append results to list to build a dataframe summarising survival analysis 
                    # for each hotspot generated from different parameter combinations across a lens function
                    hotspot_id.append(str(ps[0]) + str(int(ps[1] * 100)) + str(i))
                    hotspot_results = [ps[0], ps[1], p_success[ps][i], sum(y["hotspot"]), pvalue]
                    survival_results.append(hotspot_results)
                    signficance = True

                    
        #if any hotspots have lower p-value than 0.001 then search ends
        hotspot_df = pd.DataFrame(survival_results, columns = columns, index = hotspot_id)
        if hotspot_df.empty:
            count += 1
            print("Hotspot not significant ... continue search") 
            
            
    #### SAVE RESULTS TO FILE
    else:
        # Order results to find the hotspot with strongest survival difference between neighbourhood
        hotspot_df = hotspot_df.sort_values(by = 'logrank')
        top_parameters = list(hotspot_df[["interval","overlap"]].iloc[0])
        hotspot_df.to_csv(f"{output_path}/{dataset_name}_hotspot_dataframe.csv")
    
        # save weights and order of features from randomly generated lens function to allow us to receate it later
        np.savetxt(f"{output_path}/{dataset_name}_parameters.txt", top_parameters)
        np.savetxt(f"{output_path}/{dataset_name}_weights.txt", weights ,delimiter=",")
        np.savetxt(f"{output_path}/{dataset_name}_feature_list.txt", feature_list ,delimiter=",")
        print("Hotspot significantly impacts survival \n search ends.") 

        break


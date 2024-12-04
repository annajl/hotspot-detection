# -*- coding: utf-8 -*-
"""
This script uses the successful lens function identified through the hotspot search of the METABRIC discovery cohort, 
and searches for hotspots in the TCGA validation cohort across interval and overlap parameters. 
Graphs are generated for each parameter combination, these graphs are searched for hotspots, 
and survival analysis is performed between the hotspot group and the neighbourhood. 

The accepted graph parameters from the hotspot search are saved in 
'outputs/hotspot_search/validation/" directory, 
but the succesful parameters identified for the paper are saved in 
'outputs/hotspot_search/validaton_final_results/" 
to prevent overwriting. 


"""

import hot_mapper as hm
import numpy as np
import pandas as pd
import hdbscan
from lifelines.statistics import logrank_test



#### SET UP EXPERIMENT
project_directory = "..."
output_path = f"{project_directory}/output/hotspot_search/validation"
dataset_name = "tcga"
discovery_dataset = "metabric"

#### READ IN FILES
X = pd.read_csv(f"{project_directory}/output/processed_data/{dataset_name}_dct.csv", index_col = 0)
survival_df = pd.read_csv(f"{project_directory}/output/processed_data/{dataset_name}_survival.csv", index_col = 0)


#### USE SURVIVAL FOR ATTRIBUTE FUNCTION
# We want to search for patients experiencing survival event before 10 years
survival_df.index = list(survival_df["patient_id"])
survival_df = survival_df.drop("patient_id", axis = 1)
os = np.array((survival_df['time']<=120) & (survival_df['event'] == 1)).astype(int)



#### BUILD LENS FUNCTION IDENTIFIED ON DISCOVERY DATASET
parameters_file_path = f"{project_directory}/output/hotspot_search/discovery_final_results" 
weights = np.genfromtxt(f"{parameters_file_path}/{discovery_dataset}_weights.txt")
feature_list = np.genfromtxt(f"{parameters_file_path}/{discovery_dataset}_feature_list.txt", dtype = "i4")

#build the predefined lens from the discovery search on the validation data
discovery_lens = hm.random_lens.Lens(np.array(X), nonzero_features = len(weights), weights = weights, feature_list = feature_list) 



#### SET SEARCH PARAMETERS
#initialise the search class from the hotmapper module
search = hm.automated_parameter_search.Search(np.array(X))

#select the parameter options for the search 
parameters = {'predefined_lens': discovery_lens, #use lens identified from discovery search
              "interval_list" : range(10,32,2) , #no. of interval options
              "overlap_list" : np.linspace(0.1,0.45,8), #percentage of overlap options
              "clustering_algorithm" : hdbscan.HDBSCAN(), #keep clustering algorithm consistent 
              "attribute_function" : os, #patients who have a survival event before 10 years
              "epsilon" : 0.1,#minimum difference in attribute between hotspot and neighbourhood
              'min_samples' : 30, #minimum sample size for hotspot
              'extreme': "higher"} #hotspots with higher occurence of event before 10 years 

#


##### RUN SEARCH
print("\nSearching lens space...")

#mapper graphs are built for each lens and tested for the presence of a hotspot
search.build_graphs(parameters)

#the sucessful parameters containing hotspots
p_success = search.parameters


#### PERFORM SURVIVAL ANALYSIS ON HOTSPOTS
survival_results = []
hotspot_id = []
columns = ["interval", "overlap", "nodes", "size", "logrank"]

#multiple graphs are generated from the different successful parameter options 
#the hotspots in each graph are tested for significant survival 
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
        
        
        #if any hotspots have lower p-value than 0.001 then results are saved 
        #append results to list to build dataframe summarising survival analysis for each hotspot 
        hotspot_id.append(str(ps[0]) + str(int(ps[1] * 100)) + str(i))
        hotspot_results = [ps[0], ps[1], p_success[ps][i], sum(y["hotspot"]), pvalue]
        survival_results.append(hotspot_results)


#if any hotspots have lower p-value than 0.001 then search ends
hotspot_df = pd.DataFrame(survival_results, columns = columns, index = hotspot_id)
print(hotspot_df)



hotspot_df = hotspot_df.sort_values(by = 'logrank')
top_parameters = list(hotspot_df[["interval","overlap"]].iloc[1])
hotspot_df.to_csv(f"{output_path}/{dataset_name}_hotspot_dataframe.csv")

#save weights and order of features from randomly generated lens function to allow us to receate it later
np.savetxt(f"{output_path}/{dataset_name}_parameters.txt", top_parameters)
np.savetxt(f"{output_path}/{dataset_name}_weights.txt", weights ,delimiter=",")
np.savetxt(f"{output_path}/{dataset_name}_feature_list.txt", feature_list ,delimiter=",")





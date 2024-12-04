# -*- coding: utf-8 -*-
"""

This script uses the successful parameters identified through the hotspot search of the METABRIC discovery cohort. 
It builds the graph to visualise the output and annotates the nodes that have been identified as hotspots. 
It also appends hotspot labels to the METABRIC sample IDs and their corresponding survival information. 

To recreate the mapper graph from the paper, we are using an input directory of 
'outputs/hotspot_search/discovery_final_results/" 
for the mapper graph parameter settings. Change this to 
'outputs/hotspot_search/discovery/" directory 
if you want to see the output from a new hotspot search on the discovery group. 

"""


import pandas as pd
import numpy as np

import hot_mapper as hm
from hdbscan import HDBSCAN

from itertools import chain



#### SET UP EXPERIMENT
project_directory = "..."
output_path = f"{project_directory}/output/mapper_graphs"
dataset_name = "metabric"


#### READ IN FILES
X = pd.read_csv(f"{project_directory}/output/processed_data/{dataset_name}_dct.csv", index_col = 0)
survival_df = pd.read_csv(f"{project_directory}/output/processed_data/{dataset_name}_survival.csv", index_col = 0)



#### USE SURVIVAL FOR ATTRIBUTE FUNCTION
# We want to search for patients experiencing survival event before 10 years
survival_df.index = list(survival_df["patient_id"])
survival_df = survival_df.drop("patient_id", axis = 1)
rfs = np.array((survival_df['time']<=120) & (survival_df['event'] == 1)).astype(int)





#### USE MAPPER PARAMETERS FROM HOTSPOT SEARCH
## Change as appropriate - If you have conducted a new search, results will be saved in 
## f"{project_directory}/hotspot_search/discovery"
## The results from the paper are found in 
## f"{project_directory}/hotspot_search/discovery_final_results"
parameters_file_path = f"{project_directory}/output/hotspot_search/discovery_final_results" 
parameters = np.genfromtxt(f"{parameters_file_path}/{dataset_name}_parameters.txt")

# Use the interval and overlap parameters generating a hotspot with lowest logrank p-value
intervals = int(parameters[0])
overlap = parameters[1]
print(f"Number of intervals: {intervals}, \nOverlap percentage: {overlap * 100}%")

# Use feature weights and feature selection from successful lens 
# These two variable will be used as input to the Lens class
weights = np.genfromtxt(f"{parameters_file_path}/{dataset_name}_weights.txt")
feature_list = np.genfromtxt(f"{parameters_file_path}/{dataset_name}_feature_list.txt", dtype = "i4")

# Specify number of nonzero features in lens
feature_no = int(X.shape[1]/2)





#### BUILD MAPPER
# Generate the successful lens function
linear_lens = hm.random_lens.Lens(np.array(X), 
                                  nonzero_features = feature_no, 
                                  weights = weights, 
                                  feature_list = feature_list)


# Specify the parameters for the Mapper graph
mapper = hm.mapper.MapperGraph(data = np.array(X), 
                                lens_function = linear_lens["lens"], 
                                intervals = intervals, 
                                overlap = overlap,
                                clustering_algorithm = HDBSCAN())

# Build the network graph                
mapper.build_graph() 

# Visualise the network graph
hm.visualisation.draw_graph(mapper_graph = mapper.graph, #input the constructed mapper graph
                              attribute_function = rfs,  # colouring of nodes
                              samples_in_nodes = mapper.samples_in_nodes, # specify the samples distribution across nodes
                              size = 10,  # size of nodes
                              style = 3, # networkx layout style
                              col_legend_title = "Relapse before\n 10 years", # title of legend 
                              labels = True, # include node labes
                              tick_labels = True, # tick labels on the legend 
                              file_name = f"{output_path}/{dataset_name}_survival_mapper_labelled.png")





#### SEARCH FOR HOTSPOTS 
# Initialise the hotspot search class with the mapper graph
hotspot_search = hm.hotspot.HotspotSearch(mapper_graph = mapper.graph,
                                         attribute_function = rfs, 
                                         samples_in_nodes = mapper.samples_in_nodes)

# Perform a hotspot search on the nodes
hotspot_nodes = hotspot_search.search_graph(attribute_threshold = 0.1, # minimum difference in attribute between hotspot and neighbourhood
                                            min_sample_size = 30, # minimum hotspot sample size
                                            attribute_extreme = "higher",  # selecting hotspots in upper extreme of attribute
                                            plot_dendrogram = True ) #plot a dendogram showing connectivity of nodes

# Only one hotspot found, obtain the node numbers
hotspot_node_list = list(chain.from_iterable(hotspot_nodes))

# Plot using same settings as before
# Include hotspot_node variable
hm.visualisation.draw_graph(mapper_graph = mapper.graph, 
          attribute_function = rfs, 
          samples_in_nodes = mapper.samples_in_nodes,
          size = 10, 
          style = 3, 
          hotspot_nodes = hotspot_node_list, # Highlight nodes in hotspot
          labels = False, # remove node labels for visibility
          tick_labels = True,
          col_legend_title = "Relapse before\n 10 years",
          file_name = f"{output_path}/{dataset_name}_mapper_survival_hotspots_unlabelled.png")


print(f"hotspots found in nodes... {hotspot_node_list}")




#### SAVE THE HOTSPOT LABELS IDENTIFIED FOR SAMPLES
# save the samples assigned to each node 
node_id = mapper.samples_in_nodes
node_id.index = survival_df.index

# append hotspot label to survival information
survival_df["hotspot"] = node_id[hotspot_node_list].max(axis=1)
survival_df.to_csv(f"{project_directory}/output/processed_data/{dataset_name}_hotspot_id_survival.csv")

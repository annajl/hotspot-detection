# -*- coding: utf-8 -*-
"""

This script prepares the gene expression datasets for input into the DSGA method. Missing data is handled for GTEX by removing genes and imputing using KNN for TCGA. 

For each breast cancer cohort, a pair of dataframes is created by matching the gene features in the tumour dataset to those in the GTEX dataset. 

"""


import numpy as np 
import pandas as pd
from sklearn.impute import KNNImputer


def check_for_missing_data(list_of_datasets):
    sum_missing = [data[data.isnull().any(axis=1)].shape[0] for data in list_of_datasets]     
    return sum_missing 


def match_tumour_to_normal_data(list_of_tumour_datasets, normal_dataset):
    """Restrict features to common genes with GTEX normal tissue data for each breast tumour cohort"""
    
    for dataset_name,tumour_dataset in list_of_tumour_datasets.items():
        #transpose data
        tumour_dataset_T = tumour_dataset.T
        normal_dataset_T = normal_dataset.T
        
        #obtain list of genes in both datasets 
        tumour_genes = tumour_dataset_T.columns
        normal_genes = normal_dataset_T.columns
    
        #match the names 
        matched_genes = [i for i in normal_genes if i in tumour_genes]
        
        #subset original datasets
        tumour_matched = tumour_dataset_T[matched_genes]
        normal_matched = normal_dataset_T[matched_genes]
        
        print(F"{dataset_name}: \n tumour shape ({tumour_matched.shape} \n normal shape {normal_matched.shape})")
        
        #save to file 
        tumour_matched.to_csv(F"{output_path}/{dataset_name}_tumour_matched.csv") 
        normal_matched.to_csv(F"{output_path}/{dataset_name}_normal_matched.csv") 

        


#### SET UP EXPERIMENT
project_directory = "..."
output_path = f"{project_directory}/output/processed_data"


#### READ IN FILES
metabric = pd.read_csv(filepath_or_buffer = F"{project_directory}/data/metabric_protein_coding_genes_expression_zscore.csv", sep = ",", index_col = 0)
tcga = pd.read_csv(filepath_or_buffer = F"{project_directory}/data/tcga_protein_coding_genes_expression_zscore.csv", sep = ",", index_col = 0)
gtex = pd.read_csv(filepath_or_buffer = F"{project_directory}/data/gtex_gene_expression_zscore.csv", sep = ",", index_col = 0)

# format GTEX dataset correctly
gtex = gtex.set_index('Gene') 



#### HANDLE MISSING DATA
#check for missing data in the datasets         
missing_data_no = check_for_missing_data([metabric,tcga,gtex])
print(F"number of missing genes: {missing_data_no}")

#There is missing data for 434 genes in the TCGA dataset, these are imputed using KNN
imputer = KNNImputer(n_neighbors=10)
tcga_imputed = pd.DataFrame(imputer.fit_transform(tcga), columns = tcga.columns, index = tcga.index)

#There is a large amount of missing data in GTEX (2594 genes). These genes are removed 
gtex.isna().sum()
gtex_clean = gtex.replace(-np.inf, np.nan).dropna()



# ### MATCH GENES
# the matching genes between the tumour and normal data are identified and saved to new files 
# this has to be run seperately for metabric and tcga 

#create a dictionary of only the er+ bc tumour datasets
bc = {"metabric" : metabric, "tcga" : tcga_imputed}

#save datasets with matching genes to new dataframe
#metabric: 17903 genes, tcga: 18406 genes 
match_tumour_to_normal_data(bc, gtex_clean) 


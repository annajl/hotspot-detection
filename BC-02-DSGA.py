# -*- coding: utf-8 -*-
"""

This script performs Disease Specific Genomic Analysis (DSGA) on both tumour datasets independently. The disease component of diseased tissue data is estimated by comparison by a healthy state model. In our analysis, the GTEX data is used to build the healthy state model.

For the METABRIC data, a gene threshold is used when performing DSGA. The output is a transformed dataset representing the disease component (DcT) of the tumour samples for 575 genes with a significant deviation from the healthy data. 

For the TCGA data, a gene threshold is not used when performing DSGA. The genes in the TCGA DcT are restricted to those found in the METABRIC DcT. 

"""


import hot_mapper as hm
import pandas as pd


#### SET UP EXPERIMENT
project_directory = "..."
output_path = f"{project_directory}/output/processed_data"


#### READ IN FILES
metabric_t = pd.read_csv(filepath_or_buffer = F"{output_path}/metabric_tumour_matched.csv", sep = ",", index_col = 0)
metabric_n = pd.read_csv(filepath_or_buffer = F"{output_path}/metabric_normal_matched.csv", sep = ",", index_col = 0)
tcga_t = pd.read_csv(filepath_or_buffer = F"{output_path}/tcga_tumour_matched.csv", sep = ",", index_col = 0)
tcga_n = pd.read_csv(filepath_or_buffer = F"{output_path}/tcga_normal_matched.csv", sep = ",", index_col = 0)


#### RUN DSGA
# DSGA is performed to build the disease component transformation (dct) of a tumour dataset
# A gene threshold is used on the METABRIC dct to those genes with a significant deviation from normal tissue
metabric_dct = hm.DSGA_transformation.DSGA(df_normal = metabric_n, df_tumour = metabric_t, threshold = True) #575 genes

# DSGA is performed on TCGA but a gene threhold is not used
tcga_dct = hm.DSGA_transformation.DSGA(df_normal = tcga_n, df_tumour =  tcga_t, threshold = False) #18406

# Instead genes in TCGA  are matched to those in the METABRIC DcT
gene_list = metabric_dct.columns
tcga_dct_T = tcga_dct.T
tcga_dct_thres = tcga_dct_T[gene_list] #575 genes
print(f"TCGA genes matched to METABRIC genes: {tcga_dct_thres.shape}")


#### SAVE OUTPUT
metabric_dct.to_csv(F"{output_path}/metabric_dct.csv") 
tcga_dct_thres.to_csv(F"{output_path}/tcga_dct.csv") 

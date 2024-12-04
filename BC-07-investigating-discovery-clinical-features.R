#' In this script Chi-Squared tests are performed to test the associations between clinical categorical features 
#' and the hotspot group. Outputs are saved in "clinical_feature_investigation". 
#' The output contains sheets of the group contigency tables and 
#' chi2 results are in the "CHISQ" sheet at the end. 





library(openxlsx)
library(janitor)
library(tidyverse)

#### FUNCTIONS
#for each column create a dataframe of patient split and percentages 
split_categories <- function(column){
  hotspot_labels <- clinical_test |> 
    select(hotspot)|> 
    pull(hotspot)
  #save the table distribution of individuals 
  tbl <- table(column,hotspot_labels) 
  cat_df <- as.data.frame.matrix(tbl)
  
  #plus the percentage of the count
  #for each row take the percentage 
  cat_df_prop <- apply(cat_df, 2, function(x) round(100*prop.table(x),digits=2))
  colnames(cat_df_prop) <- c("0_%", "1_%")
  
  #combnine both dataframes 
  final_cat <- cbind(cat_df,cat_df_prop)
  return(final_cat)
}



#### READ IN FILES
setwd("...")
survival <- read_csv("./output/processed_data/metabric_hotspot_id_survival.csv") |>
  rename("id" = "...1")

clinical <- read_csv("./data/metabric_clinical_historic_version.csv") |> 
  clean_names() |>
  rename("id" = "x1")


#### CHI SQUARE TEST 
#seelct the correct patients 
clinical_subset <- clinical[clinical$id %in% survival$id,] |> 
  arrange(id)

# Run chi-square tests on the relevant clinical features against hotspot groups
# use clinical features suggested by Nick
test_factor <- c("chemotherapy", "intclust", "her2_snp6", "inferred_menopausal_state", "claudin_subtype", "threegene", "histological_subtype" )

clinical_test <- clinical_subset |>
  select(all_of(test_factor)) 


#divide clinical into hotspot and neighbourhood
clinical_test$hotspot <- as.factor(survival$hotspot)
summary(clinical_subset)

#build contigency tables 
cont_list <- apply(clinical_test[,test_factor], 2, function(x) split_categories(x))

#run chi-squared for every categorical column
chisq.stat <- apply(clinical_test[,test_factor], 2, function(x) chisq.test(x, clinical_test$hotspot)$statistic)
chisq.pval <- apply(clinical_test[,test_factor], 2, function(x) chisq.test(x, clinical_test$hotspot)$p.value)


#combine to single dataframe 
chisq.df <- cbind(chisq.stat,chisq.pval)
colnames(chisq.df) <- c("x_squared","p_val")

#### SAVE DISTIRBUTION OF CLINICAL FEATURES 
#add to list
cont_list$CHISQ <- chisq.df
cont_list

#write up excel file of lists
write.xlsx(cont_list, file = "./output/clinical_feature_investigation/clinical_features_contingency_tables.xlsx",rowNames=T)



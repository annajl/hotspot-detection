#' This script selects which type of survival information (e.g. overall survival or relapse free survival) 
#' will be used for this analysis from the clinical dataset.



library(tidyverse)
library(janitor)


#### READ IN FILES
setwd("...") # project folder

metabric_clinical <- read_csv("./data/metabric_clinical_erpos_cbioportal.csv") |> clean_names()
tcga_clinical <- read_csv("./data/tcga_clinical_erpos_cbioportal.csv") |> clean_names()


#### SELECT SURVIVAL INFORMATION

# look at the surival information available for both datasets
metabric_clinical |> 
  select(patient_id, relapse_free_status, relapse_free_status_months,overall_survival_status, overall_survival_months) |>
  glimpse()

tcga_clinical |> 
  select(patient_id_1, disease_free_status, disease_free_months,overall_survival_status, overall_survival_months) |>
  glimpse()


# Relapse free survival is used for the METABRIC data
metabric_survival <- metabric_clinical |>
  select(patient_id, relapse_free_status, relapse_free_status_months) |> # select relevant columns
  rename("event" = "relapse_free_status", "time" = "relapse_free_status_months") |> # rename to time and event
  mutate(event = recode(event, '0:Not Recurred' = 0, '1:Recurred' = 1)) # code categories as binary 

write.csv(metabric_survival, "./output/processed_data/metabric_survival.csv")


# Overall survival is used for the TCGA data
tcga_survival <- tcga_clinical |>
  select(patient_id_1, disease_free_status, disease_free_months) |> # select relevant columns
  rename("patient_id" = "patient_id_1", "event" = "disease_free_status", "time" = "disease_free_months") |> # rename to time and event
  mutate(event = recode(event, '0:DiseaseFree' = 0, '1:Recurred/Progressed' = 1)) # code categories as binary 


write.csv(tcga_survival, "./output/processed_data/tcga_survival.csv")

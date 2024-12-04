#' Perform survival analysis on METABRIC hotspot group and plot Kaplan-Meier graph. 
#' Kaplan-Meier plot needs to be manually saved due to issues with 'survminer' package




library(survminer)
library(survival)
library(tidyverse)

#### READ IN FILES
setwd("...")
survival <- read_csv("./output/processed_data/tcga_hotspot_id_survival.csv")


#looking at 10 year survival
t <- 120
survival$event[survival$time > t] <- 0
survival$Time[survival$time > t] <- t
table(survival$hotspot)
dim(survival)[1]



#### KAPLAN-MEIER ANALYSIS
#fit kaplan meier model
km_fit <- survfit(Surv(time, event) ~ hotspot, data = survival)

#summary of results
survdiff(Surv(time, event) ~ hotspot, data = survival)
table(survival$hotspot)


# manually save plot
km_plot <- ggsurvplot(km_fit, 
                      data = survival, 
                      pval = T,
                      mark.time = T,
                      xlim = c(0,120),
                      break.time.by = 24,
                      xlab = "Months", 
                      ylab = "Relapse free probability",
                      risk.table = TRUE,
                      conf.int = TRUE,# label curves directly
                      legend.labs =  c("Neighbourhood","Hotspot")) # legend instead of direct label)


km_plot




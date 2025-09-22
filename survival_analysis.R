rm(list=ls())
library(data.table)
library(ggplot2)
library(car)
library(survival)
library(survminer)
library(forestplot)
library(openxlsx)
library(data.table)


dataset= read.xlsx("/Users/cfungtammasan/Google Drive/My Drive/Chai_Shruti/WGD/2025Mar17/2168case_database.xlsx", sheet = 1)
setDT(dataset)
dataset[, WGD := fifelse(`Ploidy` >= 2.7, 
                         "positive", 
                         fifelse(`Ploidy` < 2.7, "negative", NA_character_))]
dataset[, WGD := factor(WGD, levels = c("negative", "positive"))]

dataset=dataset[!is.na(WGD)]
dim(dataset) # 2155   99

dataset[, RAS := factor(RAS, levels = c("WT", "MUT"))]
table(dataset$RAS,useNA='always')

dataset[, BRAF.V600E := factor(BRAF.V600E, levels = c("WT", "MUT"))]
dataset[, MSI := factor(MSI, levels = c("MSS", "MSI-High"))]
dataset[, PrimSite := factor(PrimSite, levels = c("Right-sided colon", "Left-sided colon","Rectum"))]
dataset[,Sex:=Gender]

## 4a
surv_object = Surv(time=dataset$DFS, event=dataset$DFS.Event)
fit=survfit(surv_object~WGD,data=dataset)
ggsurvplot(
  fit,
  data=dataset,
  pval=T,
  conf.int=T,
  risk.table=T,
  risk.table.col="strata",
  ggtheme=theme_minimal()
)
cox_model=coxph(surv_object~WGD ,data=dataset)
summary(cox_model)


## 4b
surv_object = Surv(time=dataset$DFS, event=dataset$DFS.Event)
cox_model=coxph(surv_object~WGD + Sex + Age + Stage + PrimSite + MSI + RAS + BRAF.V600E + ctDNA.Baseline + ctDNA.MRD ,data=dataset)
summary(cox_model)
vif(cox_model)
ggforest(cox_model,data=dataset)


## 5a
dataset$WGD_RAS=factor(interaction(dataset$WGD,dataset$RAS))
table(dataset$RAS,useNA='always')
table(dataset$WGD_RAS)
fit=survfit(surv_object~WGD_RAS,data=dataset)

ggsurvplot(
  fit,
  data=dataset,
  pval=T,
  conf.int=T,
  risk.table=T,
  legend.labs=c("WGD-/RAS_WT","WGD+/RAS_WT","WGD-/RAS_MUT","WGD+/RAS_MUT"),
  risk.table.col="strata",
  ggtheme=theme_minimal(),
  xlim=c(0,1400)
)
cox_model=coxph(surv_object~WGD_RAS ,data=dataset)
ggforest(cox_model)


## S5b
MSS_dataset=dataset[MSI=='MSS']
dim(MSS_dataset) #1950
surv_object2 = Surv(time=MSS_dataset$DFS, event=MSS_dataset$DFS.Event)

cox_model2=coxph(surv_object2~WGD + Sex + Age + Stage + PrimSite  + RAS + BRAF.V600E + ctDNA.Baseline + ctDNA.MRD ,data=MSS_dataset)
ggforest(cox_model2,data=MSS_dataset)
surv_object2 = Surv(time=MSS_dataset$DFS, event=MSS_dataset$DFS.Event)
fit2=survfit(surv_object2~WGD,data=MSS_dataset, conf.int=0.95,conf.type="plain")

ggsurvplot(
  fit2,
  data=MSS_dataset,
  pval=T,
  conf.int=T,
  risk.table=T,
  risk.table.col="strata",
  ggtheme=theme_minimal()
)
cox_model=coxph(surv_object2~WGD ,data=MSS_dataset)
ggforest(cox_model)


## 5b
fit=survfit(surv_object2~WGD_RAS,data=MSS_dataset)

ggsurvplot(
  fit,
  data=MSS_dataset,
  pval=T,
  conf.int=T,
  risk.table=T,
  legend.labs=c("WGD-/RAS_WT","WGD+/RAS_WT","WGD-/RAS_MUT","WGD+/RAS_MUT"),
  risk.table.col="strata",
  ggtheme=theme_minimal(),
  xlim=c(0,1400)
)
cox_model=coxph(surv_object2~WGD_RAS ,data=MSS_dataset)
ggforest(cox_model)


## 4c-4d
ctDNA_data=dataset[!is.na(`ctDNA.MRD`)]
dim(ctDNA_data) #2023
# surv_object = Surv(time=ctDNA_data$OS, event=ctDNA_data$OS.Event)
surv_object = Surv(time=ctDNA_data$DFS, event=ctDNA_data$DFS.Event)

ctDNA_data$MRD_WGD=factor(interaction(ctDNA_data$ctDNA.MRD,ctDNA_data$WGD))

fit=survfit(surv_object~MRD_WGD,data=ctDNA_data)

ggsurvplot(
  fit,
  data=ctDNA_data,
  pval=T,
  conf.int=T,
  risk.table=T,
  risk.table.col="strata",
  legend.labs=c("MRD-/WGD-","MRD+/WGD-","MRD-/WGD+","MRD+/WGD+"),
  ggtheme=theme_minimal(),
  xlim=c(0,1400)
)



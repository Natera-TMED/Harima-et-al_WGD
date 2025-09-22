rm(list=ls())
library(data.table)
library(ggplot2)
library(GGally)
library(car)
library(Hmisc)
library(gmodels)
library(pscl)
library(survival)
library(survminer)
library(ggfortify)
library(visdat)
library(forestplot)
library(dplyr)
library(broom)
library(gtsummary)
library(officer)
library(flextable)

#####################
## Helper function ##
#####################

forest_plot <- function(model){
  tidy_model=tidy(model,conf.int=T)
  tidy_model=as_tibble(tidy_model)
  tidy_model <- tidy_model %>%
    filter(term !="(Intercept)")
  tidy_model <- tidy_model %>%
    mutate(
      p_value_label=sprintf("P = %.2e", p.value)
    )
  forest_plot <- ggplot(tidy_model,aes(x=term,y=estimate)) +
    geom_point()+
    geom_errorbar(aes(ymin=conf.low,ymax=conf.high),width=0.2)+
    geom_text(aes(label=p_value_label),vjust=-0.5,size=3)+
    geom_hline(yintercept=0,linetype="dashed")+
    coord_flip()+
    xlab("Predictor")+
    ylab("Log Odds Ratio (95% CI)")+
    theme_minimal()
  print(forest_plot)
}

#######################################
## Read data and transform attribute ##
#######################################

dataset=fread("George_WGD20240402.tsv")
dataset[,23]=NULL

MSI_data=fread('MSIstatus20250210.xlsx - Sheet1.csv')
dim(dataset)
length(unique(dataset$PtsID))
dim(MSI_data)
MSI_data=MSI_data[,.(PtsID,MSI)]
dataset=MSI_data[dataset,on='PtsID']
names(dataset)
table(dataset[,.(MSI,i.MSI)]) # in current version, MSI and i.MSI are equal

names_data=names(dataset)
names_data
dim(dataset)
names_data[41]="new_RAS"
names_data[42]="new_BRAF.V600E"
names_data[44]="new_MSI" # in current version, MSI and new_MSI are equal
names(dataset)=names_data
names(dataset)
vis_dat(dataset)
table(dataset$new_MSI)
table(dataset$MSI)

summary(dataset)
describe(dataset)
describe(dataset$new_MSI)

dataset[,range(`copy number`),by=WGD]

############
## cutoff ##
############
dataset[, WGD := fifelse(`ploidy` >= 2.7, 
                         "positive", 
                         fifelse(`ploidy` < 2.7, "negative", NA_character_))]

ggplot(dataset, aes(x = ploidy)) +
  geom_histogram(bins = 200, fill='gray') +
  # geom_density()+
  geom_vline(xintercept = 2.7, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(floor(min(dataset$ploidy, na.rm=TRUE)),
                                  ceiling(max(dataset$ploidy, na.rm=TRUE)),
                                  by = 0.5))+
  labs(x = "Ploidy", y = "Number of patients")


table(dataset$WGD)

#############################
## add missing data marker ##
#############################

dataset[Stage=="荳肴・/Unknown",Stage:=NA]
dataset[Stage=="",Stage:=NA]
dataset[Stage=='0',Stage:=NA]
dataset[PrimSite=="",PrimSite:=NA]
dataset[MSI=="",MSI:=NA]
dataset[RAS=="",RAS:=NA]
dataset[BRAF.V600E=="",BRAF.V600E:=NA]
dataset[new_MSI=="",new_MSI:=NA]
dataset[new_RAS=="",new_RAS:=NA]
dataset[new_BRAF.V600E=="",new_BRAF.V600E:=NA]
dataset[ctDNA.MRD=="",ctDNA.MRD:=NA]
dataset[ctDNA.Baseline=="",ctDNA.Baseline:=NA]


##########################
## making table 1 and 2 ##
##########################
# dt.pat.table.source = dataset[, .(Sex=Gender,Age,Stage,`Tumor Location`=PrimSite,MSI,RAS,BRAF.V600E,ctDNA.MRD,NAC,ACT)]
dt.pat.table.source = dataset[, .(WGD,Sex=Gender,Age,Stage,`Tumor Location`=PrimSite,MSI,RAS,BRAF.V600E,ctDNA.MRD,NAC,ACT)]
dt.pat.table.source$WGD=factor(dt.pat.table.source$WGD,levels=c('positive','negative'))
table(dt.pat.table.source$NAC,useNA='always')

reset_gtsummary_theme()
theme_gtsummary_compact()

table <- dt.pat.table.source |>
  tbl_summary(
  )   
table

table <- dt.pat.table.source |>
  tbl_summary(
    by=WGD
  )   
table


##################
## order factor ##
##################
dataset$WGD=factor(dataset$WGD,level=c("positive","negative"))
dataset$Gender=factor(dataset$Gender,level=c("Male","Female"))
dataset$Stage=factor(dataset$Stage,level=c("I","II","III","IV"))
dataset$PrimSite=factor(dataset$PrimSite,level=c("Right-sided colon","Left-sided colon","Rectum"))
dataset$MSI=factor(dataset$MSI,level=c("MSS","MSI-HIGH"))
dataset$BRAF.V600E=factor(dataset$BRAF.V600E,level=c("WT","MUT"))
dataset$RAS=factor(dataset$RAS,level=c("WT","MUT"))
dataset$new_MSI=factor(dataset$new_MSI,level=c("MSS","MSI-High"))
dataset$new_BRAF.V600E=factor(dataset$new_BRAF.V600E,level=c("WT","MUT"))
dataset$new_RAS=factor(dataset$new_RAS,level=c("WT","MUT"))
dataset$DFS.Event.survival=dataset$DFS.Event
dataset$DFS.Event=factor(dataset$DFS.Event,level=c("TRUE","FALSE"))
dataset$ctDNA.MRD=factor(dataset$ctDNA.MRD,level=c("POSITIVE","NEGATIVE"))
dataset$ctDNA.Baseline=factor(dataset$ctDNA.Baseline,level=c("POSITIVE","NEGATIVE"))
dataset[TMB=="",TMB:=NA]

## simplify variable name
dataset$MSI=dataset$new_MSI
dataset$RAS=dataset$new_RAS
dataset$BRAF.V600E=dataset$new_BRAF.V600E

########################
## check missing data ##
########################
table(dataset$new_RAS)
vis_dat(dataset)
vis_miss(dataset)

###############################
## p-value for table 1 and 2 ##
###############################

for (i in c("Gender","Stage","PrimSite","MSI","RAS","BRAF.V600E","DFS.Event","ctDNA.MRD","ctDNA.Baseline","TMB")){
  print(i)
  WGD_t=table(dataset[[i]],dataset$WGD)
  print(WGD_t)
  print(round(100*prop.table(WGD_t,margin=2),digits=1))
  # print(chisq.test(WGD_t))
  print(fisher.test(WGD_t))
  print("")
}
dataset[,median(Age),by=dataset$WGD]

t.test(dataset[WGD=='positive',Age],dataset[WGD=='positive',Age])
)

##################
## ctDNA vs WGD ##
##################
dataset$WGD=factor(dataset$WGD,level=c("negative","positive"))
dataset$DFS.Event=factor(dataset$DFS.Event,level=c("FALSE","TRUE"))
dataset$ctDNA.MRD=factor(dataset$ctDNA.MRD,level=c("NEGATIVE","POSITIVE"))
dataset$ctDNA.Baseline=factor(dataset$ctDNA.Baseline,level=c("NEGATIVE","POSITIVE"))



## Sup 1
ori_model=glm(ctDNA.MRD ~ WGD + MSI+ Stage ,data=dataset,family="binomial"(link=logit))
forest_plot(ori_model)
summary(ori_model)
vif(ori_model)




##########CODE MODIFIED FROM STEPHANIE TULLO#################

library(RMINC)
library(lme4)
library(lmerTest)
library(ggplot2)
library(carData)
library(effects)
library(grid)
library(MRIcrotome)
library(splines)
library("RColorBrewer")
library(tidyverse)
library(magrittr) #to be able to use "%>%"

mbmdata <- read.csv("../preprocessed/steph_dbm/tp-data-DBM-analysis-202202-symptom.csv")

setwd("/data/chamal/projects/stephanie/long-asyn-PFF-project/analysis/2level_analyses_part6/")

######################START WITH CAUDOPUTAMEN##############################

###load data (most recent version)####

mask = "/data/chamal/projects/stephanie/long-asyn-PFF-project/derivatives/2level_model_niagara_202202/mask/mask_on_template_corrected.mnc"

anatFile <- "/data/chamal/projects/stephanie/long-asyn-PFF-project/derivatives/2level_model_niagara_202202/secondlevel_template0.mnc"
anatVol <- mincArray(mincGetVolume(anatFile))

#get rid of mice that did not have DBM file
mbmdata = subset(mbmdata, mbmdata$rel_log_nlin_det_blur!="")
for (row in 1:nrow(mbmdata)) {
  if(file.exists(toString(mbmdata$rel_log_nlin_det_blur[row]))==FALSE){
    # print(toString(mbmdata$rel_log_nlin_det_blur[row]))
    print(row)
    mbmdata =mbmdata[-c(row),]
    
  } 
}

#get rid of mice that did not receive an injection
mbmdata = subset(mbmdata, mbmdata$injection!="")
mbmdata = subset(mbmdata, mbmdata$injection!="NA")
mbmdata = subset(mbmdata, mbmdata$genotype!="homo")

#fix spreadsheet
mbmdata$weight[mbmdata$weight=="??"] <- NA
mbmdata$weight[mbmdata$weight==""] <- NA
mbmdata$litter_size[mbmdata$litter_size=="4/8"] <- 8

#reduce sampling -- remove cohorts with only two tps
cohortNum=c(2,3,5,6,20,22,23)
for (i in 1:length(cohortNum)){  
  mbmdata = subset(mbmdata, mbmdata$cohort!=cohortNum[i])
}

#set ref level
mbmdata$injection=as.factor(mbmdata$injection)
mbmdata$injection = relevel(mbmdata$injection, ref = "PBS")
mbmdata$sex=as.factor(mbmdata$sex)
mbmdata$sex = relevel(mbmdata$sex, ref = "F")
mbmdata$genotype=as.factor(mbmdata$genotype)
mbmdata$genotype = relevel(mbmdata$genotype, ref = "WT")
mbmdata$animal_ID = as.factor(mbmdata$animal_ID)
mbmdata$cohort = as.factor(mbmdata$cohort)
mbmdata$batch = as.factor(mbmdata$batch)
mbmdata$weight=as.numeric(as.character(mbmdata$weight))
mbmdata$X.dpi=as.numeric(as.character(mbmdata$X.dpi))
mbmdata$litter_size=as.numeric(as.character(mbmdata$litter_size))


# #count number of subjects
# nrow(subset(mbmdata, timepoint=="4"&sex=="F"))
# 
# for (i in 1:4){
#   print(nrow(subset(mbmdata, timepoint==i&sex=="M"&injection=="PBS"&genotype=="WT")))
# }

###############within sex differences##############
######################M83 PFF vs WT PBS########

#####Goal: simple t-statistic (comparable with regional-level analyses); no need to model subject as a random effect because we're looking at this cross-sectionally###

print(unique(mbmdata$X.dpi))
print(unique(mbmdata$genotype))
#print(unique(mbmdata$rel_log_nlin_det_blur))
print(unique(mbmdata$injection))

print(subset(mbmdata, genotype == "hemi" & between(X.dpi, 85, 95)))

setwd("/data/chamal/projects/natvik/sir_voxel_backup_trillium_03202026/analysis/")

MLM_cp <- mincLm(rel_log_nlin_det_blur ~  injection, subset(mbmdata, genotype=="hemi" & between(X.dpi, 85, 95)), mask = mask, parallel = c("slurm", 100))

mincWriteVolume(MLM_cp, "../derivatives/voxel_atrophy_maps/CP/M83HuPff_vs_M83PBS_t_stats_linear.mnc", column = 'tvalue-injectionHu-PFF')
mincWriteVolume(MLM_cp, "../derivatives/voxel_atrophy_maps/CP/M83MsPff_vs_M83PBS_t_stats_linear.mnc", column = 'tvalue-injectionMs-PFF')

###TODO: Convert to Allen template space (CCFv3, aligned to RAS)
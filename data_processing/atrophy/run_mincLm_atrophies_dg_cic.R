##############CODE MODIFIED FROM JANICE PARK####################
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

janice_dir <- "/data/scratch2/janice/alt_inj_site/analysis/final-dbm/"
setwd(janice_dir)

#LOAD DATA
dems <- read.csv(file = "demographics_edited.csv")
rel_jac <- read.csv(file = "rel_jac.csv", header = TRUE)

data <- merge(rel_jac,dems,by = c("ID","Timepoint"))
#set mask (converted to minc)
mask <- paste0(janice_dir, "new_mask.mnc")

#Set Male and control as default
data$Sex = relevel(as.factor(data$Sex), ref = "M")
data$Treatment = relevel(as.factor(data$Treatment), ref = "PBS")

##load template brain and its mask
anatVol <- mincArray(mincGetVolume("template_sharpen_shapeupdate.mnc"))


###############within sex differences##############
######################M83 PFF vs WT PBS########

#####Goal: simple t-statistic; no need to model subject as a random effect because we're looking at this cross-sectionally###

print(subset(data, Timepoint == 2))
setwd("/data/chamal/projects/natvik/sir_voxel_backup_trillium_03202026/analysis/")
MLM_cp <- mincLm(Relative_Jacobians ~  Treatment, subset(data, Timepoint==2), mask = mask, parallel = c("slurm", 100))

mincWriteVolume(MLM_cp, "../derivatives/voxel_atrophy_maps/DG/M83HuPff_vs_M83PBS_t_stats_linear.mnc", column = 'tvalue-TreatmentPFF')

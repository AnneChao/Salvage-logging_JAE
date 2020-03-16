library(xtable)
library(data.table)
library(reshape2)
library(cluster)
library(knitr, quietly = T)
library(ape)
library(ade4)
library(kableExtra)
library(dplyr)
library(phytools)
library(phyclust)
library(ggplot2)
tree_bayerwald <- read.nexus("Data/Trees_phylo/bayerwald_birds.tre")


bayerwald <- read.csv('Data/Abundance matrix/bayerwald_birds.csv', sep = ';')



plots=read.csv("Data/plots.csv", sep = ";")


bayerwald <- bayerwald[, c(7, 8, 10:ncol(bayerwald))]
bayerwald <- bayerwald %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
bayerwald <- bayerwald[, c(-1, -2)]




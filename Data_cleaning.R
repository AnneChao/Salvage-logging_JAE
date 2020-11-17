library(xtable)
library(knitr, quietly = T)
library(ape)
library(ade4)
library(kableExtra)
library(dplyr)
library(phytools)
library(phyclust)
tree_bayerwald <- read.nexus("Data/Trees_phylo/bayerwald_birds.tre")
tree_fontaine <- read.nexus("Data/Trees_phylo/oregon_birds_fontaine.tre")
tree_montana_hutto <- read.nexus("Data/Trees_phylo/montana_hutto_birds.tre")
tree_spain_castro <- read.nexus("Data/Trees_phylo/spain_castro_birds.tre")
tree_spain_rost <-  read.nexus("Data/Trees_phylo/spain_rost_birds.tre")
tree_zmihorski <-  read.nexus("Data/Trees_phylo/poland_zmihorski_birds.tre")
tree_cahall <- read.nexus("Data/Trees_phylo/oregon_cahall_birds.tre")
tree_lee <- read.nexus("Data/Trees_phylo/samchuck_lee_birds.tre")
tree_choi <- read.nexus("Data/Trees_phylo/samchuk_choi_birds_1.tre")
# trait <- read.csv("Traits/traits.csv", sep = ';')
# trait <- trait[, -1]
# # trait <- trait %>% arrange(species)
# rownames(trait) <- trait[, 1]
bayerwald <- read.csv('Data/Abundance matrix/bayerwald_birds.csv', sep = ';')
fontaine <- read.csv('Data/Abundance matrix/oregon_birds_fontaine.csv', sep = ';')
montana_hutto <- read.csv('Data/Abundance matrix/montana_hutto_birds.csv', sep = ';')
spain_castro <- read.csv('Data/Abundance matrix/spain_castro_birds.csv', sep = ';')
spain_rost <- read.csv('Data/Abundance matrix/spain_rost_birds.csv', sep = ';')
zmihorski <- read.csv('Data/Abundance matrix/poland_zmihorski_birds.csv', sep = ';')
cahall <- read.csv('Data/Abundance matrix/oregon_cahall_birds.csv', sep = ';')
lee <- read.csv('Data/Abundance matrix/samchuck_lee_birds.csv', sep = ';')
choi <- read.csv('Data/Abundance matrix/samchuk_choi_birds_1.csv', sep = ';')


plots=read.csv("Data/plots.csv", sep = ";")

# colnames(bayerwald)
bayerwald <- bayerwald[, c(7, 8, 10:ncol(bayerwald))]
bayerwald <- bayerwald %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
bayerwald <- bayerwald[, c(-1, -2)]
# sum(bayerwald[, -ncol(bayerwald)]%%1>0)

fontaine <- fontaine[, c(7, 8, 10:ncol(fontaine))]
fontaine <- fontaine %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
fontaine <- fontaine[, c(-1, -2)]
# sum(fontaine[, -ncol(fontaine)]%%1>0)

montana_hutto <- montana_hutto[, c(7, 8, 10:ncol(montana_hutto))]
montana_hutto <- montana_hutto %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
montana_hutto <- montana_hutto[, c(-1, -2)]
# sum(montana_hutto[, -ncol(montana_hutto)]%%1>0)

spain_castro <- spain_castro[, c(7, 8, 10:ncol(spain_castro))]
spain_castro <- spain_castro %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
spain_castro <- spain_castro[, c(-1, -2)]
# sum(spain_castro[, -ncol(spain_castro)]%%1>0)

spain_rost <- spain_rost[, c(7, 8, 10:ncol(spain_rost))]
spain_rost <- spain_rost %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
spain_rost <- spain_rost[, c(-1, -2)]
# sum(spain_rost[, -ncol(spain_rost)]%%1>0)

zmihorski <- zmihorski[, c(7, 8, 10:ncol(zmihorski))]
zmihorski <- zmihorski %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
zmihorski <- zmihorski[, c(-1, -2)]
# sum(zmihorski[, -ncol(zmihorski)]%%1>0)

cahall <- cahall[, c(7, 8, 10:ncol(cahall))]
cahall <- cahall %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
cahall <- cahall[, c(-1, -2)]

lee <- lee[, c(7, 8, 10:ncol(lee))]
lee <- lee %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
lee <- lee[, c(-1, -2)]

choi <- choi[, c(7, 8, 10:ncol(choi))]
choi <- choi %>% mutate(jahrflaeche=paste0(year_after_dist, ifelse(treatment=="salvaged", "FAE_1", "FKN_17")))
choi <- choi[, c(-1, -2)]





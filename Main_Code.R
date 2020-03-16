#Sourcing data 
source(file = 'Data_cleaning_incomplete.r')
#Sourcing functions for Phylogenetic disssimilarity 
source("Phylogenetic.r")
#Sourcing functions for taxonomic and functional disssimilarity 
source("new.R")

#Basic parameters setting
q = c(0,1,2)
conf = 0.95
nms <- c("bayerwald")
#### Example for German windstorm data (Thorn et al. 2016); see Table 1 of Georgiev et al. (2020).

#==========Taxonomic analysis=============
data_fd = list(
  bayerwald2 = data.frame(jahrflaeche = bayerwald$jahrflaeche, bayerwald[,-ncol(bayerwald)]>0)
) %>% lapply(., function(x){ process_data2(x, plots) }) 

td_diss <- lapply(1:length(nms),function(i){
  print(i)
  a <- data_fd[[i]]
  out_ <- Tax_diss(a$dat, a$mat)$output
})

names(td_diss) <- nms
all <- lapply(1:length(nms), function(i){
  # use filter(variable == "Jaccard") to plot Jaccard type instead.
  tmp <- td_diss[[i]] %>% filter(variable == "Sorensen") %>% mutate(site = nms[i])
}) %>% do.call(rbind,.)
all$year <- as.numeric(all$year)
all$site <- factor(all$site,levels = nms)
all$year[all$q=="q = 0"] <- all$year[all$q=="q = 0"] + 0.2
all$year[all$q=="q = 1"] <- all$year[all$q=="q = 1"] + 0.1

ggplot(all, aes(x = year, y = value)) +  theme_bw() + 
  geom_point( aes(color = q, shape = q),size = 2.5)+
  geom_errorbar(aes(ymin=UCL, ymax=LCL), width=0.3,alpha = 0.5)+
  facet_wrap(.~site,scales = "free_y")+ 
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+
  scale_shape_manual(values=c(16, 18, 17))+
  scale_color_manual(values=c("#009E73", "#D55E00", "darkorchid3"))



#==========Phylogeny Analysis=============
#Analyze all nine datas at the reference time 99.8305, which is the maximal height among nine phylogeny trees.

pd_diss <- lapply(1:length(nms),function(i){
  print(i)
  dat <- get(nms[i])
  tre <- get(paste0("tree_",nms[i]))
  out_ <- PHD.year(data = dat, tree = tre, plot = plots,B = 200,ref_t = 99.8305)[[1]]
  out_
})# Use B = 2 to save computation time for boostraps.
names(pd_diss) <- nms
all <- lapply(1:length(nms), function(i){
  # use filter(variable == "Jaccard") to plot Jaccard type instead.
  tmp <- pd_diss[[i]] %>% filter(variable == "Sorensen") %>% mutate(site = nms[i])
}) %>% do.call(rbind,.)
all$year <- as.numeric(all$year)
all$site <- factor(all$site,levels = nms)
all$year[all$q=="q = 0"] <- all$year[all$q=="q = 0"] + 0.2
all$year[all$q=="q = 1"] <- all$year[all$q=="q = 1"] + 0.1

ggplot(all, aes(x = year, y = value)) +  theme_bw() + 
  geom_point( aes(color = q, shape = q),size = 2.5)+
  geom_errorbar(aes(ymin=UCL, ymax=LCL), width=0.3,alpha = 0.5)+
  facet_wrap(.~site,scales = "free_y")+ 
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+
  scale_shape_manual(values=c(16, 18, 17))+
  scale_color_manual(values=c("#009E73", "#D55E00", "darkorchid3"))


#==========Fucntional Analysis=============
data_fd = list(
  bayerwald2 = data.frame(jahrflaeche = bayerwald$jahrflaeche, bayerwald[,-ncol(bayerwald)]>0)
) %>% lapply(., function(x){ process_data2(x, plots) }) 

#Read trait table and convert them into distance matrix
trait = read.csv('data/Traits/traits.csv', sep = ';')
rownames(trait) = trait$species
trait = trait[,-c(1,2)]
allsp <- sapply(data_fd, function(x){rownames(x$dat)[-1]}) %>% unlist() %>% unique() 
trait <- trait[rownames(trait) %in% allsp,]
dis <- as.matrix(daisy(trait, "gower", stand = T, weights = getWeightVector(trait), 
                       type = list(symm = getBinCol(trait))))
#Use common thresholds
knots = 50

sav <- lapply(data_fd, function(y){
  sp = rownames(y$dat)[-1]
  out <- y$dat[-1,] %>% as.data.frame() %>%
    select(which(grepl("unsalveged",colnames(.))==FALSE)) %>% cbind(sp,.)
  out$sp <- as.character(out$sp); out
}) %>% purrr::reduce(full_join, by = "sp")
sav[is.na(sav)] <- 0

# datatmp <- cbind(sav[,-1], uns[,-1]) %>% apply(.,2,function(i) i/sum(i)) %>% rowMeans()
dis <- dis[sav$sp,sav$sp]
# dmean <- sum( (datatmp %*% t(datatmp)) * dis)
taus <- c(seq(min(dis[dis>0]),max(dis),length.out = knots)) %>% sort()

fd_diss <- lapply(1:length(data_fd),function(i){
  print(i)
  a <- data_fd[[i]]
  tmp = colnames(dis) %in% rownames(a$dat)[-1]
  dis_match = dis[tmp,tmp]
  dis_match = dis_match[rownames(a$dat)[-1], rownames(a$dat)[-1]]
  out_<- FD_diss_AUC(dat = a$dat,mat = a$mat,dis = dis_match,aucboot = 200,taus_common = taus)
  out_
})# Use aucboot = 2 to save computation time for boostraps.
nms <- c("bayerwald")
names(fd_diss) <- nms
all <- lapply(1:length(nms), function(i){
  # use filter(variable == "Jaccard") to plot Jaccard type instead.
  tmp <- fd_diss[[i]] %>% filter(variable == "Sorensen") %>% mutate(site = nms[i])# use 
}) %>% do.call(rbind,.)
all$site <- factor(all$site,levels = nms)
all$year[all$q=="q = 0"] <- all$year[all$q=="q = 0"] + 0.2
all$year[all$q=="q = 1"] <- all$year[all$q=="q = 1"] + 0.1

ggplot(all, aes(x = year, y = value)) +  theme_bw() + 
  geom_point( aes(color = q, shape = q),size = 2.5)+
  geom_errorbar(aes(ymin=UCL, ymax=LCL), width=0.3,alpha = 0.5)+
  facet_wrap(.~site,scales = "free_y")+ 
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90))+
  scale_shape_manual(values=c(16, 18, 17))+
  scale_color_manual(values=c("#009E73", "#D55E00", "darkorchid3"))


all <- lapply(1:length(nms), function(i){
  tmp <- fd_diss[[i]] %>% filter(variable != "Sorensen") %>% select(q,year,value,LCL,UCL)
}) %>% do.call(rbind,.)

all_mean <- all %>% group_by(q,year) %>% summarise(value = mean(value), LCL = mean(LCL),UCL = mean(UCL))
ggplot(all_mean, aes(x = year, y = value)) +  theme_bw() + 
  geom_point(size = 2.5)+
  geom_errorbar(aes(ymin=UCL, ymax=LCL), width=0.3,alpha = 0.5)+
  facet_grid(q ~ .,scales = "free_y")+scale_x_continuous("year after disturbance", breaks=seq(1,17))+
  theme(legend.position = "bottom")+ylab("Functional Dissimilarity (Jaccard)")
  
  
  
  

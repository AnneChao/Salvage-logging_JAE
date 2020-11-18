#Sourcing data 
source(file = 'Data_cleaning.r')
#Sourcing functions for Phylogenetic disssimilarity 
source("Phylogenetic.r")
#Sourcing functions for taxonomic and functional disssimilarity 
source("new.R")

#For rost data, fix the "," issue
for (i in c(1:(ncol(spain_rost)-1))) {
  if(class(spain_rost[,i]) != "integer"){
    spain_rost[,i] = as.character(spain_rost[,i])
    tmp = grep(",",spain_rost[,i])
    spain_rost[tmp,i]="1"
    spain_rost[,i] = as.numeric(spain_rost[,i])
  }
}

#For montana_hutto, keep data years are less than 17 
myear <- as.numeric(substr(montana_hutto$jahrflaeche,start=1,
                           stop=regexpr("F",montana_hutto$jahrflaeche)-1 ))
montana_hutto = montana_hutto[myear<=17,]

#Basic parameters setting
q = c(0,1,2)
conf = 0.95
nms <- c("bayerwald","fontaine","montana_hutto","spain_castro","spain_rost", "zmihorski","cahall",
         "choi","lee")
#==========Taxonomic analysis=============
data_fd = list(
  bayerwald2 = data.frame(jahrflaeche = bayerwald$jahrflaeche, bayerwald[,-ncol(bayerwald)]>0),
  fontaine2 = data.frame(jahrflaeche = fontaine$jahrflaeche, fontaine[,-ncol(fontaine)]>0),
  montana_hutto2 = data.frame(jahrflaeche = montana_hutto$jahrflaeche, montana_hutto[,-ncol(montana_hutto)]>0),
  spain_castro2 = data.frame(jahrflaeche = spain_castro$jahrflaeche, spain_castro[,-ncol(spain_castro)]>0),
  spain_rost2 = data.frame(jahrflaeche = spain_rost$jahrflaeche, spain_rost[,-ncol(spain_rost)]>0),
  zmihorski2 = data.frame(jahrflaeche = zmihorski$jahrflaeche, zmihorski[,-ncol(zmihorski)]>0),
  cahall2 = data.frame(jahrflaeche = cahall$jahrflaeche, cahall[,-ncol(cahall)]>0),
  choi2 = data.frame(jahrflaeche = choi$jahrflaeche, choi[,-ncol(choi)]>0),
  lee2 = data.frame(jahrflaeche = lee$jahrflaeche, lee[,-ncol(lee)]>0)
) %>% lapply(., function(x){ process_data2(x, plots) }) 

td_diss <- lapply(1:9,function(i){
  print(i)
  a <- data_fd[[i]]
  out_ <- Tax_diss(a$dat, a$mat)$output
})

names(td_diss) <- nms
all <- lapply(1:8, function(i){
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

pd_diss <- lapply(1:9,function(i){
  print(i)
  dat <- get(nms[i])
  tre <- get(paste0("tree_",nms[i]))
  out_ <- PHD.year(data = dat, tree = tre, plot = plots,B = 200,ref_t = 99.8305)[[1]]
  out_
})# Use B = 2 to save computation time for boostraps.
names(pd_diss) <- nms
all <- lapply(1:8, function(i){
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
  bayerwald2 = data.frame(jahrflaeche = bayerwald$jahrflaeche, bayerwald[,-ncol(bayerwald)]>0),
  fontaine2 = data.frame(jahrflaeche = fontaine$jahrflaeche, fontaine[,-ncol(fontaine)]>0),
  montana_hutto2 = data.frame(jahrflaeche = montana_hutto$jahrflaeche, montana_hutto[,-ncol(montana_hutto)]>0),
  spain_castro2 = data.frame(jahrflaeche = spain_castro$jahrflaeche, spain_castro[,-ncol(spain_castro)]>0),
  spain_rost2 = data.frame(jahrflaeche = spain_rost$jahrflaeche, spain_rost[,-ncol(spain_rost)]>0),
  zmihorski2 = data.frame(jahrflaeche = zmihorski$jahrflaeche, zmihorski[,-ncol(zmihorski)]>0),
  cahall2 = data.frame(jahrflaeche = cahall$jahrflaeche, cahall[,-ncol(cahall)]>0),
  choi2 = data.frame(jahrflaeche = choi$jahrflaeche, choi[,-ncol(choi)]>0),
  lee2 = data.frame(jahrflaeche = lee$jahrflaeche, lee[,-ncol(lee)]>0)
) %>% lapply(., function(x){ process_data2(x, plots) }) 

#Read trait table and convert them into distance matrix
trait = read.csv('Data/Traits/traits.csv', sep = ';')
rownames(trait) = trait$species
trait = trait[,-c(1,2)]
allsp <- sapply(data_fd, function(x){rownames(x$dat)[-1]}) %>% unlist() %>% unique() 
trait <- trait[rownames(trait) %in% allsp,]
trait[,3] <- as.factor(trait[,3])
dis <- as.matrix(daisy(trait, "gower", stand = T, weights = getWeightVector(trait), 
                       type = list(symm = getBinCol(trait))))
#Use common thresholds
knots = 50
# uns <- lapply(data_fd, function(y){
#   sp = rownames(y$dat)[-1] 
#   out <- y$dat[-1,] %>% as.data.frame() %>% 
#     select(which(grepl("unsalveged",colnames(.))==TRUE)) %>% cbind(sp,.)
#   out$sp <- as.character(out$sp); out
#   }) %>% purrr::reduce(full_join, by = "sp") 
# uns[is.na(uns)] <- 0
# 
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

fd_diss <- lapply(1:9,function(i){
  print(i)
  a <- data_fd[[i]]
  tmp = colnames(dis) %in% rownames(a$dat)[-1]
  dis_match = dis[tmp,tmp]
  dis_match = dis_match[rownames(a$dat)[-1], rownames(a$dat)[-1]]
  out_<- FD_diss_AUC(dat = a$dat,mat = a$mat,dis = dis_match,aucboot = 200,taus_common = taus)
  out_
})# Use aucboot = 2 to save computation time for boostraps.
nms <- c("bayerwald","fontaine","montana_hutto","spain_castro","spain_rost", "zmihorski","cahall",
         "choi","lee")
names(fd_diss) <- nms
all <- lapply(1:8, function(i){
  tmp <- fd_diss[[i]] %>% filter(variable == "Sorensen") %>% mutate(site = nms[i])
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




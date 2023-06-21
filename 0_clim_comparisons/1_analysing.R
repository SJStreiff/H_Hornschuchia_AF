####################################################################################################
#  2023-05-03 sjrs 
# HORNSCHUCHIA and the Atlantic Forest
# 
# extract climate and do some comparisons between species
####################################################################################################
rm(list = ls())

# INIT
wd <- '~/Sync/1_Annonaceae/H_Hornschuchia_AF/0_clim_comparisons/' # the working dir
dat_file <- '~/Sync/1_Annonaceae/H_Hornschuchia_AF/Y_data/Hornschuchia_Bocagea_Trigynaea_occurence.csv'
dat1 <- read.table(file = dat_file, sep = ';', head = T)
pl <- c('ggplot2', 
        'raster',
        'vegan',
        'terra',
        'ape',
        'dplyr')
sapply(pl, require, character.only=TRUE)
source('~/Sync/serafins_functions.R')
setwd(wd)


load('new_data_extr.Rdata')
# if the Rdata does not exist, run this code:
# source('./0_extracting.R')
ab_dat$AF <- dat1$AF
colnames(ab_dat)[colnames(ab_dat) == 'name'] <- 'species'


point_count <- ab_dat %>% dplyr::count(ab_dat$species, name=NULL)

# these species only have less than 3 points! might have to drop them.
sp_too_rare <- point_count[point_count$n < 3,]

clean_dat <- ab_dat[ab_dat$species %notin% sp_too_rare$`ab_dat$species`,]


#####################################
# Choose abiotic columns to extract #
#####################################

#variables <- c('Bio12', 'Bio17', 'Bio07', 'clay', 'soil_ph')



# These are all layers we extracted data from.

# 'Bio01', 'Bio02', 'Bio03', 'Bio04', 'Bio05', 'Bio06', 'Bio07',
# 'Bio08', 'Bio09', 'Bio10', 'Bio11', 'Bio12', 'Bio13', 'Bio14',
# 'Bio15', 'Bio16', 'Bio17', 'Bio18', 'Bio19', 'Bio20', 'Bio21',
# 'bulk_dens', 'cation_exchange_capacity', 'vol_fraction_coarse',
# 'clay_prop', 'nitrogen', 'org_C_dens', 'org_C_stock', 'soil_pH',
# 'sand_prop', 'silt_prop', 'soil_org_carbon')



# for the time being only Hornschuchia
# w_dat <- na.omit(ab_dat[ab_dat$genus %in% c('Hornschuchia', 'Trigynaea', 'Bocagea'),
#                         c('species','Bio12', 'Bio17', 'Bio04', 'Bio06', 'clay_prop', 'soil_pH')])
# taxa <- w_dat[,'species']
# w_dat <- w_dat[,c('Bio12', 'Bio17', 'Bio04', 'Bio06', 'clay_prop', 'soil_pH')]

ordi_dat <- clean_dat[clean_dat$genus %in% c('Hornschuchia', 'Trigynaea', 'Bocagea'),
                   c('genus','species', 'AF', 'Bio01', 'Bio02', 'Bio03', 'Bio04', 'Bio05', 'Bio06', 'Bio07',
                     'Bio08', 'Bio09', 'Bio10', 'Bio11', 'Bio12', 'Bio13', 'Bio14',
                     'Bio15', 'Bio16', 'Bio17', 'Bio18', 'Bio19', 'Bio20', 'Bio21',
                     'bulk_dens', 'cation_exchange_capacity', 'vol_fraction_coarse',
                     'clay_prop', 'nitrogen', 'org_C_dens', 'org_C_stock', 'soil_pH',
                     'sand_prop', 'silt_prop', 'soil_org_carbon')]
ordi_dat <- na.omit(ordi_dat)

##########################################
# per species MEANS and SD               #
##########################################

for(i in 3:length(colnames(ordi_dat))){ # or for(i in 1:19){
  ggplot() +
    geom_density(data = ordi_dat, aes(x = ordi_dat[,i], group =genus, colour =genus)) +
    theme_classic() +
    xlab(colnames(ordi_dat)[i]) 
    #ggtitle(label = abc[i]) 
  print(i)
  #p2
  }


##########################################
# check different distances              #
##########################################

# #  
# w_dat <- na.omit(w_dat)
# 
# sp_vec <- w_dat[,"species"]
# 
# #nmds try 1
# t1 <- metaMDS(w_dat,
#               'bray',
#               na.rm = TRUE,
#               labels = sp_vec)
# # ordiplot(t1) %>% text('species') 
# require(vegan)
# # X-square
# 


##########################################
# PCA do figure out correlations in abiotic variables     #
##########################################



groupBY <- ordi_dat$genus
groupBY_AF <- ordi_dat$AF
pca_dat <- ordi_dat[,-1]
pca_dat <- pca_dat[,-1]
pca_dat <- pca_dat[,-1]
# this makes the correlation matrix between all variables of all genera data
correlation <- cor(pca_dat, method = 'pearson') # do a pearson correlation between all variables
library(corrplot)
corrplot(correlation,
         order = 'hclust', # order by hierarchical clusting
         addrect = 7) # add rectangle around clusters identified by hclust



#colnames(pca) <- nms
pc <- prcomp(pca_dat, scale = FALSE)
#pc <- prcomp(pc_dat, scale. = FALSE)
library(ggbiplot)
ggbiplot(pc, scale=FALSE, ellipse = TRUE, groups = groupBY) + ggtitle('PCA grouped by genus')

ggbiplot(pc, scale=TRUE, ellipse = TRUE, groups = groupBY_AF) + ggtitle('PCA grouped by biogeogr. region')
ggbiplot(pc, scale=TRUE, ellipse = TRUE, groups = groupBY, choices = c(3,4)) # pc axies 3 and 4

# just testing CUR decomposition. Might not be as relevant though...
cur_all <- dCUR::CUR(pca_dat,
                         variables =  colnames(pca_dat),
                         cur_method = 'mixture',
                         rows = 0.999,
                         columns = 0.3,
                         standardize = F)

red_pca_dat <- cur_all$C_cur
rpc <- prcomp(red_pca_dat, scale. = TRUE)
ggbiplot(rpc, scale = F, ellipse = T, groups = groupBY)
ggbiplot(rpc, scale = F, ellipse = T, groups = groupBY_AF)



#### Dimensionality reduction to identify the relevant abiotic variables.
# Doing this by genus, later maybe by presence/absence in AF

Horn_dat <- na.omit(ordi_dat[ordi_dat$genus == 'Hornschuchia',])
Trig_dat <- na.omit(ordi_dat[ordi_dat$genus == 'Trigynaea',])
Boca_dat <- na.omit(ordi_dat[ordi_dat$genus == 'Bocagea',])

#################
# just for Hornsch.
groupBY_H <- Horn_dat$species
pca_dat_H <- Horn_dat[,-1]
pca_dat_H <- pca_dat_H[,-1]



# this makes the correlation matrix between all variables of all genera data
cor_H <- cor(pca_dat_H, method = 'pearson') # do a pearson correlation between all variables
corrplot(cor_H,
         order = 'hclust', # order by hierarchical clusting
         addrect = 10) # add rectangle around clusters identified by hclust


#colnames(pca) <- nms
pc_H <- prcomp(pca_dat_H, scale. = TRUE)
ggbiplot(pc_H, scale=TRUE, ellipse = TRUE, groups = groupBY_H) + ggtitle('PCA for Hornschuchia')
ggbiplot(pc_H, scale=TRUE, ellipse = TRUE, groups = groupBY_H, choices = c(3,4)) # pc axies 3 and 4
# as seen before and in other analyses, H.bryotrophe spans across many other species values...

cur_H <- dCUR::CUR(pca_dat_H,
                     variables =  colnames(pca_dat),
                     cur_method = 'mixture',
                     rows = 0.999,
                     columns = 0.3,
                     standardize = T)

red_pca_dat_H <- cur_H$C_cur
rpc_H <- prcomp(red_pca_dat_H, scale. = TRUE)
ggbiplot(rpc_H, scale = F, ellipse = T, groups = groupBY_H)


#################
# just for Bocagea
groupBY_B <- Boca_dat$species
pca_dat_B <- Boca_dat[,-1]
pca_dat_B <- pca_dat_B[,-1]

cor_B <- cor(pca_dat_B, method = 'pearson') # do a pearson correlation between all variables
corrplot(cor_B,
         order = 'hclust', # order by hierarchical clusting
         addrect = 10) # add rectangle around clusters identified by hclust


#colnames(pca) <- nms
pc_B <- prcomp(pca_dat_B, scale. = TRUE)
ggbiplot(pc_B, scale=TRUE, ellipse = TRUE, groups = groupBY_B) + ggtitle('PCA for Bocagea')
ggbiplot(pc_B, scale=TRUE, ellipse = TRUE, groups = groupBY_B, choices = c(3,4)) # pc axies 3 and 4
# as seen before and in other analyses, H.bryotrophe spans across many other species values...

#################
# just for Trigynaea
groupBY_T <- Trig_dat$species
pca_dat_T <- Trig_dat[,-1]
pca_dat_T <- pca_dat_T[,-1]

cor_T <- cor(pca_dat_T, method = 'pearson') # do a pearson correlation between all variables
corrplot(cor_T,
         order = 'hclust', # order by hierarchical clusting
         addrect = 10) # add rectangle around clusters identified by hclust




#colnames(pca) <- nms
pc_T <- prcomp(pca_dat_T, scale. = TRUE)
ggbiplot(pc_T, scale=FALSE, ellipse = TRUE, groups = groupBY_T) #+ ggtitle('PCA for Trigynaea')
ggbiplot(pc_T, scale=FALSE, ellipse = TRUE, groups = groupBY_T, choices = c(3,4)) # pc axies 3 and 4
# as seen before and in other analyses, H.bryotrophe spans across many other species values...
cur_T <- dCUR::CUR(pca_dat_T,
                   variables =  colnames(pca_dat),
                   cur_method = 'mixture',
                   rows = 0.999,
                   columns = 0.3,
                   standardize = T)

red_pca_dat_T <- cur_T$C_cur
rpc_T <- prcomp(red_pca_dat_T, scale. = F)
ggbiplot(rpc_T, scale = TRUE, ellipse = T, groups = groupBY_T)




##########################################
# whats the deal with hypervolumes?      #
##########################################

#just the climate variables we want

abio <- c('species','AF', 'Bio01', 'Bio07','soil_pH', 'Bio12', 'nitrogen', 'Bio14', 'vol_fraction_coarse', 'cation_exchange_capacity')
# Round 1: A niche best describing 
# hyp_dat <- na.omit(ab_dat[ab_dat$genus %in% c('Hornschuchia', 'Trigynaea', 'Bocagea'), 
#                            abio])

hyp_dat <- na.omit(clean_dat[clean_dat$genus %in% c('Hornschuchia', 'Trigynaea', 'Bocagea'), 
                          abio])
AF <- hyp_dat$AF
pc_hyp_dat <- hyp_dat[,-1]

pc_red_dim <- prcomp(pc_hyp_dat[,-1], scale = T)
ggbiplot(pc_red_dim, scale = T, ellipse = T, groups = AF) + ggtitle('reduced variables grouped by AF')
#library(hypervolume)


library(MASS)
da_dat <- hyp_dat[,-1]
af_vector <- da_dat[,'AF']
da_dat <- da_dat[,-1]

lda1 <- lda(da_dat, grouping = af_vector )
lda1
lda1pred <- predict(lda1)
ldahist(data = lda1pred$x[,1], g=af_vector)
ldahist(data = lda1pred$x[,2], g=af_vector)


#install.packages('dynRB')
library(dynRB)

r <- dynRB_VPa(hyp_dat[,-2])
overview(r)


om <- reshape(r$result[,c(1,2,4)], direction="wide", idvar="V1", timevar="V2") 
mantel(as.dist(om[2:ncol(om)]), as.dist(t(om[2:ncol(om)])), permutations = 1000)
plot(as.dist(om[2:ncol(om)]), as.dist(t(om[2:ncol(om)]))) 

result <- r$result
Overlap <- as.numeric(ifelse(result$V1 == result$V2, "NA", result$port_prod))  
is.numeric(Overlap)
Result2 <- cbind(result, Overlap)


breaks <- seq(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE),  
              by=round(max(Overlap, na.rm=TRUE)/10, digits=3))
col1 <- colorRampPalette(c("white", "navyblue")) #define color gradient


ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black") +
  scale_fill_gradientn(colours=col1(8), breaks=breaks, guide="colorbar",  
                       limits=c(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE))) +
  ggtitle('N-hypervolume overlap for reduced # of variables') +
  theme(
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(colour="black", size = rel(1.5), angle=35, hjust = 1),
    axis.text.y = element_text(colour="black", size = rel(1.5)),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

### do the same but over pc axes (inbuilt in dynRB)
r <- dynRB_VPa(hyp_dat, pca.corr = T)
#overview(r)


om <- reshape(r$result[,c(1,2,4)], direction="wide", idvar="V1", timevar="V2") 

mantel(as.dist(om[2:ncol(om)]), as.dist(t(om[2:ncol(om)])), permutations = 1000)
plot(as.dist(om[2:ncol(om)]), as.dist(t(om[2:ncol(om)]))) 


result <- r$result
Overlap <- as.numeric(ifelse(result$V1 == result$V2, "NA", result$port_prod))  
is.numeric(Overlap)
Result2 <- cbind(result, Overlap)
ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black")

breaks <- seq(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE),  
              by=round(max(Overlap, na.rm=TRUE)/10, digits=3))
col1 <- colorRampPalette(c("white", "navyblue")) #define color gradient


ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black") +
  scale_fill_gradientn(colours=col1(8), breaks=breaks, guide="colorbar",  
                       limits=c(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE))) +
  ggtitle('N-hypervolume overlap for \n reduced # of variables over  \nPC axes') +
  
  theme(
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(colour="black", size = rel(1.5), angle=35, hjust = 1),
    axis.text.y = element_text(colour="black", size = rel(1.5)),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )




overlap_per_dimension <- dynRB_Vn(hyp_dat, steps = 201)

cpca 




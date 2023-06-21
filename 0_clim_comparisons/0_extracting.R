####################################################################################################
#  2023-04-25 sjrs 
# HORNSCHUCHIA and the Atlantic Forest
# 
# extract climate for all occurrence records
####################################################################################################
rm(list = ls())

# INIT
wd <- '~/Sync/1_Annonaceae/H_Hornschuchia_AF/0_clim_comparisons/' # the working dir
dat_file <- '~/Sync/1_Annonaceae/H_Hornschuchia_AF/Y_data/Hornschuchia_Bocagea_Trigynaea_occurence.csv'
#dat_file <- '~/Sync/1_Annonaceae/8_SDM/1_data/Bocageeae_records_clean.csv'
pl <- c('ggplot2', 'raster', 'terra')
sapply(pl, require, character.only=TRUE)
source('~/Sync/serafins_functions.R')

setwd(wd)
dat <- read.csv(dat_file, header = T, sep =';')

# to subset it is easier if i have a genus and a species
dat$genus        <- sapply(dat$name, function(i) stringr::str_split(i, '_')[[1]][1])
dat$specific_epi <- sapply(dat$name, function(i) stringr::str_split(i, '_')[[1]][2])

# sanity check
print(unique(dat$genus))
# now we can filter out all the unwanted genera :-)
# dat <- dat[dat$genus == 'Hornschuchia',]

# now we have a nice reduced expert dataset


###########################
# load abiotic variables  #
###########################

abio_dir <- '~/Sync/1_Annonaceae/Y_DATA/1_abiotic_rasters/'

clim_stk <- raster::stack(paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO01.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO02.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO03.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO04.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO05.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO06.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO07.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO08.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO09.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO10.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO11.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO12.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO13.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO14.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO15.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO16.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO17.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO18.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO19.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO20.tif'),
                          paste0(abio_dir, '1_bioclim/MOD11C3v6.0-CHIRPSv2.0_BIOCLIMS_03m/BIO21.tif'))

soil_stk <- raster::stack(paste0(abio_dir, '2_soilgrids/bdod_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/cec_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/cfvo_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/clay_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/nitrogen_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/ocd_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/ocs_0-30cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/phh2o_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/sand_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/silt_5-15cm_mean_5000.tif'),
                          paste0(abio_dir, '2_soilgrids/soc_5-15cm_mean_5000.tif'))
clim_stk <- terra::rast(clim_stk)                      
clim_stk[clim_stk == -9999] <- NA # add NAs instead of -9999

soil_stk <- terra::rast(soil_stk)                
soil_stk <- terra::project(soil_stk, clim_stk) # project soil data onto climate extend and resolution        

env_stk <- c(clim_stk, soil_stk)
names(env_stk) <- c('Bio01', 'Bio02', 'Bio03', 'Bio04', 'Bio05', 'Bio06', 'Bio07',
                    'Bio08', 'Bio09', 'Bio10', 'Bio11', 'Bio12', 'Bio13', 'Bio14',
                    'Bio15', 'Bio16', 'Bio17', 'Bio18', 'Bio19', 'Bio20', 'Bio21',
                    'bulk_dens', 'cation_exchange_capacity', 'vol_fraction_coarse',
                    'clay_prop', 'nitrogen', 'org_C_dens', 'org_C_stock', 'soil_pH',
                    'sand_prop', 'silt_prop', 'soil_org_carbon')


#################################
# prepare coordinates, extract  #
#################################

# extract abitoic values for all layers and points
coords <- dat[, c("long", "lat")]

abio <- NULL
for(i in names(env_stk)){
  print(i)
  # loop thorugh abiotic layers
  lay <- env_stk[[i]]
  abio_tmp <- as.data.frame(terra::extract(lay, coords))
  colnames(abio_tmp[2]) <- i
  if(i == names(env_stk[[1]])){abio <- abio_tmp[2]}
  else{abio <- cbind(abio, abio_tmp[2])}
}



ab_dat <- cbind(dat, abio)
#save(ab_dat, file = 'new_data_extr.Rdata')
###################################



##############################################
# and additionally with soil maps from Brazil#
##############################################

s_map <- raster::stack('../Y_data/LC08_EOS_Maps_1155/data/Brazil_Soils_Map_19class.tif',
                       '../Y_data/LC08_EOS_Maps_1155/data/Brazil_Soils_Map_70class.tif',
                       '../Y_data/LC08_EOS_Maps_1155/data/Brazil_Soils_Map_249class.tif')
s_map <- terra::rast(s_map)   
names(s_map) <- c('soilmap19', 'soilmap70', 'soilmap249')


coords <- dat[, c("decimalLongitude", "decimalLatitude")]

abio <- NULL
for(i in names(s_map)){
  print(i)
  # loop thorugh abiotic layers
  lay <- s_map[[i]]
  abio_tmp <- as.data.frame(terra::extract(lay, coords))
  colnames(abio_tmp[2]) <- i
  if(i == names(s_map[[1]])){abio <- abio_tmp[2]}
  else{abio <- cbind(abio, abio_tmp[2])}
}

ab_dat <- cbind(ab_dat, abio)
save(ab_dat, file = 'new_data_extr.Rdata')



  
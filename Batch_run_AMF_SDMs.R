library(devtools)
install_github("bobmuscarella/ENMeval@Version-0.2.2")

library(ENMeval)
library(raster)
library(dismo)

############################################
#### PART 1: Set your working directory and load data
############################################

# setwd("/Users/au529793/Projects/Postdoc/AMF_distributions")
setwd("E:/Bob/AMF_SDMs")

### These are the raster stacks from the dropbox folder
env_all <- stack('DATA/SDMLayers/allenv')
env_climate <- stack('DATA/SDMLayers/climateenv')
env_resources <- stack('DATA/SDMLayers/resourcesenv')

### This is the most recent occurrence data
occ_recs <- as.data.frame(readxl::read_excel('DATA/Locations_2018_03_17.xlsx'), sheet=1)


############################################
#### PART 2: Build models for units of various taxonomic ranks
############################################

### The function below will help to do this easily ### 
taxloop <- function(occ_recs=NULL, eall=env_all, eclim=env_climate, eres=env_resources,
                    minN=20, rank="MAARJAM_ID", method='checkerboard2', algorithm='maxent.jar',
                    background="target", parallel=F, ncores=12){

  namelist <- sort(unique(occ_recs[,rank]))
  namelist <- namelist[!is.na(namelist)]
  
  mods_list <- list()
  for(i in seq_along(namelist)){
    
    foctax <- namelist[i]
    print(paste('now working on:', foctax, ':', i, 'of', length(namelist)), immediate. = T)
    
    fococc <- occ_recs[occ_recs[,rank] %in% foctax,]
    bgocc <- occ_recs[!occ_recs[,rank] %in% foctax,]

    fococcrows <- !is.na(rowSums(extract(eall, fococc[,c('Longitude','Latitude')])))
    occ <- fococc[fococcrows, c('Longitude','Latitude')]
    
    if(nrow(occ) >= minN){
      
      if(background == "target"){
        bg <- bgocc[!is.na(rowSums(extract(eall, bgocc[,c('Longitude','Latitude')]))),c('Longitude','Latitude')]
      } else {
        mask <- sum(!is.na(eall)) == nlayers(eall)
        values(mask)[values(mask)==F] <- NA
        bg <- randomPoints(mask, 10000)
      }

        mod_all <- ENMevaluate(occ, eall, bg.coords=bg, 
                               method=method, parallel=parallel, 
                               numCores=ncores, algorithm=algorithm, RM=1, fc='L')
        mod_climate <- ENMevaluate(occ, eclim, bg.coords=bg, 
                                   method=method, parallel=parallel, 
                                   numCores=ncores, algorithm=algorithm, RM=1, fc='L')
        mod_resources <- ENMevaluate(occ, eres, bg.coords=bg, 
                                     method=method, parallel=parallel, 
                                     numCores=ncores, algorithm=algorithm, RM=1, fc='L')
      }

      mods_list[[i]] <- list(all=mod_all, climate=mod_climate, resources=mod_resources)
      names(mods_list)[i] <- foctax
  }
  mods_list <- Filter(length, mods_list)
  return(mods_list)
}


genus_models <- taxloop(occ_recs=fiocc, rank='Genus')
family_models <- taxloop(occ_recs=fiocc, rank='Family')
order_models <- taxloop(occ_recs=fiocc, rank='Order')




### TAXA MODELS (SAVE ONE AT A TIME)

namelist <- sort(unique(occ_recs[,rank]))
namelist <- namelist[!is.na(namelist)]

for(i in seq_along(namelist)){
  
  foctax <- namelist[i]
  print(paste('now working on:', foctax, ':', i, 'of', length(namelist)), immediate. = T)
  
  fococc <- occ_recs[occ_recs[,rank] %in% foctax,]
  bgocc <- occ_recs[!occ_recs[,rank] %in% foctax,]

  fococcrows <- !is.na(rowSums(extract(eall, fococc[,c('Longitude','Latitude')])))
  occ <- fococc[fococcrows, c('Longitude','Latitude')]
  
  if(nrow(occ) >= minN){
    
    if(background == "target"){
      bg <- bgocc[!is.na(rowSums(extract(eall, bgocc[,c('Longitude','Latitude')]))),c('Longitude','Latitude')]
    } else {
      mask <- sum(!is.na(eall)) == nlayers(eall)
      values(mask)[values(mask)==F] <- NA
      bg <- randomPoints(mask, 10000)
    }
    
    mod_all <- ENMevaluate(occ, eall, bg.coords=bg, 
                           method=method, parallel=parallel, 
                           numCores=ncores, algorithm=algorithm)
    mod_climate <- ENMevaluate(occ, eclim, bg.coords=bg, 
                               method=method, parallel=parallel, 
                               numCores=ncores, algorithm=algorithm)
    mod_resources <- ENMevaluate(occ, eres, bg.coords=bg, 
                                 method=method, parallel=parallel, 
                                 numCores=ncores, algorithm=algorithm)
  }
  
  mods <- list()
  mods[[1]] <- list(all=mod_all, climate=mod_climate, resources=mod_resources)
  names(mods) <- foctax
  
  filename <- paste0('RESULTS/', foctax, '.RDS')
  saveRDS(mods, file=filename)
  
}



########################
### Part 3. Explore the output (we can work on developing this as necessary)
########################
# e.g.,
order_models$Paraglomerales$all@results

# Predict the model across you study extent and visualize it
plot(order_models$Paraglomerales$all@predictions[[1]])
points(order_models$Paraglomerales$all@occ.pts, pch=16, cex=.4)
# (Note that there are lots of gaps in the map because there are gaps in the soil P map)

# Look at response curves
response(order_models$Paraglomerales$all@models[[1]])


library(devtools)
install_github("bobmuscarella/ENMeval@Version-0.2.2")

library(ENMeval)
library(raster)
library(dismo)

############################################
#### PART 1: Set your working directory and load data
############################################

setwd("/Users/au529793/Projects/Postdoc/AMF_distributions")

### These are the raster stacks from the dropbox folder
env_all <- stack('DATA/SDMLayers/allenv')
env_climate <- stack('DATA/SDMLayers/climateenv')
env_resources <- stack('DATA/SDMLayers/resourcesenv')

### This is the most recent occurrence data
occ_recs <- as.data.frame(readxl::read_excel('DATA/Locations_2018_03_17.xlsx'), sheet=1)


############################################
#### PART 2: Build models for units of various taxonomic ranks
############################################

### The function below will help to do this easily
### 
taxloop <- function(occ_recs=NULL, envall=env_all, envclim=env_climate, envres=env_resources,
                    minN=20, rank="MAARJAM_ID", method='checkerboard2', algorithm='maxent.jar'){
  
  namelist <- sort(unique(occ_recs[,rank]))
  namelist <- namelist[!is.na(namelist)]
  
  mods_list <- list()
  for(i in seq_along(namelist)){
    
    foctax <- namelist[i]
    fococc <- occ_recs[occ_recs[,rank] %in% foctax,]
    bgocc <- occ_recs[!occ_recs[,rank] %in% foctax,]
    
    if(nrow(fococc) >= minN){
      
      fococcrows <- !is.na(rowSums(extract(envall, fococc[,c('Longitude','Latitude')])))
      occ <- fococc[fococcrows,c('Longitude','Latitude')]
      
      ### The following line selects background as all 'non-focal' occurrences
      ### (a problem here is that there will be few background points for some models [e.g., the Glomerales order model])
      focbgrows <- !is.na(rowSums(extract(envall, bgocc[,c('Longitude','Latitude')])))
      bg <- bgocc[focbgrows,c('Longitude','Latitude')]
      
      mod_all <- ENMevaluate(occ, envall, bg.coords=bg, method=method, parallel=T, algorithm=algorithm)
      mod_climate <- ENMevaluate(occ, envclim, bg.coords=bg, method=method, parallel=T, algorithm=algorithm)
      mod_resources <- ENMevaluate(occ, envres, bg.coords=bg, method=method, parallel=T, algorithm=algorithm)
      
      mods_list[[i]] <- list(all=mod_all, climate=mod_climate, resources=mod_resources)
      names(mods_list)[i] <- foctax
    }
  }
  return(mods_list)
}


### Now use the function to make models of each taxon in each rank
taxa_models <- taxloop(occ_recs=occ_recs, rank='MAARJAM_ID')
genus_models <- taxloop(occ_recs=occ_recs, rank='Genus')
family_models <- taxloop(occ_recs=occ_recs, rank='Family')
order_models <- taxloop(occ_recs=occ_recs, rank='Order')


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

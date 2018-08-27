############################################
#### AUTOMATED RUNS OF AM FUNGI SDMS
############################################

### Load an older version of ENMeval - this should be changed to work with the newest version.
library(devtools)
install_github("bobmuscarella/ENMeval@Version-0.2.2")

library(ENMeval)
library(raster)
library(dismo)

############################################
#### PART 1: Set the working directory and load data
############################################

setwd("E:/Bob/AMF_SDMs")

### Load raster stacks (from the dropbox folder)
env_all <- stack('DATA/SDMLayers/allenv')
env_climate <- stack('DATA/SDMLayers/climateenv')
env_resources <- stack('DATA/SDMLayers/resourcesenv')

### Load the occurrence data
fiocc <- as.data.frame(readxl::read_excel('DATA/Locations_2018_03_17.xlsx'), sheet=1)

### Load helper functions
# The 'taxloop' function is to automate running models 
taxloop <- function(occ_recs=NULL, eall=env_all, eclim=env_climate, eres=env_resources,
                    minN=20, rank="MAARJAM_ID", method='checkerboard2',
                    background="random", parallel=T, ncores=12, ...){
  
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
                             numCores=ncores, ...)
      mod_climate <- ENMevaluate(occ, eclim, bg.coords=bg, 
                                 method=method, parallel=parallel, 
                                 numCores=ncores, ...)
      mod_resources <- ENMevaluate(occ, eres, bg.coords=bg, 
                                   method=method, parallel=parallel, 
                                   numCores=ncores, ...)
    }
    
    mods_list[[i]] <- list(all=mod_all, climate=mod_climate, resources=mod_resources)
    names(mods_list)[i] <- foctax
  }
  mods_list <- Filter(length, mods_list)
  return(mods_list)
}

# The 'summarize.aic.mods' function summarizes results of models after selecting the top model with AIC
summarize.aic.mods <- function(mod, allvarnames, rank){
  taxa <- rep(names(mod), each=3)
  model.name <- rep(c('all','climate','resources'), times=length(mod))
  minAIClist <- lapply(mod, function(x) {
    lapply(x, function(x){
      unique(x@results$AICc[x@results$delta.AICc==0 & !is.na(x@results$delta.AICc)])
    })
  })
  
  featureslist <- lapply(mod, function(x) {
    lapply(x, function(x){
      as.character(unique(x@results$features[x@results$delta.AICc==0 & !is.na(x@results$delta.AICc)]))[1]
    })
  })
  
  RMlist <- lapply(mod, function(x) {
    lapply(x, function(x){
      max(unique(x@results$rm[x@results$delta.AICc==0 & !is.na(x@results$delta.AICc)]))
    })
  })
  
  df <- data.frame(taxa=taxa, rank=rank, model.name=model.name, 
                   min.AICc=unlist(minAIClist), features=unlist(featureslist), 
                   rm=unlist(RMlist))
  
  ORminlist <- lapply(mod, function(x) {
    lapply(x, function(x){
      max(unique(x@results$Mean.ORmin[x@results$delta.AICc==0 & !is.na(x@results$delta.AICc)]))
    })
  })
  
  OR10list <- lapply(mod, function(x) {
    lapply(x, function(x){
      max(unique(x@results$Mean.OR10[x@results$delta.AICc==0 & !is.na(x@results$delta.AICc)]))
    })
  })
  
  AUClist <- lapply(mod, function(x) {
    lapply(x, function(x){
      max(unique(x@results$full.AUC[x@results$delta.AICc==0 & !is.na(x@results$delta.AICc)]))
    })
  })
  
  df$ormin <- unlist(ORminlist)
  df$or10 <- unlist(OR10list)
  df$auc <- unlist(AUClist)
  df$index <- match(paste0(df$features, '_', df$rm), as.character(mod[[1]]$all@results$settings))
  
  
  tmpdf <- data.frame()
  for(i in 1:length(unlist(mod))){
    tmptab <- unlist(mod)[[i]]@models[[df$index[i]]]@results
    var.names <- rownames(tmptab)[grepl('permutation',rownames(tmptab))]
    var.names <- gsub('.permutation.importance','',var.names)
    pi.vals <- setNames(tmptab[grepl('permutation',rownames(tmptab))], var.names)
    pi.vals <- pi.vals[match(allvarnames, names(pi.vals))]
    vi.vals <- setNames(tmptab[grepl('.contribution',rownames(tmptab))], var.names)
    vi.vals <- vi.vals[match(allvarnames, names(vi.vals))]
    tmpdf <- rbind(tmpdf,
                   setNames(c(pi.vals, vi.vals), c(paste0('pi.',allvarnames), paste0('vi.',allvarnames))))
    names(tmpdf) <-  c(paste0('pi.',allvarnames), paste0('vi.',allvarnames))
  }
  
  df <- cbind(df, tmpdf)
  
  rownames(df) <- NULL
  return(df)
}


############################################
#### PART 2: Build models for units of various taxonomic ranks
############################################

############
### TARGET-GROUP BACKGROUND
############

### Run order, family and genus models at once and save the results together
order_models <- taxloop(occ_recs=fiocc, rank='Order', background='target')
saveRDS(order_models, file="RESULTS/_order_targetBG_20180719.RDS")

family_models <- taxloop(occ_recs=fiocc, rank='Family', background='target')
saveRDS(family_models, file="RESULTS/_family_targetBG_20180719.RDS")

genus_models <- taxloop(occ_recs=fiocc, rank='Genus', background='target')
saveRDS(genus_models, file="RESULTS/_genus_targetBG_20180719.RDS")

### Run taxon-level models and save one at a time
background <- 'target'
minN <- 20
namelist <- sort(unique(fiocc[,'MAARJAM_ID']))
namelist <- namelist[!is.na(namelist)]

for(i in seq_along(namelist)){
  
  foctax <- namelist[i]
  print(paste('now working on:', foctax, ':', i, 'of', length(namelist)), immediate. = T)
  
  fococc <- fiocc[fiocc[,'MAARJAM_ID'] %in% foctax,]
  bgocc <- fiocc[!fiocc[,'MAARJAM_ID'] %in% foctax,]
  
  fococcrows <- !is.na(rowSums(extract(env_all, fococc[,c('Longitude','Latitude')])))
  occ <- fococc[fococcrows, c('Longitude','Latitude')]
  
  if(nrow(occ) >= minN){
    
    if(background == "target"){
      bg <- bgocc[!is.na(rowSums(extract(env_all, bgocc[,c('Longitude','Latitude')]))),c('Longitude','Latitude')]
    } else {
      mask <- sum(!is.na(env_all)) == nlayers(env_all)
      values(mask)[values(mask)==F] <- NA
      bg <- randomPoints(mask, 10000)
    }
    
    mod_all <- ENMevaluate(occ, env_all, bg.coords=bg, 
                           method='checkerboard2', parallel=T, 
                           numCores=12)
    mod_climate <- ENMevaluate(occ, env_climate, bg.coords=bg, 
                               method='checkerboard2', parallel=T, 
                               numCores=12)
    mod_resources <- ENMevaluate(occ, env_resources, bg.coords=bg, 
                                 method='checkerboard2', parallel=T, 
                                 numCores=12)
  }
  
  mods <- list()
  mods[[1]] <- list(all=mod_all, climate=mod_climate, resources=mod_resources)
  names(mods) <- foctax
  
  filename <- paste0('RESULTS/', foctax, '_targetBG_20180719.RDS')
  saveRDS(mods, file=filename)
  
}


############
### RANDOM BACKGROUND
############

### Run order, family and genus models at once and save the results together
order_models <- taxloop(occ_recs=fiocc, rank='Order')
saveRDS(order_models, file="RESULTS/_order_randomBG_20180719.RDS")

family_models <- taxloop(occ_recs=fiocc, rank='Family')
saveRDS(family_models, file="RESULTS/_family_randomBG_20180719.RDS")

genus_models <- taxloop(occ_recs=fiocc, rank='Genus')
saveRDS(genus_models, file="RESULTS/_genus_randomBG_20180719.RDS")

### Run taxon-level models and save one at a time
background <- 'random'
minN <- 20
namelist <- sort(unique(fiocc[,'MAARJAM_ID']))
namelist <- namelist[!is.na(namelist)]

for(i in seq_along(namelist)){
  
  foctax <- namelist[i]
  print(paste('now working on:', foctax, ':', i, 'of', length(namelist)), immediate. = T)
  
  fococc <- fiocc[fiocc[,'MAARJAM_ID'] %in% foctax,]
  bgocc <- fiocc[!fiocc[,'MAARJAM_ID'] %in% foctax,]
  
  fococcrows <- !is.na(rowSums(extract(env_all, fococc[,c('Longitude','Latitude')])))
  occ <- fococc[fococcrows, c('Longitude','Latitude')]
  
  if(nrow(occ) >= minN){
    
    if(background == "target"){
      bg <- bgocc[!is.na(rowSums(extract(env_all, bgocc[,c('Longitude','Latitude')]))),c('Longitude','Latitude')]
    } else {
      mask <- sum(!is.na(env_all)) == nlayers(env_all)
      values(mask)[values(mask)==F] <- NA
      bg <- randomPoints(mask, 10000)
    }
    
    mod_all <- ENMevaluate(occ, env_all, bg.coords=bg, 
                           method='checkerboard2', parallel=T, 
                           numCores=12)
    mod_climate <- ENMevaluate(occ, env_climate, bg.coords=bg, 
                               method='checkerboard2', parallel=T, 
                               numCores=12)
    mod_resources <- ENMevaluate(occ, env_resources, bg.coords=bg, 
                                 method='checkerboard2', parallel=T, 
                                 numCores=12)
  }
  
  mods <- list()
  mods[[1]] <- list(all=mod_all, climate=mod_climate, resources=mod_resources)
  names(mods) <- foctax
  
  filename <- paste0('RESULTS/', foctax, '_randomBG_20180719.RDS')
  saveRDS(mods, file=filename)
  
}



############################################
#### PART 3: SUMMARIZE AIC SELECTED TARGET BACKGROUND MODELS ###
############################################
list.files('RESULTS/20180719_targetBG')

### Summarize order, family, and genus models
mod <- readRDS('RESULTS/20180719_targetBG/_order_targetBG_20180719.RDS')
order.df <- summarize.aic.mods(mod, names(env_all), rank='order')

mod <- readRDS('RESULTS/20180719_targetBG/_family_targetBG_20180719.RDS')
family.df <- summarize.aic.mods(mod, names(env_all), rank='family')

mod <- readRDS('RESULTS/20180719_targetBG/_genus_targetBG_20180719.RDS')
genus.df <- summarize.aic.mods(mod, names(env_all), rank='genus')

### Summarize OTU models
otus <- list.files('RESULTS/20180719_targetBG')[grepl('VTX',list.files('RESULTS/20180719_targetBG'))]

otu.df <- data.frame()
for(i in seq_along(otus)){
  print(paste(i, 'out of', length(otus)))
  mod <- readRDS(paste0('RESULTS/20180719_targetBG/', otus[i]))
  otu.df <- rbind(otu.df, summarize.aic.mods(mod, names(env_all), rank='otu'))
  rm(mod)
}

df1 <- rbind(order.df, family.df, genus.df, otu.df)
df1$bg.method <- 'target'

### Save target background results
saveRDS(df1, file='RESULTS/20180719_targetBG_summary_table.RDS')


#######################################################
### SUMMARIZE AIC SELECTED RANDOM BACKGROUND MODELS ###
#######################################################
list.files('RESULTS/20180719_randomBG')

### Summarize order, family, and genus models
order_mod <- readRDS('RESULTS/20180719_randomBG/_order_randomBG_20180719.RDS')
order.df <- summarize.aic.mods(order_mod, names(env_all), rank='order')

family_mod <- readRDS('RESULTS/20180719_randomBG/_family_randomBG_20180719.RDS')
family.df <- summarize.aic.mods(family_mod, names(env_all), rank='family')

genus_mod <- readRDS('RESULTS/20180719_randomBG/_genus_randomBG_20180719.RDS')
genus.df <- summarize.aic.mods(genus_mod, names(env_all), rank='genus')

### Summarize OTU models
otus <- list.files('RESULTS/20180719_randomBG')[grepl('VTX',list.files('RESULTS/20180719_randomBG'))]

otu.df <- data.frame()
for(i in seq_along(otus)){
  print(paste(i, 'out of', length(otus)))
  otu_mod <- readRDS(paste0('RESULTS/20180719_randomBG/',otus[i]))
  otu.df <- rbind(otu.df, summarize.aic.mods(otu_mod, names(env_all), rank='otu'))
  rm(otu_mod)
}

df2 <- rbind(order.df, family.df, genus.df, otu.df)
df2$bg.method <- 'random'

### Save random background results
saveRDS(df2, file='RESULTS/20180719_randomBG_summary_table.RDS')

### Compile results from random and target background
df1 <- readRDS('RESULTS/20180719_targetBG_summary_table.RDS')
df2 <- readRDS('RESULTS/20180719_randomBG_summary_table.RDS')

### Write compiled data as .RDS and .csv files
df <- rbind(df1, df2)
write.csv(df, 'RESULTS/20180810_AIC_selected_summary_table.csv')

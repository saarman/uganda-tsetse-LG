library(raster)
library(rgdal)

#function to create custom biovars for wettest/warmest season. Format: mean [var1] of max[var2] season (ex: mean temp of wettest season)
bio_custom_max <- function(var1,var2){
  #average var2 values by season (Uganda Seasons: WS1:April-May,DS1:June-July,WS2:August-September,DS2:October-March)
  season.mean2 <- stackApply(var2, indices=c(4,4,4,1,1,2,2,3,3,4,4,4), fun=mean) 
  #extract and sum all of the var2 values for each pixel in each layer
  r.sum <- cellStats(season.mean2,'sum')
  #identify the layer with the max var2 value
  r.max <- which.max(r.sum)
  #take the mean var1 by season for the final output
  season.mean1 <- stackApply(var1, indices=c(4,4,4,1,1,2,2,3,3,4,4,4), fun=mean)
  #subset the var1 data set to only include the layer with the corresponding max var2 value 
  stack(subset(season.mean1,as.numeric(r.max)))
}

#function to create custom biovars for driest/coldest. Format: mean [var1] of min[var2] season (ex: mean temp of driest season)
bio_custom_min <- function(var1,var2){
  #average var2 values by season
  season.mean2 <- stackApply(var2, indices=c(4,4,4,1,1,2,2,3,3,4,4,4), fun=mean)
  #extract and sum all of the var2 values for each pixel in each layer
  r.sum <- cellStats(season.mean2,'sum')
  #identify the layer with the min var2 value
  r.min <- which.min(r.sum)
  #take the mean var1 by season for the final output
  season.mean1 <- stackApply(var1, indices=c(4,4,4,1,1,2,2,3,3,4,4,4), fun=mean)
  #subset the var1 data set to only include the layer with the corresponding min var2 value 
  stack(subset(season.mean1,as.numeric(r.min)))
}

#function to create seasonal bioclim variables
#temp and prec objects should be RasterStacks with 12 layers (months)
biovars_custom <- function(temp,prec){
  #BIO8S = Mean Temperature of Wettest Season
  BIO8S<-bio_custom_max(temp,prec)
  #BIO9S = Mean Temperature of Driest Season
  BIO9S<-bio_custom_min(temp,prec)
  #BIO10S = Mean Temperature of Warmest Season
  BIO10S<-bio_custom_max(temp,temp)
  #BIO11S = Mean Temperature of Coldest Season
  BIO11S<-bio_custom_min(temp,temp)
  #BIO16S = Precipitation of Wettest Season
  BIO16S<-bio_custom_max(prec,prec)
  #BIO17S = Precipitation of Driest Season
  BIO17S<-bio_custom_min(prec,prec)
  #BIO18S = Precipitation of Warmest Season
  BIO18S<-bio_custom_max(prec,temp)
  #BIO19S = Precipitation of Coldest Season
  BIO19S<-bio_custom_min(prec,temp) 
  stack(BIO8S,BIO9S,BIO10S,BIO11S,BIO16S,BIO17S,BIO18S,BIO19S)
}

#for loop to generate seasonal bioclim variables by year
for (year in c(2008:2013)){
  prec=stack(paste0("/home/fas/caccone/apb56/project/CHELSA/prec/uganda_kenya_clips/CHELSA_prec_",year,"_",c("01","02","03","04","05","06","07","08","09","10","11","12"),"_V1.2.1_UgandaKenyaClip.tif"))
  temp=stack(paste0("/home/fas/caccone/apb56/project/CHELSA/tmean/uganda_kenya_clips/CHELSA_tmean_",year,"_",c("01","02","03","04","05","06","07","08","09","10","11","12"),"_V1.2.1_UgandaKenyaClip.tif"))
  #make sure extents are the same (Uganda/Kenya)
  prec<-setExtent(prec,extent(28.6,42.5,-4.8,5))
  temp<-setExtent(temp,extent(28.6,42.5,-4.8,5))
  #create custom biovars
  biovars<-biovars_custom(temp,prec)
  writeRaster(biovars,paste0("/home/fas/caccone/apb56/project/CHELSA/biovars/UgandaBiovarsSeasonal",year),format="GTiff",overwrite=TRUE)
}


#Average seasonal bioclim variables for all years (current)
biostack<-stack(paste0("/home/fas/caccone/apb56/project/CHELSA/biovars/UgandaBiovarsSeasonal",c(2008:2013),".tif"))
indices<-rep(c(1:8), times = 6) #combine corresponding biovar layers by mean (take mean of all biovars for all years)
biostack_mean<-stackApply(biostack, indices, fun=mean)
names(biostack_mean)<-c("BIO8S","BIO9S","BIO10S","BIO11S","BIO16S","BIO17S","BIO18S","BIO19S")
writeRaster(biostack_mean,paste0("/home/fas/caccone/apb56/project/CHELSA/biovars/UgandaBiovarsSeasonalAllYears"),format="GTiff",overwrite=TRUE)

#Average seasonal bioclim variables for all years (future)
prec.f=stack(paste0("/home/fas/caccone/apb56/project/CHELSA/future/prec/uganda_clips/CHELSA_prec_",c("1","2","3","4","5","6","7","8","9","10","11","12"),"_rcp45_2041-2060_UgandaClip.tif"))
temp.f=stack(paste0("/home/fas/caccone/apb56/project/CHELSA/future/tmean/uganda_clips/CHELSA_tmean_",c("1","2","3","4","5","6","7","8","9","10","11","12"),"_rcp45_2041-2060_UgandaClip.tif"))
#convert temp to same units as current
temp.f<-temp.f*10
#create custom biovars
biovars.f<-biovars_custom(temp.f,prec.f)
names(biovars.f)<-c("BIO8S","BIO9S","BIO10S","BIO11S","BIO16S","BIO17S","BIO18S","BIO19S")
writeRaster(biovars.f,paste0("/home/fas/caccone/apb56/project/CHELSA/biovars/UgandaBiovarsSeasonalFuture"),format="GTiff",overwrite=TRUE)
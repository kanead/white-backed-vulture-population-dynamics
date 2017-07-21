rm(list=ls(all=TRUE)) 

library(maptools)
library(SDMTools)
library(raster)
library(sp)
library(dplyr)
setwd("C:\\Users\\akane\\Desktop\\Science\\Manuscripts\\White-backed Vulture Pop Dynamics\\Code\\white-backed-vulture-population-dynamics")

mydata<-read.csv("resights.csv",header = T,sep = ",")
#mydata<-head(mydata)


crswgs84=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
Kruger=readShapePoly("WDPA_July2017_protected_area_873-shapefile-polygons.shp",proj4string=crswgs84,verbose=TRUE)

#plot(Kruger)
#KrugerPoints <- mydata[mydata$Location=="Kruger" , ] 
#KrugerPoints<-droplevels(KrugerPoints)
#KZNPoints <- mydata[mydata$Location=="KZN" , ] 
#KZNPoints<-droplevels(KZNPoints)
#points(KrugerPoints$Latitude ~ KrugerPoints$Longitude, col = "red", cex = 1)
#points(KZNPoints$Latitude ~ KZNPoints$Longitude, col = "blue", cex = 1)

# One useful function that we will see in this post is the gContains function from the rgeos package. 
# This function checks whether a polygon contains a point or more generally whether a geometry contains another geometry.

mydata <- mydata[mydata$Location=="KZN" , ] 
mydata<-droplevels(mydata)


xy <- mydata[,c(3,2)]
xy<-droplevels(xy)

spdf <- SpatialPointsDataFrame(coords = xy, data = xy,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))


proj4string(spdf)
# [1] " +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
proj4string(Kruger)
# [1] " +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
proj4string(spdf) == proj4string(Kruger)

plot(Kruger)
plot(spdf, col="red" , add=TRUE)
res <- over(spdf, Kruger)
summary(res$NAME)

####################################################################################

#mydata<- mydata %>%
#  group_by(Location) %>% 
  
#  mydata[,c(3,2)]  %>%
#  droplevels(xy)

#spdf <- SpatialPointsDataFrame(coords = xy, data = xy,
         #                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))



#proj4string(spdf)
# [1] " +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
#proj4string(Kruger)
# [1] " +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
#proj4string(spdf) == proj4string(Kruger)

#plot(Kruger)
#plot(spdf, col="red" , add=TRUE)
#res <- over(spdf, Kruger)
#summary(res$NAME)




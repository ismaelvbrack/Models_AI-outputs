
library(sf)
library(terra)
library(ggplot2)

root <- file.path("Whitails_DeLuca")

# Import camera deployments ----------------------------------
deploys <- read.csv(file.path(root,"data","deployments.csv"))

# Transform to date
deploys$start_date <- as.Date(deploys$start_date, "%m/%d/%Y")
deploys$end_date <- as.Date(deploys$end_date, "%m/%d/%Y")

# Total time interval for each camera deployment
deploys$time.int <- as.numeric(deploys$end_date - deploys$start_date)

#* Define each camera deployment as a site
cameras <- unique(deploys[,c("placename","deployment_id","start_date","end_date","time.int",
                             "longitude","latitude")])

cameras <- cameras[order(cameras$placename),] # order

# Filter sites that ended before July/2022
#cameras <- cameras[which(as.Date(cameras$end_date) < "2022-07-01"),] 
#cameras <- cameras[-which(cameras$deployment_id == "CT-6B_Spring1"),]
#* Exclude short-time camera deployments
cameras <- cameras[which(cameras$time.int>20),]
cameras <- cameras[which(cameras$end_date < "2022-07-01"),] 
cameras <- cameras[-which(cameras$placename == "CT-Preserve1" |
                            cameras$placename == "CT-Preserve2"),]

# Problematic deployment (no images folder
cameras <- cameras[-which(cameras$deployment_id == "CT-1A_April28-May19"),]


# Define sites from camera deployments ------------------------------------
sites <- data.frame(
  site=unique(cameras$placename),
  deployment_id=NA,
  start_date=NA,
  end_date=NA
)
# Reference data to filter days
ref.date <- as.Date("2022-05-03")

# Number of days
ndays <- 15

# Split into a list of DF per site
cam.ls <- split(cameras, cameras$placename)

# for each site...
sites <- lapply(cam.ls, function(y){
  # Get the line position (deployment)
  y[which(apply(y, 1,
                # That the start date is BEFORE ref.date
                function(x) x["start_date"]<ref.date & 
                  # end date is AFTER ref.date
                  x["end_date"]>ref.date &
                  # there are at least 20 days from ref.date and end date
                  (ref.date+ndays) < x["end_date"])),]
})

for(i in 1:length(sites)){
  # If there is one deployment that fullfill the previous conditions:
  if(nrow(sites[[i]])==1){
    # Start date is the ref.date
    sites[[i]][,"start_date"] <- ref.date
    # Sum 9 days to define the end date
    sites[[i]][,"end_date"] <- ref.date+(ndays-1)
  }
  # If no deployment has the conditions
  if(nrow(sites[[i]])==0){
    # Get the maximum end date
    sites[[i]] <- cam.ls[[i]][which.max(cam.ls[[i]]$end_date),]
    # Subtract 10 to define start date
    sites[[i]][,"start_date"] <- sites[[i]]$end_date - (ndays)
    # Exclude last day 
    sites[[i]][,"end_date"] <- sites[[i]]$end_date - 1
  }
}

# Combine
sites <- do.call(rbind, sites)


# Map layers ----------------------------------------------------------------------------------------

sites <- st_as_sf(sites, coords=c("longitude","latitude"), crs=4326)

#st_write(sites,file.path(root,"GIS","DeLuca_27camtrap sites.shp"))

sites500 <- st_buffer(sites, 500)

lulc <- st_read(file.path(root,"GIS","DeLuca_Land_Cover.shp"))

forest <- rasterize(vect(lulc[which(lulc$LEVEL1_L_1 == "Upland Forest"),]), rast(vect(lulc),res=9.259259e-05) )

hardforest <- rasterize(vect(lulc[which(lulc$LEVEL2_L_1 == "Upland Hardwood Forests" | lulc$LEVEL2_L_1 == "Upland Mixed Forests"),]),
                    rast(vect(lulc),res=9.259259e-05) )

# open.classes <- c("Agriculture","Rangeland","Vegetated Non-Forested Wetlands")

open <- rasterize(vect(lulc[which(lulc$LEVEL1_L_1 == "Agriculture" | lulc$LEVEL1_L_1 == "Rangeland" |
                                    lulc$LEVEL2_L_1 == "Vegetated Non-Forested Wetlands"),]),
                            rast(vect(lulc),res=9.259259e-05) )

wetland <- rasterize(vect(lulc[which(lulc$LEVEL1_L_1 == "Wetlands"),]), rast(vect(lulc),res=9.259259e-05) )

unique(lulc$LEVEL1_L_1)

twi <- rast(file.path(root,"GIS","TWI_DeLuca.tif"))

plot(twi)
plot(st_geometry(lulc[which(lulc$LEVEL1_L_1 == "Upland Forest"),]),
     border="green",add=T)
plot(st_geometry(lulc[which(lulc$LEVEL1_L_1 == "Wetlands"),]),
     border="cyan",add=T)
plot(sites$geometry,,pch=19,col="red",add=T)
plot(sites500$geometry, border="orange", add=T)

hist(siteCovs$twi)

# Extract LULC and TWI --------------------------------------------------------------------------------------------

sites$twi <- extract(twi, sites500, mean)[,2]

sites$forest <- propMap(forest, vect(sites500))

sites$hardforest <- propMap(hardforest, vect(sites500))

sites$wetland <- propMap(wetland, vect(sites500))

sites$open <- propMap(open, vect(sites500))

pairs(cbind(rowSums(Y),st_drop_geometry(sites[c("twi","wetland","forest","open")])))

dets <- Y
Y <- ifelse(Y>0,1,0)

occu <- ifelse(rowSums(Y)>0,1,0)

tab <- cbind(st_drop_geometry(sites),occu=occu)

ggplot(data=tab, aes(x=twi,y=occu)) +
  geom_point() +
  geom_smooth(method = glm, method.args= list(family="binomial"))

ggplot(data=tab, aes(x=forest,y=occu)) +
  geom_point() +
  geom_smooth(method = glm, method.args= list(family="binomial"))

# Temporal covars -------------------------------------------------------------------------------------------------

weather <- read.csv(file.path(root,"data","Okeechobee_455_2022.csv"))
  
weather$Date.Time <- as.POSIXct(weather$Date.Time, format="%Y-%m-%d %H:%M")

weather <- weather[-which(is.na(weather$Date.Time)),]

rain <- temps <- matrix(NA, nrow=nrow(sites),ncol=ndays)

for(s in 1:nrow(sites)){
  
  mydates <- seq(sites$start_date[s], sites$end_date[s], by=1)
  for(j in 1:ndays){
    aqui <- which(as.Date(weather$Date.Time) %in% mydates[j])
    
    temps[s,j] <- mean(weather[aqui, "Temp...60cm..C."])
    rain[s,j] <- sum(weather[aqui, "Rainfall.Amount..in."])
  }
}

# Export ----------------------------------------------------------------------------------------------------------

write.table(st_drop_geometry(sites), file.path(root,"data","siteCovars.txt"))

write.table(temps, file.path(root,"data","Temperatures60cm.txt"))
write.table(rain, file.path(root,"data","Rainfall.txt"))


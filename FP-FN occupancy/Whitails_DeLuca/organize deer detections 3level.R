
library(ggplot2)

root <- file.path("Whitails_DeLuca")

# Reference data to filter days
ref.date <- as.Date("2022-05-03")

# Import camera deployments ----------------------------------
deploys <- read.csv(file.path(root,"data","deployments.csv"))

# Transform to date
deploys$start_date <- as.Date(deploys$start_date, "%m/%d/%Y")
deploys$end_date <- as.Date(deploys$end_date, "%m/%d/%Y")

# Total time interval for each camera deployment
deploys$time.int <- as.numeric(deploys$end_date - deploys$start_date)

#* Define each camera deployment as a site
cameras <- unique(deploys[,c("placename","deployment_id","start_date","end_date","time.int")])

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

# See time frame for each camera site
ggplot(data=cameras,
       aes(x=start_date,xend=end_date,y=placename,yend=placename,color=placename)) +
  geom_segment(linewidth=2) +
  theme_bw() +
  theme(legend.position="none")

# Define sites from camera deployments ------------------------------------
sites <- data.frame(
  site=unique(cameras$placename),
  deployment_id=NA,
  start_date=NA,
  end_date=NA
)

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

ggplot(data=sites,
       aes(x=start_date,xend=end_date,y=placename,yend=placename,color=placename)) +
  geom_segment(linewidth=2) +
  theme_bw() +
  theme(legend.position="none")

range(c(sites$start_date, sites$end_date))

# Import images information -----------------------------------------------
images <- read.csv(file.path(root,"data","images_2004061.csv"))

images$timestamp <- as.POSIXct(images$timestamp)

images <- images[order(images$deployment_id, images$timestamp),]

# any deplyment ID in images do not coincide with camera deployments?
any(!images$deployment_id %in% deploys$deployment_id)

# Filter images for the selected deployments
images <- images[which(images$deployment_id %in% sites$deployment_id),]

#* Select images in for the selected time frame for each deployment 
aqui <- which(apply(images, 1, # each line of images:
                    function(x) {
                      # match the deployment_id and get true or false for images AFTER the defined start_date
                      x["timestamp"] >= sites[match(x["deployment_id"],sites$deployment_id),"start_date"] &
                        # AND
                        # match the deployment_id and get true or false for images BEFORE the defined end_date
                        x["timestamp"] <= sites[match(x["deployment_id"],sites$deployment_id),"end_date"] 
                    }) # the 'which' will get the positions
)
images <- images[aqui,]

# Include camera site name
images$site <- sites[match(images$deployment_id, sites$deployment_id),"placename"]

# Group 3-image sequences into a single ID
images$imseq.id <- NA


images$imseq.id[1] <- paste0(images$site[1],"_","01")

for(i in 2:nrow(images)){
  # Condition! Previous images in same site AND <5s apart?
  cond <- images$deployment_id[1:(i-1)] == images$deployment_id[i] &
    (abs(difftime(images$timestamp[1:(i-1)], images$timestamp[i], units="secs")) < 5)
  
  # If condition = FALSE: it's a new sequence!
  if(!any(cond)){
    newid <- length(grep(images$site[i], unique(images$imseq.id)))+1
    
    giveID <- paste0(images$site[i],"_",sprintf("%02d", newid))
  }
  
  # If condition = TRUE: 
  if(any(cond)){
    giveID <- images[which(cond),"imseq.id"][1]
  }
  
  images[i,"imseq.id"] <- giveID
}

#* Filter only the second photo. If there is only one, get it.
seqids <- do.call(c, sapply(table(images$imseq.id), 
                  function(x) 1:x))

simages <- images[which(seqids==2),]
#simages <- images[match(unique(images$imseq.id), images$imseq.id),]

# Import Algorithm detections ----------------------------------------------

## All detected objects in EcoAssist platform
alldets <- read.csv(file.path(root,"data","detections_EcoAssist.csv"))

# Photo filename
alldets$filename <- sapply(strsplit(alldets$relative_path, "/"), function(x) tail(x,1))

# Filter AI detections only for the selected images by site and timeframe
aqui <- which(paste0(alldets$deployment_id, alldets$filename) %in%
                paste0(simages$deployment_id,simages$filename))

dets <- alldets[aqui,]

# Find the site in image data.frame and include it in dets
dets$site <- simages[match(paste0(dets$deployment_id, dets$filename),
                           paste0(simages$deployment_id,simages$filename)),
                     "site"]

# Any detection has a deployment that is not in sites?
any(!sites$deployment_id %in% dets$deployment_id)
# to check if the images were not passed through the AI

# Any AI detection has no corresponding image?
any(!paste0(dets$deployment_id,dets$filename) %in% 
      paste0(simages$deployment_id,simages$filename))

#table(dets$label,dets$site)
table(table(simages$common_name,simages$site)["White-tailed Deer",])
table(table(dets$label,dets$site)["deer",])
table(simages$site)

# Arrange data for modelling --------------------------------------------------
nsites <- nrow(sites)
ndays
max.imgs <- max(table(simages$site))

#* Create the scores table
deer <- dets[which(dets$label=="deer"),c("site","deployment_id","filename","label","confidence")]

deer$site.id <- as.integer(factor(deer$site, sites$placename))
deer$image.id <- deer$day.id <- NA

deer$confidence[which(deer$confidence==1)] <- 0.999999
deer$confidence <- qlogis(deer$confidence)

#* Create empty objects
ndets <- array(NA,dim=c(nsites,ndays,max.imgs),
                     dimnames=list(sites$placename,1:ndays,1:max.imgs))

nimgs <- matrix(NA, nrow=nsites, ncol=ndays,
                dimnames=list(sites$placename,1:ndays))

site.imgs <- split(simages, simages$site)

# Any site in images different from the deployments?
any(names(site.imgs) != sites$placename) # FALSE!


for(s in 1:nsites){
  
  # Get the Day ID for each image of this site
  days <- as.POSIXct(seq(sites$start_date[s], sites$end_date[s], 1), tz="America/New_York") + 4*3600
  site.imgs[[s]]$day.id <- findInterval(site.imgs[[s]]$timestamp,days)
  
  # Define the image ID for each day with images
  im.ids <- sapply(table(site.imgs[[s]]$day.id), 
                   function(x) 1:x)
  if(is.list(im.ids)){site.imgs[[s]]$img.id <- do.call(c, im.ids)}
  if(!is.list(im.ids)){site.imgs[[s]]$img.id <- im.ids}
  
  #* Get the number of images per day for each site
  nimgs[s,] <- table(factor(site.imgs[[s]]$day.id,levels=1:ndays))
  
  # Empty ndets
  site.imgs[[s]]$ndets <- NA
  
  # Loop over each image of each site
  for(i in 1:nrow(site.imgs[[s]])){
    img <- site.imgs[[s]][i,] # image info
    
    # find the corresponding images in the data.frame with deer AI detections
    findi <- which(paste0(deer$deployment_id, deer$filename) %in% 
                     paste0(img$deployment_id,img$filename))
    
    # if there is a deer detected
    if(length(findi)>0){
      # Fill  deer AI detections with day ID and image ID
      deer[findi,"day.id"] <- img[,"day.id"]
      deer[findi,"image.id"] <- img[,"img.id"]
    }
    # Fill 
    site.imgs[[s]][i,"ndets"] <- length(findi)
    ndets[s,img[,"day.id"],img[,"img.id"]] <- length(findi)
    
  } # i
  
} # s

site.imgs <- do.call(rbind, site.imgs)

table(ndets)


# Verified image data -----------------------------------------------------
nverif <- 200

deer$obs.class <- NA

# verifs <- sample(which(site.imgs$ndets>0),nverif)
# 
# verified <- deer[which(paste0(deer$deployment_id, deer$filename) %in%
#             paste0(site.imgs$deployment_id[verifs], site.imgs$filename[verifs])),]
# 
# write.csv(verified, file.path(root,"data","deer_200-verify.csv"),row.names = F)

#* Import verified images with labeled classes
verified <- read.csv(file.path(root,"data","deer_200-verify.csv"))

table(verified$obs.class)
ggplot(verified, aes(x=confidence, color=obs.class)) +
  geom_density(linewidth=1)

# Include observed classes in the data.frame containing all objects
aqui <-  match(paste0(deer$deployment_id, deer$filename), 
                 paste0(verified$deployment_id, verified$filename))

deer[,"obs.class"] <- ifelse(verified[aqui,"obs.class"]=="TP",1,0)

#* Create the partially observed u
# u.obs <- array(NA,dim=c(nsites,ndays,max(nimgs)),
#                dimnames=list(sites$placename,1:ndays,1:max(nimgs)))
# 
# verif.ls <- split(verified, list(verified$deployment_id, verified$filename), drop=T)
# 
# for(i in 1:length(verif.ls)){
#   s <- unique(which(rownames(u) %in% verif.ls[[i]]$site))
#   im <- unique(verif.ls[[i]]$image.id)
#   
#   u[s,im] <- ifelse(any(verif.ls[[i]]$obs.class=="TP"),1,0)
# }


# Export data -----------------------------------------------------------------------------------------------------

write.csv(site.imgs, file.path(root,"data","sites_images_annotations.csv"),row.names=F)
write.table(deer, file.path(root,"data","deer_detections_scores.txt"),row.names=F)
write.table(nimgs, file.path(root,"data","nimgs_site-visit.txt"),row.names=F)
saveRDS(ndets, file.path(root,"data","3Darray_ndets_deer.rds"))








simul_occuFPFN <- function(
    nsites = 30, # number of sites
    nvisits = 10, # number of occasions
    
    mean.nimgs = 6, # mean number of images
    
    psi = 0.6, # site occupancy prob
    phi = 0.4, # occasion availability prob
    
    lam = 1, # rate of true positives per images
    ome = 0.4, # rate of false positives per image
    
    prop.verif = 1/10, # proportion of images with verification
    
    mu.scores = c(qlogis(0.6), qlogis(.8)), # mean score for FPs and TPs
    
    sd.scr = c(0.5, 0.5) # sd for FP and TP scores
){
  
  stopifnot(mu.scores[2] > mu.scores[1] )
  
  tp.mu = mu.scores[2]
  fp.mu = mu.scores[1]
  
  # Simulate  ---------------------------------------------------------------
  
  # Occupancy state
  z <- rbinom(nsites,1,psi)
  
  # Number of images per occasion per site
  nimgs <- matrix(rpois(nsites*nvisits, mean.nimgs),nsites,nvisits)
  
  # Empty
  u <- matrix(NA, nsites, nvisits)
  K <- Q <- array(NA, dim=c(nsites, nvisits, max(nimgs)))
  
  # Availability process
  for(s in 1:nsites){
    #* True positives process
    # Species availablity in the site for a visit
    u[s,] <- rbinom(nvisits,1,z[s]*phi)
    for(j in 1:nvisits){
      # jump loop if there is no image in this visit
      if(nimgs[s,j]==0){next} 
      
      # Species presence in the image (at least one TP)
      #a[s,j,1:nimgs[s,j]] <- rbinom(nimgs[s,j],1,u[s,j]*theta)
      
      # Number of TP detected objects in the image
      K[s,j,1:nimgs[s,j]] <- rpois(nimgs[s,j],lam*u[s,j])
      # ztPoisson: E(Ni|a=1) = lambda / (1-exp(-lambda)) ; mean(Ni[a==1])
      
      #* False positives process
      Q[s,j,1:nimgs[s,j]] <- rpois(nimgs[s,j],ome) # number of FP detections in image
    } #j
  } #s
  
  
  # Simulate scores for each detection
  ndets <- K+Q # total number of detections
  ndets <- ifelse(is.na(ndets),0,ndets)
  
  
  # Object-level data -----------------------------------------------------------------------------------------------
  yscr <- list()
  
  oo=0
  for(s in 1:nsites){
    for(j in 1:nvisits){
      
      if(nimgs[s,j]==0){next}
      
      for(i in 1:nimgs[s,j]){
        oo=oo+1
        # True/False detection
        if(ndets[s,j,i]>0){
          tab <- data.frame(
            class=c(rep(1,K[s,j,i]),rep(0,Q[s,j,i])),
            score=NA,
            site=s,
            visit=j,
            image=i
          )
          # Scores for each detected object
          if(K[s,j,i]>0){
            tab[which(tab$class==1),"score"] <- rnorm(K[s,j,i], tp.mu, sd.scr[2]) # TP scores
          }
          if(Q[s,j,i]>0){
            tab[which(tab$class==0),"score"] <- rnorm(Q[s,j,i], fp.mu, sd.scr[1]) # FP scores
          }
          
          yscr[[oo]] <- tab
          
        } # ndets>0
      } #i
    } #j
  } #s
  
  obj.data <- do.call(rbind, yscr)
  
  obj.data$class.obs <- NA
  
  # Create data.frame containing the site, visit, image and ndets
  images <- as.data.frame(table(obj.data[,c("site","visit","image")]),responseName="ndets")
  
  # Exclude site-visits with no images or no detected objects
  no.obj <- which(is.na(ndets), arr.ind=T)
  images <- images[which(!paste(images$site,images$visit,images$image) %in% 
                           paste(no.obj[,1], no.obj[,2], no.obj[,3])),]
  
  
  # Sample images to verify and observe true classes ----------------------------
  
  nverif <- round(sum(images$ndets>0)*prop.verif) # number of verified images
  
  # Sample images with at least one detection to verify
  verifs <- sample((1:nrow(images))[which(images$ndets>0)],nverif)
  
  TP.obs <- FP.obs <-
    array(NA,c(nsites,nvisits,max(nimgs))) # empty partially observed variables
  
  
  for(i in 1:nverif){
    # Get the verified image
    im <- as.numeric(images[verifs[i],])
    
    
    TP.obs[im[1],im[2],im[3]] <- K[im[1],im[2],im[3]]
    FP.obs[im[1],im[2],im[3]] <- Q[im[1],im[2],im[3]]
    # Get position of the objects in the image
    pos <- which(obj.data$site==im[1] & 
                   obj.data$visit==im[2] &
                   obj.data$image==im[3])
    # Get the true class for these objects
    t.class <- obj.data[pos,"class"]
    # Fill the observed class with the true class
    obj.data[pos,"class.obs"] <- t.class
  }
  
  return(list(
    true = list(
      params=c(psi=psi, phi=phi, lam=lam, ome=ome),
      scores= rbind(mu.scores, sd.scr)
    ),
    nimgs = nimgs,
    ndets = ndets,
    obj.data = obj.data
  ))
  
} # function

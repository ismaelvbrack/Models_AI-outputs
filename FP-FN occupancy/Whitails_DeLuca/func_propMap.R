propMap <-  function(inRaster,inPoly,returnMap=FALSE){
  
  # Errors ------------------------------------------------------------------
  if(!is(inRaster,"SpatRaster")){
    stop("inRaster must be a SpatRaster")
  }
  if(!is(inPoly,"SpatVector")){
    stop("inPoly must be a SpatVector")
  }
  if(any(values(inRaster)!=1, na.rm=T)){
    stop("inRaster must contain cells with 1 values for the class intended to calculate the proportion.\n
         All other cells must equal 0 or NA")
  }
  
  
  # Calculate proportion ----------------------------------------------------
  pix.vals <- terra::extract(inRaster,inPoly,weights=T)
  
  pix.vals$var.w <- pix.vals[,2]*pix.vals[,3]
  
  prop.vals <- tapply(pix.vals[,2], pix.vals$ID, sum, na.rm=T) / 
    tapply(pix.vals[,3], pix.vals$ID, sum, na.rm=T)
  
  # Return ------------------------------------------------------------------
  
  # Only values
  if(returnMap==FALSE){
    return(prop.vals)
  }
  
  # SpatVect with values
  if(returnMap==TRUE){
    inPoly$prop.vals <- prop.vals
    return(inPoly)
  }
  
} #function
